abstract type AbstractAlgorithm end 

"""
    Forward(L::Int)

Simulated ``L`` independent lineages. When a cell divides, it is replaced by one of its daughter cells.
Does not currently support cell death. This algorithm is *not* equivalent to [`Strict`](@ref), except in
the case ``L = 1``.

**Note:** This algorithm does not allow estimation of growth rates.
"""
struct Forward <: AbstractAlgorithm
    L::Int
end 

get_thinning_rate(int, ::Forward) = 0.
function init!(alg::Forward, int)
    length(int.queue) == alg.L || throw(DimensionMismatch("Using Chemostats.Forward($(alg.L)) on a population of size $(length(int.queue))"))
end 

is_parallel(::Forward) = true
filter_offspring(cells, ::Forward, rng) = ((rand(rng, cells),), NaN)
update_algorithm!(::Forward, int) = nothing
update_queue!(queue, ::Forward, t) = nothing

### 

""" 
    Thin(rate::Float64)

Simulates a population where each cell is subject to a constant death rate ``δ \\geq 0`` (equal to `rate`). 
If the original population grows with rate ``Λ``, the simulated population grows with rate ``Λ - δ``, 
from which ``Λ`` can be estimated. 

Faster than [`Direct`](@ref), but still slow. 

**Note:** This algorithm increases the probability that the population will die out. 
To reduce the chances of this happening, increase the starting size of the population.

This algorithm requires ``δ < Λ`` to be stable, otherwise the population will eventually die out. 
This approach has exponential complexity with rate ``0 < Λ - δ < Λ``. 
"""
struct Thin <: AbstractAlgorithm
    rate::Float64 
end 

get_thinning_rate(int, alg::Thin) = alg.rate
init!(::Thin, int) = nothing
is_parallel(::Thin) = true
filter_offspring(cells, ::Thin, rng) = (cells, 0)
update_algorithm!(::Thin, int) = nothing
update_queue!(int, ::Thin, t) = nothing

""" 
    Direct()

Simulates an entire population. The standard approach. Slow. If the population grows with rate 
``Λ``, simulating it for time ``t`` will take ``e^{Λt}`` computational effort, 
while providing accuracy that only scales as ``1 / t``.

**Note:** This algorithm is implemented as an alias of `Thin(0)`.
"""
Direct() = Thin(0)

### 

"""
    Strict(L::Int)

Implements the strict chemostat, similar to the cloning algorithm of Giardinà *et al.* (2006) and Lecomte & Tailleur (2007). 
Simulates a population of exactly ``L`` cells. When a cell divides in two, one of the resulting
``L+1`` cells is discarded at random; when a cell dies, another is cloned into its place. 

Since the population size is kept constant, this algorithm has stable simulation times. Introduces an error of order 
``1 / L`` into growth rate estimates. For many problems, ``L = 10-100`` should be enough to obtain good estimates.
For ``L \\geq 100`` we recommend [`Lax`](@ref) instead.

If `clone = false`, dead cells are not replaced by clones, so the population can shrink below L. 

**Note:** This algorithm does not support parallelisation.

**Note 2:** Lineage recording currently does not support clones.
"""
struct Strict <: AbstractAlgorithm
    L::Int
    clone::Bool
end 

Strict(L::Int) = Strict(L, true)

get_thinning_rate(int, ::Strict) = 0.
function init!(alg::Strict, int)
    alg.clone && !isnothing(int.chem.tree) && @warn "Population trees are not currently supported with cloning. Use Strict(L, false) to disable cloning" maxlog=1
    update_queue!(int, alg, int.t)
end
is_parallel(alg::Strict) = false
filter_offspring(cells, ::Strict, rng) = (cells, 0)
update_algorithm!(::Strict, int) = nothing

function update_queue!(int, alg::Strict, t)
    # Strict is serial, so the heap layout is reproducible
    _shrink_pop_unsafe!((queue, n) -> _pop_random!(queue, n, int.rng), int, alg.L, t)
    alg.clone && _clone_pop_unsafe!(int, alg.L, t)
    nothing
end


###

"""
    Lax(L::Int, t_adapt::Float64; tol = 2, smoothing = 0.5, L_min = 25)

Implements the relaxed chemostat algorithm, combining [`Thin`](@ref) and [`Strict`](@ref).
The population is kept around ``L`` by introducing a time-varying death rate ``δ`` that
is adapted to the current growth rate. The death rate is updated at intervals of length `t_adapt`.

This algorithm trades the strict population size guarantees of [`Strict`](@ref) for better parallelisation.
As the population size is only approximately ``L``, runtimes may be somewhat less stable than [`Strict`](@ref).
In particular, if the population size hits ``0`` within an interval `t_adapt`, the population dies out.
For this reason, the algorithm sets ``δ = 0`` when the population drops below ``L_{min}``.
To reduce the chances of this happening, increase ``L`` or decrease `t_adapt`. We recommend ``L \\geq 50-100``.

This algorithm determines ``δ`` on the fly by estimating the instantaneous growth rate of the population as 

``\\hat Λ_n = \\frac{1 - β}{1 - β^n} \\sum_{i=0}^{n-1} β^i \\hat \\Lambda(t_{n-i-1}, t_{n-i})``

and requiring the expected population size after the next interval to equal ``L``. Here

``\\hat \\Lambda(t_{n-i-1}, t_{n-i}) = \\frac{\\log \\hat N(t_{n-i}) - \\log \\hat N(t_{n-i-1})}{\\tau}``

is the growth rate in the `n-i`th interval. The parameter `smoothing` (``β`` above, ``0 \\leq β \\leq 1``) controls smoothing:
Larger `smoothing` produces more stable estimates, but does not capture sudden variations in growth rate,
whereas smaller `smoothing` produces less stable estimates that can better adapt to fast fluctuations.
If `smoothing = 0`, only the growth rate in the last interval is used. The default is `smoothing = 0.5`.
"""
struct Lax <: AbstractAlgorithm
    L::Int
    t_adapt::Float64
    tol::Float64
    smoothing::Float64
    L_min::Int

    function Lax(L::Int, t_adapt; tol=2., smoothing=0.5, L_min=25)
        @argcheck tol >= 1
        @argcheck t_adapt > 0
        @argcheck 0 <= smoothing <= 1
        @argcheck L_min >= 1

        new(L, t_adapt, tol, smoothing, L_min)
    end
end

# Snapshot closest to `t` if within `t_adapt/2`, searching back from index `i`. Also returns
# the index to continue from, so decreasing times take a single pass
function get_snapshot_near(snaps::AbstractVector{Snapshot}, t, t_adapt, i = lastindex(snaps))
    # Last snapshot before `t`, then the closer of it and the next one
    while i > firstindex(snaps) && snaps[i].t > t
        i -= 1
    end
    j = i < lastindex(snaps) && snaps[i + 1].t - t < t - snaps[i].t ? i + 1 : i
    (abs(snaps[j].t - t) <= t_adapt / 2 ? snaps[j] : nothing), i
end

function est_Λ_curr(snaps::AbstractVector{Snapshot}, t, alg::Lax)
    rnd = floor(Int, (t - get_start_t(snaps)) / alg.t_adapt + 1e-6)
    rnd <= 0 && return NaN

    # Snapshots at t - k * t_adapt for k = 1, …, rnd; a missing one makes the estimate undefined
    snap_2 = _get_snapshot(snaps, t)
    i = lastindex(snaps)
    ret = 0.0
    for d in 0:rnd-1
        snap_1, i = get_snapshot_near(snaps, t - (d + 1) * alg.t_adapt, alg.t_adapt, i)
        isnothing(snap_1) && return NaN
        ret += (d > 0 ? alg.smoothing^d : one(alg.smoothing)) * est_growth_rate(snap_1, snap_2)
        snap_2 = snap_1
    end

    norm = alg.smoothing < 1 ? ((1 - alg.smoothing^rnd) / (1 - alg.smoothing)) : rnd * alg.smoothing
    ret / norm
end

init!(alg::Lax, int) = update_algorithm!(alg, int)

get_thinning_rate(int, alg::Lax) = get_thinning_rate(int.chem, int.t, alg)
get_thinning_rate(chem::Chemostat, t, alg::Lax) = get_thinning_rate(chem.snaps, t, alg)

function get_thinning_rate(snaps::AbstractVector{Snapshot}, t, alg::Lax)
    snap = _get_snapshot(snaps, t)

    N = snap.N > alg.L * alg.tol ? alg.L : snap.N

    # Small population
    N < alg.L_min && return 0.
    t - get_start_t(snaps) < alg.t_adapt && return 0.

    # Estimate current growth rate
    Λ̂ = est_Λ_curr(snaps, t, alg)
    isfinite(Λ̂) || return 0.

    # Expected population size after time t_adapt is
    #   `length(queue) * exp((Λ - δ) * t_adapt)`
    # Choose δ to make this equal to L
    ΔlogN = log(alg.L) - log(N)
    ret = Λ̂ - ΔlogN / alg.t_adapt
    max(0, ret)
end

function update_algorithm!(alg::Lax, int)
    if length(int.queue) > alg.L * alg.tol
        # After a multithreaded step, the heap layout depends on scheduling
        _shrink_pop_unsafe!(_pop_by_keys!, int, alg.L, int.t)
    end

    add_tstop!(int, int.t + alg.t_adapt)
end 

is_parallel(::Lax) = true
filter_offspring(cells, ::Lax, rng) = (cells, 0)
update_queue!(int, ::Lax, t) = nothing

###

# Population control, used by `Strict` and `Lax`

# `pop_random!(queue, n)` removes `n` cells chosen at random
function _shrink_pop_unsafe!(pop_random!, int, L::Int, t)
    N_start = length(int.queue)
    N_start > L || return

    for cell in pop_random!(int.queue, N_start - L)
        data = get_data(cell)
        kill!(cell, t)
        @lock int.tree_lock record_end!(int.chem.tree, cell, t, data)
    end

    int.log_f += log(N_start) - log(length(int.queue))
end

function _clone_pop_unsafe!(int, L::Int, t)
    isempty(int.queue) && return 

    N_start = length(int.queue)
    N_start < L || return
    # Trees cannot represent clones
    isnothing(int.chem.tree) || throw(ArgumentError(
        "Cannot record a population tree when cloning cells, use `Strict(L, false)`, `Lax` or `save_mode=:none`"))

    while length(int.queue) < L
        _clone_random!(int.queue, t, int.rng)
    end

    int.log_f += log(N_start) - log(length(int.queue))
end
