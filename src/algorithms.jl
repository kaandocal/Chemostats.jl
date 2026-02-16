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

get_δ(int, ::Forward) = 0.
function init!(alg::Forward, int)
    length(int.queue) == alg.L || throw(DimensionMismatch("Using Chemostats.Forward($(alg.L)) on a population of size $(length(int.queue))"))
end 

is_parallel(::Forward) = true
filter_offspring(cells, ::Forward) = (Iterators.take(cells, 1), NaN)
update_algorithm!(::Forward, int) = nothing
update_queue!(queue, ::Forward, t) = nothing

### 

""" 
    Thin(δ::Float64)

Simulates a population where each cell is subject to a constant death rate ``δ > 0``. 
If the original population grows with rate ``Λ``, the simulated population grows with rate ``Λ - δ``, 
from which ``Λ`` can be estimated. 

Faster than [`Direct`](@ref), but still slow. 

**Note:** This algorithm increases the chances that the population will die out. 
To reduce the chances of this happening, increase the starting size of the population.

This algorithm requires ``δ < Λ`` to be stable, otherwise the population will eventually die out. 
This approach has exponential complexity with rate ``0 < Λ - δ < Λ``. 
"""
struct Thin <: AbstractAlgorithm
    δ::Float64 
end 

get_δ(int, alg::Thin) = alg.δ
init!(::Thin, int) = nothing
is_parallel(::Thin) = true
filter_offspring(cells, ::Thin) = (cells, 0)
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

Implements the cloning algorithm of Giardinà *et al.* (2006) and Lecomte & Tailleur (2007). 
Simulates a population of exactly ``L`` cells. When a cell divides in two, one of the resulting
``L+1`` cells is discarded at random; when a cell dies, another is cloned into its place. 

Since the population size is kept constant, it has stable simulation times. Introduces an error of order 
``1 / L`` into growth rate estimates. For many problems, ``L = 10-100`` should be enough to obtain good estimates.
For ``L \\geq 100`` we recommend [`Lax`](@ref) instead.

**Note:** This algorithm does not support parallelisation.
"""
struct Strict <: AbstractAlgorithm
    L::Int
end 

get_δ(int, ::Strict) = 0.
init!(alg::Strict, int) = _resize_pop!(int, alg.L, int.t)
is_parallel(alg::Strict) = false
filter_offspring(cells, ::Strict) = (cells, 0)
update_algorithm!(::Strict, int) = nothing
update_queue!(int, alg::Strict, t) = _resize_pop!(int, alg.L, t)

###

"""
    Lax(L::Int, τ::Float64; tol = 2, β = 0.5)

Implements the relaxed chemostat algorithm, combining [`Thin`](@ref) and [`Strict`](@ref). 
The population is kept around ``L`` by introducing a time-varying death rate ``δ`` that 
is adapted to the current growth rate. The death rate is updated at intervals of length ``τ``. 

This algorithm trades the strict population size guarantees of [`Strict`](@ref) for better parallelisation. 
As the population size is only approximately ``L``, runtimes may be somewhat less stable than [`Strict`](@ref). 
In particular, if the population size hits ``0`` within an interval ``τ``, the population dies out. 
To reduce the chances of this happening, increase ``L`` or decrease ``τ``. We recommend ``L \\geq 50-100``.

This algorithm determines ``δ`` on the fly by estimating the instantaneous growth rate of the population as 

``\\hat Λ_n = \\sum_{i=0}^\\infty \\frac {β^i} {1 - β}  \\hat \\Lambda(t_{n-i-1}, t_{n-i})``

and requiring the expected population size after the next interval to equal ``L``. Here

``\\hat \\Lambda(t_{n-i-1}, t_{n-i}) = \\frac{\\log \\hat N(t_{n-i}) - \\log \\hat N(t_{n-i-1})}{\\tau}``

is the growth rate in the `n-i`th interval. The parameter ``0 \\leq β \\leq 1`` controls smoothing:
Larger ``β`` produces more stable estimates, but does not capture sudden variations in growth rate,
whereas smaller ``β`` produces less stable estimates that can better adapt to fast fluctuations.
If ``β = 0``, only the growth rate in the last interval is used.
"""
struct Lax <: AbstractAlgorithm
    L::Int
    τ::Float64
    tol::Float64
    β::Float64

    function Lax(L::Int, τ; tol=2., β=0.5)
        @argcheck tol >= 1
        @argcheck τ > 0
        @argcheck 0 <= β <= 1

        new(L, τ, tol, β)
    end 
end 

function est_Λ_curr(snaps::AbstractVector{Snapshot}, t, alg::Lax)
    rnd = floor(Int, t / alg.τ + 1e-6)
    rnd <= 0 && return NaN

    ret = sum(0:rnd-1) do d 
        snap_1 = get_snapshot(snaps, t - (d + 1) * alg.τ)
        snap_2 = get_snapshot(snaps, t - d * alg.τ)
        Λ_est = est_Λ(snap_1, snap_2)
        (d > 0 ? alg.β^d : one(alg.β)) * Λ_est
    end 

    norm = alg.β < 1 ? ((1 - alg.β^rnd) / (1 - alg.β)) : rnd * alg.β
    ret / norm
end 

init!(alg::Lax, int) = update_algorithm!(alg, int)

get_δ(int, alg::Lax) = get_δ(int.chem, int.t, alg)
get_δ(chem::Chemostat, t, alg::Lax) = get_δ(chem.snaps, t, alg)

function get_δ(snaps::AbstractVector{Snapshot}, t, alg::Lax)
    # Small population 
    snap = get_snapshot(snaps, t)
    snap.N < 50 && return 0.
    t < alg.τ && return 0.

    # Estimate current growth rate
    Λ̂ = est_Λ_curr(snaps, t, alg)
    isfinite(Λ̂) || return 0.

    # Expected population size after time τ is 
    #   `length(queue) * exp((Λ - δ) * τ)` 
    # Choose δ to make this equal to L
    ΔlogN = log(alg.L) - log(snap.N)
    ret = Λ̂ - ΔlogN / alg.τ
    max(0, ret)
end 

function update_algorithm!(alg::Lax, int)
    if length(int.queue) > alg.L * alg.tol
        _resize_pop!(int, alg.L, int.t)
    end 

    add_tstop!(int, int.t + alg.τ)
end 

is_parallel(::Lax) = true
filter_offspring(cells, ::Lax) = (cells, 0)
update_queue!(int, ::Lax, t) = nothing