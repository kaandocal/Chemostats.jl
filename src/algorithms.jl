abstract type AbstractAlgorithm end 

struct Forward <: AbstractAlgorithm
    L::Int
end 

get_δ(int, ::Forward) = 0.
init!(alg::Forward, int) = resize_pop!(int, alg.L, int.t)
is_parallel(::Forward) = true
filter_offspring(cells, ::Forward) = Iterators.take(cells, 1)
update_algorithm!(::Forward, int) = nothing
update_queue!(queue, ::Forward, t) = nothing

### 

struct Thin <: AbstractAlgorithm
    δ::Float64 
end 

get_δ(int, alg::Thin) = alg.δ
init!(::Thin, int) = nothing
is_parallel(::Thin) = true
filter_offspring(cells, ::Thin) = cells 
update_algorithm!(::Thin, int) = nothing
update_queue!(int, ::Thin, t) = nothing

Direct() = Thin(0)

### 

struct Strict <: AbstractAlgorithm
    L::Int
end 

get_δ(int, ::Strict) = 0.
init!(alg::Strict, int) = resize_pop!(int, alg.L, int.t)
is_parallel(alg::Strict) = false
filter_offspring(cells, ::Strict) = cells 
update_algorithm!(::Strict, int) = nothing
update_queue!(int, alg::Strict, t) = resize_pop!(int, alg.L, t)

###

struct Lax <: AbstractAlgorithm
    L::Int
    t_adapt::Float64
    tol::Float64

    function Lax(L::Int, t_adapt, tol=0.5)
        @argcheck 0 < tol < 1
        @argcheck t_adapt > 0

        new(L, t_adapt, tol)
    end 
end 

function est_Λ_curr(chem::Chemostat, t, alg::Lax)
    snap1 = get_snapshot(chem, t)

    @assert issorted(chem.saved; by = snap -> snap.t)
    i = searchsortedlast(map(snap -> snap.t, chem.saved), t - alg.t_adapt)
    isempty(i) && return NaN
    snap0 = chem.saved[i]

    est_Λ(snap0, snap1)
end 

function register_tstop!(int, alg::Lax) 
    int.t_next < int.t + alg.t_adapt && add_tstop!(int, int.t + alg.t_adapt)
end 

init!(alg::Lax, int) = register_tstop!(int, alg)

function get_δ(int, alg::Lax)
    register_tstop!(int, alg)

    get_δ(int.chem, int.t, alg)
end

function get_δ(chem::Chemostat, t, alg::Lax)
    snap = get_snapshot(chem, t)

    # Small population 
    snap.N < 50 && return 0.

    # Estimate current growth rate
    Λ̂ = est_Λ_curr(chem, t, alg)
    isfinite(Λ̂) || return 0.

    # Expected population size after time t_adapt is 
    #   `length(queue) * exp((Λ - δ) * t_adapt)` 
    # Choose δ to make this equal to L
    ΔlogN = log(alg.L) - log(snap.N)
    ret = Λ̂ - ΔlogN / alg.t_adapt
    max(0, ret)
end 

function update_algorithm!(alg::Lax, int)
    if length(int.queue) > alg.L * (1 + alg.tol)
        resize_pop!(int, alg.L, int.t)
    end 
end 

is_parallel(::Lax) = true
filter_offspring(cells, ::Lax) = cells
update_queue!(int, ::Lax, t) = nothing