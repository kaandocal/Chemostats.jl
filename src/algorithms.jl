abstract type AbstractAlgorithm end 

struct Forward <: AbstractAlgorithm
    L::Int
end 

get_δ(int, ::Forward) = 0.
init!(alg::Forward, int) = resize_pop!(int, alg.L)
is_parallel(::Forward) = true
update_algorithm!(::Forward, int) = nothing
filter_offspring(cells, ::Forward) = Iterators.take(cells, 1)
should_sync(queue, ::Forward) = false
sync!(queue, ::Forward) = nothing

### 

struct Thin <: AbstractAlgorithm
    δ::Float64 
end 

get_δ(int, alg::Thin) = alg.δ
init!(::Thin, int) = nothing
is_parallel(::Thin) = true
update_algorithm!(::Thin, int) = nothing
filter_offspring(cells, ::Thin) = cells 
should_sync(queue, alg::Thin) = false
sync!(int, ::Thin) = nothing

Direct() = Thin(0)

### 

struct Strict <: AbstractAlgorithm
    L::Int
end 

get_δ(int, ::Strict) = 0.
init!(alg::Strict, int) = resize_pop!(int, alg.L)
is_parallel(alg::Strict) = false
update_algorithm!(::Strict, int) = nothing
filter_offspring(cells, ::Strict) = cells 
should_sync(queue, ::Strict) = true
sync!(int, alg::Strict) = resize_pop!(int, alg.L)

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
    snap.N < 10 && return 0.

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

update_algorithm!(::Lax, int) = nothing  
is_parallel(::Lax) = true
filter_offspring(cells, ::Lax) = cells
should_sync(queue, alg::Lax) = length(queue) > alg.L * (1 + alg.tol)
sync!(int, alg::Lax) = resize_pop!(int, alg.L)
