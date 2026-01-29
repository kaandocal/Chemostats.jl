abstract type AbstractAlgorithm end 

### Do not get proper growth rate estimates!!!
struct Forward <: AbstractAlgorithm
    L::Int
end 

get_δ(int, ::Forward) = 0.
init!(alg::Forward, int) = _resize_pop!(int, alg.L, int.t)
is_parallel(::Forward) = true
filter_offspring(cells, ::Forward) = (Iterators.take(cells, 1), NaN)
update_algorithm!(::Forward, int) = nothing
update_queue!(queue, ::Forward, t) = nothing

### 

struct Thin <: AbstractAlgorithm
    δ::Float64 
end 

get_δ(int, alg::Thin) = alg.δ
init!(::Thin, int) = nothing
is_parallel(::Thin) = true
filter_offspring(cells, ::Thin) = (cells, 0)
update_algorithm!(::Thin, int) = nothing
update_queue!(int, ::Thin, t) = nothing

Direct() = Thin(0)

### 

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
    rnd = get_round(t, alg)
    rnd == 0 && return NaN

    ret = sum(0:rnd-1) do d 
        snap_1 = get_snapshot(snaps, (rnd - d - 1) * alg.τ)
        snap_2 = get_snapshot(snaps, (rnd - d) * alg.τ)
        Λ_est = est_Λ(snap_1, snap_2)
        (d > 0 ? alg.β^d : one(alg.β)) * Λ_est
    end 

    norm = alg.β < 1 ? ((1 - alg.β^rnd) / (1 - alg.β)) : rnd * alg.β
    ret / norm
end 

get_round(t, alg::Lax) = floor(Int, div(t, alg.τ))

function register_next_tstop!(int, alg::Lax) 
    rnd = get_round(int.t, alg)
    add_tstop!(int, (rnd + 1) * alg.τ)
end 

init!(alg::Lax, int) = register_next_tstop!(int, alg)

function get_δ(int, alg::Lax)
    get_δ(int.chem, int.t, alg)
end

function get_last_snapshot(snaps::AbstractVector{Snapshot}, t, alg::Lax)
    rnd = get_round(t, alg)
    get_snapshot(snaps, rnd * alg.τ)
end 

get_δ(chem::Chemostat, t, alg::Lax) = get_δ(chem.saved, t, alg)

function get_δ(snaps::AbstractVector{Snapshot}, t, alg::Lax)
    snap = get_last_snapshot(snaps, t, alg)

    # Small population 
    snap.N < 50 && return 0.

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

    register_next_tstop!(int, alg)
end 

is_parallel(::Lax) = true
filter_offspring(cells, ::Lax) = (cells, 0)
update_queue!(int, ::Lax, t) = nothing