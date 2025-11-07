# Adapted from DataStructures.jl 
_getindex(queue::BinaryHeap, i::Integer) = queue.valtree[i]

function _force_up!(queue::BinaryHeap, i::Integer)
    x = queue.valtree[i]

    @inbounds while i > 1
        j = DataStructures.heapparent(i)
        queue.valtree[i] = queue.valtree[j]
        i = j
    end

    queue.valtree[i] = x
end 

function _popat!(queue::BinaryHeap, idx::Integer)
    _force_up!(queue, idx)
    pop!(queue)
end 

function duplicate_random!(queue::BinaryHeap, t)
    i = sample(1:length(queue))
    push!(queue, duplicate_cell(queue.valtree[i], t))
end 

function truncate_queue!(queue, N)
    while length(queue) > N
        j = sample(1:length(queue))
        _popat!(queue, j)
    end 
end 

function ensure_size!(queue, L::Int, t)
    while length(queue) < L
        duplicate_random!(queue, t)
    end

    truncate_queue!(queue, L)
end 

### 

abstract type AbstractAlgorithm end 

struct Forward <: AbstractAlgorithm
    L::Int
end 

get_δ(int, ::Forward) = 0.
init!(alg::Forward, int) = ensure_size!(int.queue, alg.L, int.t)
update_algorithm!(::Forward, int) = nothing

function add_offspring!(queue, cells, ::Forward)
    isempty(cells) && return 0
    push!(queue, first(cells))
    return 1
end

### 

struct Thin <: AbstractAlgorithm
    δ::Float64 
end 

get_δ(int, alg::Thin) = alg.δ
init!(::Thin, int) = nothing
update_algorithm!(::Thin, int) = nothing

function add_offspring!(queue, cells, ::Thin)
    for cell in cells
        push!(queue, cell)
    end

    return length(cells)
end

### 

struct Strict <: AbstractAlgorithm
    L::Int
end 

get_δ(int, ::Strict) = 0.
init!(alg::Strict, int) = resize_queue!(int, alg)
update_algorithm!(alg::Strict, int) = resize_queue!(int, alg)

function resize_queue!(int, alg::Strict) 
    if isempty(int.queue)
        @warn "No cells left in strict chemostat"
        int.retcode = ReturnCode.Unstable
        return
    end 

    # Ensure we have exactly L clones 
    N_start = length(int.queue)
    ensure_size!(int.queue, alg.L, int.t)
    int.log_f += log(N_start) - log(alg.L)
end 

function add_offspring!(queue, cells, alg::Strict)
    @check length(queue) == alg.L - 1

    nsim = 0

    for cell in cells
        push!(queue, cell)
        nsim += 1
    end 

    while length(queue) > alg.L
        j = sample(1:length(queue))
        if _getindex(queue, j) in cells
            nsim -= 1
        end 

        _popat!(queue, j)
    end 

    return nsim
end

###

struct Lax <: AbstractAlgorithm
    L::Int
    t_adapt::Float64
    tol::Float64

    function Lax(L::Int, t_adapt, tol=0.5)
        @argcheck 0 < tol < 1
        @argcheck t_adapt > 0

        new(L,t_adapt, tol)
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
    # This should == L
    ΔlogN = log(alg.L) - log(snap.N)
    ret = Λ̂ - ΔlogN / alg.t_adapt
    max(0, ret)
end 

update_algorithm!(alg::Lax, int) = nothing  

function add_offspring!(queue, cells, alg::Lax)
    for cell in cells
        push!(queue, cell)
    end

    if length(queue) > alg.L * (1 + alg.tol)
        truncate_queue!(queue, alg.L)
    end

    return length(cells)
end
