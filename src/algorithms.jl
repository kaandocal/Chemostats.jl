### Queue management

struct ThreadedBinaryHeap{HT <: BinaryHeap}
    heap::HT 
    lock::ReentrantLock
end

ThreadedBinaryHeap(heap::BinaryHeap) = ThreadedBinaryHeap(heap, ReentrantLock())
Base.eltype(heap::ThreadedBinaryHeap) = eltype(heap.heap)

lockqueue(f::Function, ::BinaryHeap) = f()
lockqueue(f::Function, queue::ThreadedBinaryHeap) = lock(f, queue.lock)

const QueueType = Union{BinaryHeap,ThreadedBinaryHeap}

function Base.first(queue::ThreadedBinaryHeap) 
    lockqueue(queue) do 
        first(queue.heap)
    end 
end

function Base.popfirst!(queue::ThreadedBinaryHeap) 
    lockqueue(queue) do 
        popfirst!(queue.heap)
    end 
end

function Base.push!(queue::ThreadedBinaryHeap, v)
    lockqueue(queue) do 
        push!(queue.heap, v)
    end 
end 

function Base.append!(queue::QueueType, vals)
    lockqueue(queue) do 
        for v in vals 
            push!(queue, v)
        end 
    end
end

function resize_int!(int, L::Int)
    lockqueue(int.queue) do 
        if isempty(int.queue)
            @warn "No cells left in chemostat, terminating..."
            int.retcode = ReturnCode.Unstable
            return
        end  

        N_start = length(int.queue)
        
        while length(int.queue) < L
            _duplicate_random!(int.queue, int.t)
        end

        _truncate_queue!(int.queue, L)

        int.log_f += log(N_start) - log(L)
    end
end

# Adapted from DataStructures.jl 
_getindex(queue::BinaryHeap, i::Integer) = queue.valtree[i]

# THESE FUNCTIONS ARE NOT THREAD SAFE 
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

function _duplicate_random!(queue::BinaryHeap, t)
    i = sample(1:length(queue))
    push!(queue, duplicate_cell(queue.valtree[i], t))
end 

function _truncate_queue!(queue, N)
    while length(queue) > N
        j = sample(1:length(queue))
        _popat!(queue, j)
    end 
end 

### 

abstract type AbstractAlgorithm end 

struct Forward <: AbstractAlgorithm
    L::Int
end 

get_δ(int, ::Forward) = 0.
init!(alg::Forward, int) = resize_int!(int, alg.L)
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

### 

struct Strict <: AbstractAlgorithm
    L::Int
end 

get_δ(int, ::Strict) = 0.
init!(alg::Strict, int) = resize_int!(int, alg.L)
is_parallel(alg::Strict) = false
update_algorithm!(::Strict, int) = nothing
filter_offspring(cells, ::Strict) = cells 
should_sync(queue, ::Strict) = true
sync!(int, alg::Strict) = resize_int!(int, alg.L)

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
    # This should == L
    ΔlogN = log(alg.L) - log(snap.N)
    ret = Λ̂ - ΔlogN / alg.t_adapt
    max(0, ret)
end 

update_algorithm!(::Lax, int) = nothing  
is_parallel(::Lax) = true
filter_offspring(cells, ::Lax) = cells
should_sync(queue, alg::Lax) = length(queue) > alg.L * (1 + alg.tol)
sync!(int, alg::Lax) = resize_int!(int, alg.L)
