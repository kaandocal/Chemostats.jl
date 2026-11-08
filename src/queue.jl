function get_curr_t end
const TimeOrder = Base.By(get_curr_t)

"""
    ThreadedQueue{H}

A thread-safe priority queue (min-heap ordered by `get_curr_t`). Knows nothing
about "workers" or when work is exhausted -- see [`WorkerPool`](@ref) for that.
"""
mutable struct ThreadedQueue{H <: MutableBinaryHeap}
    heap::H
    lock::ReentrantLock
    cond_wait::Threads.Condition

    function ThreadedQueue(heap::MutableBinaryHeap)
        lock = ReentrantLock()
        new{typeof(heap)}(heap, lock, Threads.Condition(lock))
    end
end

ThreadedQueue(vals) = ThreadedQueue(MutableBinaryHeap(TimeOrder, vals))

Base.lock(f::Function, queue::ThreadedQueue) = lock(f, queue.lock)
Base.length(queue::ThreadedQueue) = @lock queue.lock length(queue.heap)
Base.isempty(queue::ThreadedQueue) = length(queue) == 0
Base.first(queue::ThreadedQueue) = first(queue.heap)

# Must be called while already holding `queue.lock` (i.e. from within `lock(queue) do ... end`).
Base.wait(queue::ThreadedQueue) = wait(queue.cond_wait)
Base.notify(queue::ThreadedQueue; all=true) = notify(queue.cond_wait; all)

function Base.push!(queue::ThreadedQueue, v)
    @lock queue.lock begin
        push!(queue.heap, v)
        notify(queue.cond_wait; all=false)
    end
end

"""
    trypop!(queue::ThreadedQueue)

Pops the minimum element, or returns `nothing` immediately if empty. Never blocks.
"""
function trypop!(queue::ThreadedQueue)
    @lock queue.lock begin
        isempty(queue.heap) ? nothing : pop!(queue.heap)
    end
end

"""
    iter_unsafe(queue::ThreadedQueue)

Lazy iterable over `queue`'s current contents -- no copy, no wrapper allocation.
NOT thread safe: no locking. Only valid where nothing else can be concurrently
mutating `queue` (e.g. after `@sync` in `step!` has already rejoined every
worker). A lock per element wouldn't give a true snapshot anyway (another task
could still mutate between steps), and holding the lock across the whole
iteration would mean holding it across arbitrary caller code in between.
"""
iter_unsafe(queue::ThreadedQueue) = Iterators.map(node -> node.value, queue.heap.nodes)

function _append!(queue::ThreadedQueue, vals)
    @lock queue.lock begin
        for v in vals
            push!(queue.heap, v)
        end
        notify(queue.cond_wait; all=true)
    end
end

function extract_queue!(pop, queue::ThreadedQueue)
    @lock queue.lock begin
        while !isempty(queue.heap)
            push!(pop, pop!(queue.heap))
        end
    end
end

# Removes the values at positions `idxs` of the heap's nodes
function _pop_nodes_unsafe!(heap::MutableBinaryHeap, idxs)
    nodes = heap.nodes[idxs]
    for node in nodes
        delete!(heap, node.handle)
    end
    [ node.value for node in nodes ]
end

# Removes `n` values chosen uniformly at random
function _pop_random!(queue::ThreadedQueue, n::Int, rng)
    @lock queue.lock _pop_nodes_unsafe!(queue.heap, sample(rng, eachindex(queue.heap.nodes), n; replace=false))
end

# As `_pop_random!`, but independent of the heap layout, which varies under multithreading:
# each value draws a key from its own RNG, and the `n` smallest keys are removed
function _pop_by_keys!(queue::ThreadedQueue, n::Int)
    @lock queue.lock begin
        keys = [ rand(get_cell_rng(node.value)) for node in queue.heap.nodes ]
        _pop_nodes_unsafe!(queue.heap, partialsortperm(keys, 1:n))
    end
end

function _clone_random!(heap::MutableBinaryHeap, t, rng)
    source = rand(rng, heap.nodes).value
    cell = clone_cell(source, t)
    push!(heap, cell)
    source, cell
end

_clone_random!(queue::ThreadedQueue, t, rng) = @lock queue.lock _clone_random!(queue.heap, t, rng)
