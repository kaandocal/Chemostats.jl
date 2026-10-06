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

function _pop_random_unsafe!(heap::MutableBinaryHeap)
    i = rand(1:length(heap.nodes))
    v = heap.nodes[i].value
    delete!(heap, heap.nodes[i].handle)
    v
end

_pop_random!(queue::ThreadedQueue) = @lock queue.lock _pop_random_unsafe!(queue.heap)

function _clone_random!(heap::MutableBinaryHeap, t)
    i = rand(1:length(heap.nodes))
    source = heap.nodes[i].value
    cell = clone_cell(source, t)
    push!(heap, cell)
    source, cell
end

_clone_random!(queue::ThreadedQueue, t) = @lock queue.lock _clone_random!(queue.heap, t)
