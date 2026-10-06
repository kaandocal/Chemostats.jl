function get_curr_t end;
const TimeOrder = Base.By(get_curr_t)

mutable struct ThreadedQueue{H <: BinaryHeap}
    heap::H
    lock::ReentrantLock
    cond_wait::Threads.Condition
    @atomic nwork::Int

    function ThreadedQueue(heap::BinaryHeap)
        lock = ReentrantLock()
        new{typeof(heap)}(heap, lock, Threads.Condition(lock), 0)
    end
end 

Base.lock(f::Function, queue::ThreadedQueue) = lock(f, queue.lock)
function Base.length(queue::ThreadedQueue) 
    @lock queue.lock length(queue.heap) 
end

Base.isempty(queue::ThreadedQueue) = length(queue) == 0
ThreadedQueue(vals) = ThreadedQueue(BinaryHeap(TimeOrder, vals))

function register_listener!(queue::ThreadedQueue)
    # Race condition?
    @atomic queue.nwork += 1
    @debug "Thread $(Threads.threadid()): register (# $(queue.nwork))..."
end

function Base.push!(queue::ThreadedQueue, v)
    @lock queue.lock begin 
        @debug "Thread $(Threads.threadid()): put..."
        push!(queue.heap, v)
        notify(queue.cond_wait; all=false)
    end
end 

function fetch!(queue::ThreadedQueue)
    @debug "Thread $(Threads.threadid()): fetching..."
    @lock queue.cond_wait begin
        @atomic queue.nwork -= 1
        while isempty(queue.heap)
            if queue.nwork == 0
                @debug "Thread $(Threads.threadid()): detecting done..."
                notify(queue.cond_wait; all=true)
                return nothing
            end 

            @debug "Thread $(Threads.threadid()): wait..."
            wait(queue.cond_wait)
            @debug "Thread $(Threads.threadid()): wake..."
        end

        # Two different locks here
        @debug "Thread $(Threads.threadid()): take ($(queue.nwork) waiting)..."
        @atomic queue.nwork += 1
        ret = pop!(queue.heap)
        notify(queue.cond_wait; all=false)
        ret
    end
end 

###

function release_and_notify!(queue::ThreadedQueue)
    @lock queue.cond_wait begin
        @atomic queue.nwork -= 1
        notify(queue.cond_wait; all=true)
    end
end

Base.first(queue::ThreadedQueue) = first(queue.heap)

"""
    iter_unsafe(queue::ThreadedQueue)

Returns `queue`'s backing array directly -- already a fully-featured
iterable, so no wrapper type is needed. NOT thread safe: no locking, no
copy. Only valid where nothing else can be concurrently mutating `queue`
(e.g. after `@sync` in `step!` has already rejoined every worker). A lock
per element wouldn't give a true snapshot anyway (another task could still
mutate between steps), and holding the lock across the whole iteration would
mean holding it across arbitrary caller code in between.
"""
iter_unsafe(queue::ThreadedQueue) = queue.heap.valtree

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
        while !isempty(queue)
            push!(pop, pop!(queue.heap))
        end
    end
end 

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

function _popat!(queue::ThreadedQueue, idx::Integer)
    @lock queue.lock begin
        _popat!(queue.heap, idx)
    end
end 

function _clone_random!(queue::BinaryHeap, t)
    i = rand(1:length(queue))
    source = queue.valtree[i]
    cell = clone_cell(source, t)
    push!(queue, cell)
    source, cell
end

_clone_random!(queue::ThreadedQueue, t) = _clone_random!(queue.heap, t)
