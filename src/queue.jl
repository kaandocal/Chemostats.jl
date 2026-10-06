function get_curr_t end;
const TimeOrder = Base.By(get_curr_t)

mutable struct ThreadedQueue{H <: MutableBinaryHeap}
    heap::H
    lock::ReentrantLock
    cond_wait::Threads.Condition
    @atomic nwork::Int

    function ThreadedQueue(heap::MutableBinaryHeap)
        lock = ReentrantLock()
        new{typeof(heap)}(heap, lock, Threads.Condition(lock), 0)
    end
end

Base.lock(f::Function, queue::ThreadedQueue) = lock(f, queue.lock)
function Base.length(queue::ThreadedQueue)
    @lock queue.lock length(queue.heap)
end

Base.isempty(queue::ThreadedQueue) = length(queue) == 0
ThreadedQueue(vals) = ThreadedQueue(MutableBinaryHeap(TimeOrder, vals))

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
        while !isempty(queue)
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

function _pop_random!(queue::ThreadedQueue)
    @lock queue.lock begin
        _pop_random_unsafe!(queue.heap)
    end
end

function _clone_random!(heap::MutableBinaryHeap, t)
    i = rand(1:length(heap.nodes))
    source = heap.nodes[i].value
    cell = clone_cell(source, t)
    push!(heap, cell)
    source, cell
end

_clone_random!(queue::ThreadedQueue, t) = _clone_random!(queue.heap, t)
