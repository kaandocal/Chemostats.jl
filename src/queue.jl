lockqueue(f::Function, ::BinaryHeap) = f()

###

struct ThreadedBinaryHeap{HT <: BinaryHeap}
    heap::HT 
    lock::ReentrantLock
end

ThreadedBinaryHeap(heap::BinaryHeap) = ThreadedBinaryHeap(heap, ReentrantLock())
Base.eltype(heap::ThreadedBinaryHeap) = eltype(heap.heap)
lockqueue(f::Function, queue::ThreadedBinaryHeap) = lock(f, queue.lock)

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

const QueueType = Union{BinaryHeap,ThreadedBinaryHeap}

function get_t end;
const TimeOrder = Base.By(get_t)

create_queue(vals, ::EnsembleSerial) = BinaryHeap(TimeOrder, vals)
create_queue(vals, ::EnsembleThreads) = ThreadedBinaryHeap(BinaryHeap(TOrder, vals))

function Base.append!(queue::QueueType, vals)
    lockqueue(queue) do 
        for v in vals 
            push!(queue, v)
        end 
    end
end

function extract_queue!(f::Function, pop, queue::QueueType)
    lockqueue(queue) do 
        while !isempty(queue)
            push!(pop, f(pop!(queue)))
        end
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

function _truncate_queue!(queue::BinaryHeap, N)
    while length(queue) > N
        j = sample(1:length(queue))
        _popat!(queue, j)
    end 
end 
