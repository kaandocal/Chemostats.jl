"""
    WorkerPool(queue)

Tracks how many workers sharing `queue` are currently active, so [`fetch!`](@ref)
can tell "temporarily empty, another worker may still push more" apart from
"empty and nobody left who could ever push again".
"""
mutable struct WorkerPool{Q}
    queue::Q
    @atomic active::Int
end

WorkerPool(queue) = WorkerPool(queue, 0)

function register!(wp::WorkerPool)
    @atomic wp.active += 1
    nothing
end

"""
    fetch!(wp::WorkerPool)

Pops the next item, waiting if `wp`'s queue is momentarily empty but some other
worker is still active. Returns `nothing` once the queue is empty *and* no
worker is active, i.e. nothing more can ever arrive.
"""
function fetch!(wp::WorkerPool)
    queue = wp.queue
    lock(queue) do
        @atomic wp.active -= 1
        while isempty(queue)
            if wp.active == 0
                notify(queue; all=true)
                return nothing
            end
            wait(queue)
        end
        @atomic wp.active += 1
        v = trypop!(queue)
        notify(queue; all=false)
        v
    end
end

"""
    release!(wp::WorkerPool)

Marks this worker as no longer active (e.g. aborting early) and wakes every
other worker so they can re-check the completion condition.
"""
function release!(wp::WorkerPool)
    lock(wp.queue) do
        @atomic wp.active -= 1
        notify(wp.queue; all=true)
    end
end

"""
    run!(f, wp::WorkerPool)

Julia has no deterministic destructors, so the RAII-style guarantee here is
built the idiomatic way: a `try`/`finally` hidden inside a do-block-taking
function (the same shape `open(file) do io ... end` uses to guarantee `close`).

Registers this worker with `wp`, then calls `f(cell)` for each item
[`fetch!`](@ref) returns until there's none left or `f` returns `false`.
Releases the worker's slot on every exit -- normal, early stop via `false`, or
exception -- so callers never call `register!`/`release!` themselves.
"""
function run!(f, wp::WorkerPool)
    register!(wp)
    cell = nothing
    try
        cell = fetch!(wp)
        while !isnothing(cell)
            f(cell) === false && return
            cell = fetch!(wp)
        end
    finally
        isnothing(cell) || release!(wp)
    end
end
