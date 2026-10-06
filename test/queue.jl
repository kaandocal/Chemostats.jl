using Test
using Random
using Chemostats

# A minimal value type for exercising ThreadedQueue/WorkerPool in isolation,
# with no dependency on DECell or the rest of the simulation engine.
struct TimedVal
    t::Float64
    id::Int
end

Chemostats.get_curr_t(v::TimedVal) = v.t
Chemostats.clone_cell(v::TimedVal, t) = TimedVal(t, v.id)

mk_queue(ts) = Chemostats.ThreadedQueue([ TimedVal(t, i) for (i, t) in enumerate(ts) ])

@testset "ThreadedQueue: single-threaded basics" begin
    q = mk_queue([3.0, 1.0, 2.0])
    @test length(q) == 3
    @test !isempty(q)
    @test first(q).t == 1.0

    @test Chemostats.trypop!(q).t == 1.0
    @test Chemostats.trypop!(q).t == 2.0
    @test Chemostats.trypop!(q).t == 3.0
    @test isnothing(Chemostats.trypop!(q))
    @test isempty(q)
end

@testset "ThreadedQueue: push! respects heap order" begin
    q = mk_queue(Float64[])
    for t in [5.0, 1.0, 4.0, 2.0, 3.0]
        push!(q, TimedVal(t, 0))
    end

    popped = [ Chemostats.trypop!(q).t for _ in 1:5 ]
    @test popped == [1.0, 2.0, 3.0, 4.0, 5.0]
end

@testset "ThreadedQueue: iter_unsafe / _append! / extract_queue!" begin
    q = mk_queue([1.0, 2.0])
    @test Set(v.t for v in Chemostats.iter_unsafe(q)) == Set([1.0, 2.0])

    Chemostats._append!(q, [ TimedVal(3.0, 0), TimedVal(4.0, 0) ])
    @test length(q) == 4

    out = TimedVal[]
    Chemostats.extract_queue!(out, q)
    @test isempty(q)
    @test Set(v.t for v in out) == Set([1.0, 2.0, 3.0, 4.0])
end

@testset "ThreadedQueue: _pop_random!/_clone_random! preserve the rest" begin
    q = mk_queue(1.0:10.0)

    popped = Chemostats._pop_random!(q)
    @test length(q) == 9
    remaining = Set(v.t for v in Chemostats.iter_unsafe(q))
    @test popped.t ∉ remaining
    @test remaining == setdiff(Set(1.0:10.0), [popped.t])

    source, clone = Chemostats._clone_random!(q, 99.0)
    @test length(q) == 10
    @test clone.t == 99.0
    @test clone.id == source.id
end

@testset "ThreadedQueue: concurrent push!/trypop! lose nothing, duplicate nothing" begin
    q = mk_queue(Float64[])
    nproducers, per_producer = 8, 500
    total = nproducers * per_producer

    @sync for p in 1:nproducers
        Threads.@spawn for i in 1:per_producer
            push!(q, TimedVal(rand(), (p - 1) * per_producer + i))
        end
    end
    @test length(q) == total

    seen = zeros(Bool, total)
    nconsumers = 8
    @sync for _ in 1:nconsumers
        Threads.@spawn begin
            while true
                v = Chemostats.trypop!(q)
                isnothing(v) && break
                @test !seen[v.id]   # no duplicate delivery
                seen[v.id] = true
            end
        end
    end
    @test all(seen)
    @test isempty(q)
end

@testset "WorkerPool: fetch! hands out pushed items" begin
    q = mk_queue([1.0, 2.0])
    wp = Chemostats.WorkerPool(q)
    Chemostats.register!(wp)

    @test Chemostats.fetch!(wp).t == 1.0
    @test Chemostats.fetch!(wp).t == 2.0
end

@testset "WorkerPool: fetch! returns nothing once empty and nobody else active" begin
    q = mk_queue(Float64[])
    wp = Chemostats.WorkerPool(q)
    Chemostats.register!(wp)

    @test isnothing(Chemostats.fetch!(wp))
end

@testset "WorkerPool: fetch! blocks while another worker is active, wakes on push!" begin
    q = mk_queue(Float64[])
    wp = Chemostats.WorkerPool(q)
    Chemostats.register!(wp)   # worker A: about to go idle
    Chemostats.register!(wp)   # worker B: still "active" from A's perspective

    done = Threads.Atomic{Bool}(false)
    task = Threads.@spawn begin
        v = Chemostats.fetch!(wp)   # should block: queue empty, wp.active > 0
        done[] = true
        v
    end

    sleep(0.2)
    @test !done[]   # still blocked

    push!(q, TimedVal(42.0, 1))
    v = fetch(task)
    @test done[]
    @test v.t == 42.0
end

@testset "WorkerPool: fetch! blocks then wakes via release! (no more work coming)" begin
    q = mk_queue(Float64[])
    wp = Chemostats.WorkerPool(q)
    Chemostats.register!(wp)
    Chemostats.register!(wp)

    task = Threads.@spawn Chemostats.fetch!(wp)
    sleep(0.2)
    @test !istaskdone(task)

    Chemostats.release!(wp)   # the other worker gives up; queue still empty
    @test isnothing(fetch(task))
end

@testset "WorkerPool: dynamic work (like cell division) terminates with no loss" begin
    q = mk_queue([0.0])
    wp = Chemostats.WorkerPool(q)

    produced = Threads.Atomic{Int}(1)   # the seed item above
    consumed = Threads.Atomic{Int}(0)

    nworkers = 8
    @sync for _ in 1:nworkers
        Threads.@spawn begin
            Chemostats.register!(wp)
            while true
                v = Chemostats.fetch!(wp)
                isnothing(v) && break
                Threads.atomic_add!(consumed, 1)

                # ~40% chance of producing 2 children, else none -- a finite process w.h.p.
                if v.id < 12 && rand() < 0.4
                    Threads.atomic_add!(produced, 2)
                    push!(q, TimedVal(0.0, 2v.id))
                    push!(q, TimedVal(0.0, 2v.id + 1))
                end
            end
        end
    end

    @test consumed[] == produced[]
    @test isempty(q)
end

@testset "WorkerPool.run!: processes every item, releases normally" begin
    q = mk_queue([1.0, 2.0, 3.0])
    wp = Chemostats.WorkerPool(q)

    seen = Float64[]
    Chemostats.run!(wp) do cell
        push!(seen, cell.t)
    end

    @test seen == [1.0, 2.0, 3.0]
    @test wp.active == 0   # register!'s +1 was given back by fetch!'s own "done" path
end

@testset "WorkerPool.run!: releases on early stop (f returns false)" begin
    q = mk_queue([1.0, 2.0, 3.0])
    wp = Chemostats.WorkerPool(q)

    seen = Float64[]
    Chemostats.run!(wp) do cell
        push!(seen, cell.t)
        cell.t == 2.0 ? false : nothing
    end

    @test seen == [1.0, 2.0]   # stopped after the second item, third left untouched
    @test length(q) == 1
    @test wp.active == 0       # released despite stopping early
end

@testset "WorkerPool.run!: releases on exception, which still propagates" begin
    q = mk_queue([1.0, 2.0])
    wp = Chemostats.WorkerPool(q)

    seen = Float64[]
    @test_throws ErrorException Chemostats.run!(wp) do cell
        push!(seen, cell.t)
        cell.t == 1.0 && error("boom")
    end

    @test seen == [1.0]
    @test wp.active == 0   # released even though f threw
end

@testset "WorkerPool.run!: dynamic work (like cell division) terminates with no loss" begin
    q = mk_queue([0.0])
    wp = Chemostats.WorkerPool(q)

    produced = Threads.Atomic{Int}(1)
    consumed = Threads.Atomic{Int}(0)

    nworkers = 8
    @sync for _ in 1:nworkers
        Threads.@spawn Chemostats.run!(wp) do v
            Threads.atomic_add!(consumed, 1)
            if v.id < 12 && rand() < 0.4
                Threads.atomic_add!(produced, 2)
                push!(q, TimedVal(0.0, 2v.id))
                push!(q, TimedVal(0.0, 2v.id + 1))
            end
        end
    end

    @test consumed[] == produced[]
    @test isempty(q)
    @test wp.active == 0
end
