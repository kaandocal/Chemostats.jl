using Test
using Random
using SciMLBase
using Chemostats

include("models/exponential.jl")
include("models/jumpcell.jl")

@testset "clone_cell independence" begin
    ref = JumpCell(; k = 10.0, V_d = 2.0)
    Chemostats.init_cell!(ref)
    clones = [ Chemostats.clone_cell(ref, 0.0) for _ in 1:5 ]

    @test all(c.prob === ref.prob for c in clones)   # cloning is cheap: prob is shared

    aggs = [ jumpcell_aggregator(c.int) for c in [ref; clones] ]
    @test allunique(objectid.(aggs))

    set_V_d = SciMLBase.setp(ref.int, :V_d)
    for (i, c) in enumerate(clones)
        set_V_d(c.int, 2.0 + i)
    end
    @test [ c.int.ps[:V_d] for c in clones ] == 2.0 .+ (1:5)
    @test allunique(objectid.([ c.int.p for c in clones ]))
    @test ref.int.ps[:V_d] == 2.0   # mutating clones must never touch ref

    # A clone made after its siblings were mutated must still see ref's original value.
    late_clone = Chemostats.clone_cell(ref, 0.0)
    @test late_clone.int.ps[:V_d] == 2.0

    u_before = [ copy(c.int.u) for c in clones[2:end] ]
    Chemostats.init_cell!(clones[1])
    Chemostats.step!(clones[1], 0.1, nothing)
    @test all(c.int.u == u for (c, u) in zip(clones[2:end], u_before))
end

@testset "clone_cell independence (remake = true)" begin
    ref = JumpCell(; k = 10.0, V_d = 2.0, remake = true)
    Chemostats.init_cell!(ref)
    clones = [ Chemostats.clone_cell(ref, 0.0) for _ in 1:5 ]

    @test isnothing(ref.prob)   # remake = true cells never hold the original prob

    aggs = [ jumpcell_aggregator(c.int) for c in [ref; clones] ]
    @test allunique(objectid.(aggs))

    set_V_d = SciMLBase.setp(ref.int, :V_d)
    for (i, c) in enumerate(clones)
        set_V_d(c.int, 2.0 + i)
    end
    @test [ c.int.ps[:V_d] for c in clones ] == 2.0 .+ (1:5)
    @test allunique(objectid.([ c.int.p for c in clones ]))
    @test ref.int.ps[:V_d] == 2.0

    u_before = [ copy(c.int.u) for c in clones[2:end] ]
    Chemostats.init_cell!(clones[1])
    Chemostats.step!(clones[1], 0.1, nothing)
    @test all(c.int.u == u for (c, u) in zip(clones[2:end], u_before))
end

@testset "clone_cell u independence" begin
    for ctor in (() -> ExponentialCell((; id = 1)), () -> JumpCell(; k = 20.0, V_d = 2.0))
        ref = ctor()
        Chemostats.init_cell!(ref)
        clone1 = Chemostats.clone_cell(ref, 0.0)
        clone2 = Chemostats.clone_cell(ref, 0.0)

        u2_before = copy(clone2.int.u)
        Chemostats.init_cell!(clone1)
        Chemostats.step!(clone1, 0.1, nothing)

        @test clone2.int.u == u2_before
        @test clone2.int.t == 0.0
    end
end

@testset "get_children aliasing" begin
    for ctor in (() -> ExponentialCell((; id = 1)), () -> JumpCell(; k = 20.0, V_d = 2.0))
        cell = ctor()
        Chemostats.init_cell!(cell)
        Chemostats.step!(cell, 20.0, nothing)
        @test Chemostats.get_state(cell) == Chemostats.CellState.EndOfLife

        children = Chemostats.get_children(cell)
        @test length(children) == 2
        d1, d2 = children

        u2_before = copy(d2.int.u)
        if !(d1.int.p isa Union{NamedTuple, Number})
            @test d1.int.p !== d2.int.p
        end

        # Regression check: `reinit` used to alias `old_int.u` (the
        # *parent's* live state array) as the base for each child's `u0`.
        Chemostats.init_cell!(d1)
        Chemostats.step!(d1, 0.1, nothing)
        @test d2.int.u == u2_before
    end
end

@testset "get_children independence" begin
    cell = JumpCell(; k = 50.0, V_d = 2.0)
    Chemostats.init_cell!(cell)
    Chemostats.step!(cell, 20.0, nothing)
    @test Chemostats.get_state(cell) == Chemostats.CellState.EndOfLife

    for daughter in Chemostats.get_children(cell)
        Chemostats.init_cell!(daughter)
        Chemostats.step!(daughter, 20.0, nothing)
        # Regression check: cells produced via `get_children`/`reinit` on a
        # `JumpProblem` used to silently lose their jump schedule and callback.
        @test Chemostats.get_state(daughter) == Chemostats.CellState.EndOfLife
        @test daughter.int[:P] >= 0
    end
end

