using Test
using ArgCheck
using OrdinaryDiffEqTsit5
using Random
using SciMLBase
using SymbolicIndexingInterface: SymbolicIndexingInterface
using Distributions: Uniform
using Chemostats
using Chemostats: Offspring, ExactPartition

include("models/exponential.jl")
include("models/jumpcell.jl")

const iV = CartesianIndex(1)

@testset "die!/divide!" begin
    cell = ExponentialCell((; id = 1))
    @test Chemostats.get_state(cell) == Chemostats.CellState.Newborn

    @test_logs (:warn, r"Called `die!`") Chemostats.die!(cell)
    @test Chemostats.get_state(cell) == Chemostats.CellState.Newborn

    @test_logs (:warn, r"tried to divide") Chemostats.divide!(cell)
    @test Chemostats.get_state(cell) == Chemostats.CellState.Newborn
end

@testset "kill! always kills, regardless of prior state" begin
    # Alive
    prob1 = ODEProblem(f_exp, 5.0, (0., 0.), (; id = 1); callback = cb_exp)
    cell1 = DECell(prob1, Tsit5(), divide_exp)
    Chemostats.init_cell!(cell1)
    Chemostats.step!(cell1, 2.0, nothing)
    @test Chemostats.get_state(cell1) == Chemostats.CellState.Alive

    Chemostats.kill!(cell1, Chemostats.get_curr_t(cell1))
    @test Chemostats.get_state(cell1) == Chemostats.CellState.Killed
    @test SciMLBase.done(cell1.int)   # integrator terminated

    prob2 = ODEProblem(f_exp, 0.1, (0., 0.), (; id = 2); callback = cb_exp)
    cell2 = DECell(prob2, Tsit5(), divide_exp)
    Chemostats.init_cell!(cell2)
    Chemostats.step!(cell2, 1.0, nothing)
    @test Chemostats.get_state(cell2) == Chemostats.CellState.EndOfLife

    Chemostats.kill!(cell2, Chemostats.get_curr_t(cell2))
    @test Chemostats.get_state(cell2) == Chemostats.CellState.Killed
    @test SciMLBase.done(cell2.int)
end

@testset "clone independence" begin
    prob = ODEProblem(f_exp, 5.0, (0., 0.), (; id = 1); callback = cb_exp)
    cell = DECell(prob, Tsit5(), divide_exp)
    Chemostats.init_cell!(cell)
    Chemostats.step!(cell, 2.0, nothing)

    clone = Chemostats.clone_cell(cell, 2.0)

    @test clone !== cell
    @test Chemostats.get_curr_t(clone) == Chemostats.get_curr_t(cell) == 2.0
    @test Chemostats.get_state(clone) == Chemostats.CellState.Newborn   # clone starts fresh

    # Stepping the clone further must not retroactively affect the original.
    Chemostats.init_cell!(clone)
    Chemostats.step!(clone, 1.0, nothing)
    @test Chemostats.get_curr_t(cell) == 2.0
    @test Chemostats.get_curr_t(clone) == 3.0
end

@testset "clone_cell rejects out-of-range t" begin
    prob = ODEProblem(f_exp, 5.0, (0., 0.), (; id = 1); callback = cb_exp)
    cell = DECell(prob, Tsit5(), divide_exp)
    Chemostats.init_cell!(cell)
    Chemostats.step!(cell, 2.0, nothing)
    t_curr = Chemostats.get_curr_t(cell)   # 2.0 -- not yet simulated past this

    @test_throws ArgCheck.CheckError Chemostats.clone_cell(cell, t_curr + 1.0)  # future
    @test_throws ArgCheck.CheckError Chemostats.clone_cell(cell, -1.0)          # before trajectory start

    # Within range: must still succeed (not an off-by-one on the boundary).
    clone = Chemostats.clone_cell(cell, t_curr)
    @test Chemostats.get_curr_t(clone) == t_curr
end

@testset "reset_t = true" begin
    prob = ODEProblem(f_exp, randexp(), (0., 0.), (; id = 1); callback = cb_exp)
    cell = DECell(prob, Tsit5(), divide_exp; reset_t = true)

    chem = Chemostat([cell])
    tmax = 5.0
    Chemostats.simulate!(chem, tmax, Chemostats.Direct())

    @test chem.snaps[end].t == tmax
    @test all(c -> Chemostats.get_curr_t(c) == tmax, chem.pop)
end

@testset "Integer start times with a JumpProblem (reset_t=$reset_t)" for reset_t in (true, false)
    chem = Chemostat([ JumpCell(; k = 5.0, reset_t, t0 = 0) for _ in 1:3 ])
    Chemostats.simulate!(chem, 3.0, Chemostats.Direct())

    @test chem.status == SciMLBase.ReturnCode.Success
    @test all(c -> Chemostats.get_curr_t(c) == 3.0, chem.pop)
end

@testset "Offspring count edge cases" begin
    @testset "one child" begin
        divide_single = Chemostats.Offspring(1, CartesianIndex() => Exponential())
        prob = ODEProblem(f_exp, randexp(), (0., 0.), (; id = 1); callback = cb_exp)
        chem = Chemostat([ DECell(prob, Tsit5(), divide_single) ])

        Chemostats.simulate!(chem, 5.0, Chemostats.Direct())
        @test chem.snaps[end].N == 1   # each division replaces 1 cell with 1
    end
end

@testset "step! warns and no-ops on a non-alive cell" begin
    cell = ExponentialCell((; id = 1))
    @test Chemostats.get_state(cell) == Chemostats.CellState.Newborn

    # step! requires Alive; calling it on a Newborn cell should warn and
    # leave the cell untouched, not error.
    @test_logs (:warn, r"Tried to simulate cell in state") Chemostats.step!(cell, 1.0, nothing)
    @test Chemostats.get_state(cell) == Chemostats.CellState.Newborn
end

@testset "step! error kills the cell" begin
    fail_cb = SciMLBase.DiscreteCallback(
        (u, t, int) -> t > 0.5, int -> SciMLBase.terminate!(int, SciMLBase.ReturnCode.Unstable)
    )
    prob = ODEProblem(f_exp, 10.0, (0., 0.), (; id = 1); callback = CallbackSet(cb_exp, fail_cb))
    cell = DECell(prob, Tsit5(), divide_exp)
    Chemostats.init_cell!(cell)

    @test_logs (:warn, r"Cell solver errored") Chemostats.step!(cell, 1.0, nothing)
    @test Chemostats.get_state(cell) == Chemostats.CellState.Killed
end

@testset "clone_cell independence (remake = $remake)" for remake in (false, true)
    ref = JumpCell(; k = 10.0, V_d = 2.0, remake)
    Chemostats.init_cell!(ref)
    clones = [ Chemostats.clone_cell(ref, 0.0) for _ in 1:5 ]

    @test all(c.prob === ref.prob for c in clones)   # cloning is cheap: prob is shared

    aggs = [ jumpcell_aggregator(c.int) for c in [ref; clones] ]
    @test allunique(objectid.(aggs))

    set_V_d = SymbolicIndexingInterface.setp(ref.int, :V_d)
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

@testset "Division of a DECell" begin
    # The volume grows exponentially, and the cell divides once it has doubled
    grow!(du, u, p, t) = (du[1] = u[1])
    cb = Chemostats.DivideCallback((u, t, int) -> u[1] - 2.0)
    prob = ODEProblem(grow!, [1.0], (0., Inf); callback = cb)

    function children(divide; seed = 1)
        cell = DECell(prob, Tsit5(), divide; rng = Xoshiro(seed))
        Chemostats.init_cell!(cell)
        Chemostats.step!(cell, 10.0, nothing)
        Chemostats.get_children(cell)
    end
    volume(cell) = cell.int.u[1]

    kids = children(Offspring(2, iV => ExactPartition()))
    @test length(kids) == 2
    @test all(k -> isapprox(volume(k), 1; atol = 1e-6), kids)

    @test length(children(int -> Offspring(3, iV => ExactPartition()))) == 3
    @test isempty(children(int -> nothing))
    @test isempty(children(Offspring(0)))

    # Random rules use the RNG of the dividing cell
    sample(seed) = volume.(children(Offspring(2, iV => Uniform(0.5, 1.5)); seed))
    @test sample(1) == sample(1)
    @test sample(1) != sample(2)
end
