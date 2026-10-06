using Test
using ArgCheck
using OrdinaryDiffEqTsit5
using Random
using SciMLBase
using Chemostats

include("models/exponential.jl")

@testset "Snapshot math" begin
    s1 = Chemostats.Snapshot(0.0, 10, 10, 0.0)
    s2 = Chemostats.Snapshot(1.0, 20, 20, 0.0)

    @test Chemostats.est_logN(s1) == log(10)
    @test Chemostats.est_N(s1) ≈ 10
    @test Chemostats.est_Λ(s1, s2) ≈ log(2)
    @test isnan(Chemostats.est_Λ(s2, s1))

    s3 = Chemostats.Snapshot(0.0, 10, 10, log(2))
    @test Chemostats.est_N(s3) ≈ 20
end

@testset "get_snapshot" begin
    snaps = [
        Chemostats.Snapshot(0.0, 10, 10, 0.0),
        Chemostats.Snapshot(1.0, 10, 10, 0.0),
        Chemostats.Snapshot(1.0, 15, 15, 0.0),
    ]

    @test_throws ArgumentError Chemostats.get_snapshot(snaps, 5.0)
    @test Chemostats.get_snapshot(snaps, 1.0).N == 15
    @test Chemostats.get_snapshot(snaps).t == 1.0
end

@testset "Lax constructor" begin
    @test_throws ArgumentError Chemostats.Lax(10, 1.0; tol = 0.5)
    @test_throws ArgumentError Chemostats.Lax(10, 1.0; β = 1.5)
    @test_throws ArgumentError Chemostats.Lax(10, 1.0; β = -0.1)
    @test_throws ArgumentError Chemostats.Lax(10, -1.0)
    @test Chemostats.Lax(10, 1.0) isa Chemostats.Lax
end

@testset "Lax est_Λ_curr matches a known growth rate across β" begin
    Λ_true = 0.3
    τ = 0.5
    N0 = 1.0e6   # large, so rounding N to an Int is negligible

    ts = 0:τ:10
    snaps = [ Chemostats.Snapshot(t, round(Int, N0 * exp(Λ_true * t)), 0, 0.0) for t in ts ]
    t = ts[end]

    Λ_direct = Chemostats.est_Λ(snaps[1], snaps[end])
    @test Λ_direct ≈ Λ_true atol = 1e-6

    for β in (0.0, 0.3, 0.5, 0.8, 1.0)
        alg = Chemostats.Lax(50, τ; β)
        Λ_curr = Chemostats.est_Λ_curr(snaps, t, alg)

        @test Λ_curr ≈ Λ_true atol = 1e-6
        @test Λ_curr ≈ Λ_direct atol = 1e-6
    end
end

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
    @test isnothing(cell1.int)   # properly finalized

    prob2 = ODEProblem(f_exp, 0.1, (0., 0.), (; id = 2); callback = cb_exp)
    cell2 = DECell(prob2, Tsit5(), divide_exp)
    Chemostats.init_cell!(cell2)
    Chemostats.step!(cell2, 1.0, nothing)
    @test Chemostats.get_state(cell2) == Chemostats.CellState.EndOfLife

    Chemostats.kill!(cell2, Chemostats.get_curr_t(cell2))
    @test Chemostats.get_state(cell2) == Chemostats.CellState.Killed
    @test isnothing(cell2.int)
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

@testset "Chemostat option copy" begin
    cells = make_population(3)
    chem_shared = Chemostat(cells)             # default: copy=false
    chem_copied = Chemostat(cells; copy = true)

    @test chem_shared.pop === cells
    @test chem_copied.pop !== cells
end

@testset "Offspring count edge cases" begin
    @testset "one child" begin
        divide_single(int) = [ (u0 = randexp(), p = int.p) ]
        prob = ODEProblem(f_exp, randexp(), (0., 0.), (; id = 1); callback = cb_exp)
        chem = Chemostat([ DECell(prob, Tsit5(), divide_single) ])

        Chemostats.simulate!(chem, 5.0, Chemostats.Direct())
        @test chem.snaps[end].N == 1   # each division replaces 1 cell with 1
    end
end

@testset "Snapshot/Chemostat show methods" begin
    snap = Chemostats.Snapshot(1.5, 20, 25, 0.1)
    str = sprint(show, snap)
    @test occursin("t=1.5", str)
    @test occursin("N=20", str)
    @test occursin("nsim=25", str)

    chem = Chemostat(make_population(4))
    chem_str = sprint(show, chem)
    @test occursin("Chemostat(", chem_str)
    @test occursin("t=0.0", chem_str)
end

@testset "est_logN(chem) matches the Snapshot-based version" begin
    chem = Chemostat(make_population(4))
    chem.snaps[1] = Chemostats.Snapshot(0.0, 4, 4, log(2))   # pretend partial sampling

    @test Chemostats.est_logN(chem) ≈ Chemostats.est_logN(chem.snaps[end])
    @test Chemostats.est_logN(chem) ≈ log(4) + log(2)
end

@testset "step! warns and no-ops on a non-alive cell" begin
    cell = ExponentialCell((; id = 1))
    @test Chemostats.get_state(cell) == Chemostats.CellState.Newborn

    # step! requires Alive; calling it on a Newborn cell should warn and
    # leave the cell untouched, not error.
    @test_logs (:warn, r"Tried to simulate cell in state") Chemostats.step!(cell, 1.0, nothing)
    @test Chemostats.get_state(cell) == Chemostats.CellState.Newborn
end

@testset "step! rejects unsupported ensemble algorithms" begin
    chem = Chemostat(make_population(4))
    alg = Chemostats.Direct()
    int = Chemostats.PopIntegrator(chem, alg, EnsembleSerial())
    Chemostats.init!(alg, int)

    @test_throws "not supported" Chemostats.step!(int, 5.0, SciMLBase.EnsembleDistributed())
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

@testset "default_ensalg with is_parallel and nthreads" begin
    for alg in (Chemostats.Forward(10), Chemostats.Strict(10), Chemostats.Thin(0.1), Chemostats.Lax(10, 1.0))
        ensalg = Chemostats.default_ensalg(alg)
        expect_threaded = Chemostats.is_parallel(alg) && Threads.nthreads() > 1
        @test (ensalg isa EnsembleThreads) == expect_threaded
    end
end
