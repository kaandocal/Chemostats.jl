using Test
using ArgCheck
using OrdinaryDiffEqTsit5
using Random
using Chemostats

include("models/exponential.jl")

@testset "Snapshot math" begin
    s1 = Chemostats.Snapshot(0.0, 10, 10, 0.0)
    s2 = Chemostats.Snapshot(1.0, 20, 20, 0.0)

    @test Chemostats.est_logN(s1) == log(10)
    @test Chemostats.est_N(s1) ≈ 10
    @test Chemostats.est_Λ(s1, s2) ≈ log(2)
    @test isnan(Chemostats.est_Λ(s2, s1))   # before.t >= after.t is undefined

    # log_f != 0 means we're only tracking a fraction of the true population.
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

@testset "die!/divide!" begin
    cell = ExponentialCell((; id = 1))
    @test Chemostats.get_state(cell) == Chemostats.CellState.Newborn

    Chemostats.die!(cell)
    @test Chemostats.get_state(cell) == Chemostats.CellState.Newborn

    Chemostats.divide!(cell)
    @test Chemostats.get_state(cell) == Chemostats.CellState.Newborn
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
