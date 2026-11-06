using Test
using ArgCheck
using OrdinaryDiffEqTsit5
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

@testset "_get_snapshot" begin
    snaps = [
        Chemostats.Snapshot(0.0, 10, 10, 0.0),
        Chemostats.Snapshot(1.0, 10, 10, 0.0),
        Chemostats.Snapshot(1.0, 15, 15, 0.0),
    ]

    @test_throws ArgumentError Chemostats._get_snapshot(snaps, 5.0)
    @test Chemostats._get_snapshot(snaps, 1.0).N == 15
    @test Chemostats._get_snapshot(snaps).t == 1.0
end

@testset "Chemostat start time" begin
    cell(t0) = DECell(ODEProblem(f_exp, 1.0, (0., 0.), (; id = 1); callback = cb_exp), Tsit5(), divide_exp; t0)

    chem = Chemostat([ cell(5.) for _ in 1:3 ])
    @test Chemostats.get_curr_t(chem) == 5. && only(chem.snaps).t == 5.
    @test_throws ArgumentError Chemostat([ cell(0.), cell(5.) ])
    @test_throws ArgumentError Chemostat(DECell[])
end

@testset "Chemostat option copy" begin
    cells = make_population(3)
    chem_shared = Chemostat(cells)             # default: copy=false
    chem_copied = Chemostat(cells; copy = true)

    @test chem_shared.pop === cells
    @test chem_copied.pop !== cells
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

@testset "Simulation resuming" begin
    K = 4
    chem = Chemostat(make_population(K))

    Chemostats.simulate!(chem, 1.0, Chemostats.Direct())
    @test chem.snaps[end].t == 1.0

    Chemostats.simulate!(chem, 2.0, Chemostats.Direct())
    @test chem.snaps[end].t == 2.0
    @test length(chem.snaps) >= 3
end

@testset "saveat records snapshots" begin
    K = 4
    chem = Chemostat(make_population(K))
    tmax = 3.0

    Chemostats.simulate!(chem, tmax, Chemostats.Direct(); saveat = [1.0, 2.0, 3.0])
    ts = [ s.t for s in chem.snaps ]

    @test issorted(ts)
    @test all(t -> t in ts, (1.0, 2.0, 3.0))
end

@testset "simulate creates and simulates a Chemostat" begin
    cells = make_population(3)
    chem = Chemostats.simulate(cells, 1.0, Chemostats.Direct(); copy = true, save_mode = :full)

    @test chem.snaps[end].t == 1.0
    @test !isnothing(chem.tree)
    @test all(c -> Chemostats.get_curr_t(c) == 0.0, cells)     # copied cells are left untouched
end
