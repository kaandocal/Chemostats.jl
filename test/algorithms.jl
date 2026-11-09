using Test
using OrdinaryDiffEqTsit5
using SciMLBase
using Chemostats

include("models/deterministic.jl")
include("models/exponential.jl")

ENSALGS = (EnsembleSerial(), EnsembleThreads())

# Exact bookkeeping with deterministic cells, which double at rate Λ_gt
Λ_gt = 1.4
prob = Deterministic(Λ_gt)
alg = Tsit5()

for ensalg in [ EnsembleSerial(), EnsembleThreads() ]
    @testset "Direct() with $ensalg" begin 
        tmax = 5. / Λ_gt

        chem = Chemostat([ DECell(prob, alg, divide_det) ])
        Chemostats.simulate!(chem, tmax, Chemostats.Direct(), ensalg) 
        snap = chem.snaps[end]

        @test snap.t == tmax
        @test snap.N == 2 ^ floor(Int, tmax * Λ_gt)
        @test snap.nsim == 2 * snap.N - 1
        @test Chemostats.est_N(snap) == snap.N
    end 
end 

### 

tmax = 20.5 

@testset "Strict($L, clone=$clone)" for L in [ 1, 10 ], clone in (true, false)
    chem = Chemostat([ DECell(prob, alg, divide_det) for i in 1:L ])

    @test_throws "Strict does not support parallelisation" Chemostats.simulate!(chem, tmax, Chemostats.Strict(L, clone), EnsembleThreads())

    Chemostats.simulate!(chem, tmax, Chemostats.Strict(L, clone))
    snap = chem.snaps[end]

    @test snap.t == tmax
    @test snap.N == L

    if L == 1
        @test Chemostats.est_logN(snap) ≈ floor(Int, tmax * Λ_gt) * log(1 + 1/L) * L + log(L)
    end
end

for ensalg in [ EnsembleSerial(), EnsembleThreads() ]
    @testset "Forward($L) with $ensalg" for L in [ 1, 10 ]
        chem = Chemostat([ DECell(prob, alg, divide_det) for i in 1:L ])
        
        @test_throws DimensionMismatch Chemostats.simulate!(chem, tmax, Chemostats.Forward(L+1), ensalg)

        Chemostats.simulate!(chem, tmax, Chemostats.Forward(L), ensalg)
        snap = chem.snaps[end]

        @test snap.t == tmax
        @test snap.N == L
        @test isnan(snap.log_f)
        @test snap.nsim == ceil(Int, tmax * Λ_gt) * L
    end
end 

### Bookkeeping with stochastic cells

@testset "Direct()/Thin(0) exact bookkeeping" for ensalg in ENSALGS
    K = 5
    chem = Chemostat(make_population(K))
    tmax = 3.0

    Chemostats.simulate!(chem, tmax, Chemostats.Direct(), ensalg)
    snap = chem.snaps[end]

    @test snap.t == tmax
    @test length(chem.pop) == snap.N
    @test snap.N >= K                       # no cell dies in this model
    @test snap.nsim == 2*snap.N - K         # forest of K full binary trees
    @test snap.log_f == 0.
    @test Chemostats.est_N(snap) ≈ snap.N
    @test all(cell -> Chemostats.get_curr_t(cell) == tmax, chem.pop)

    # Every starting lineage must have descendants
    ids = Set(cell_params(cell).id for cell in chem.pop)
    @test ids == Set(1:K)
end

@testset "Forward(L) preserves lineages" for ensalg in ENSALGS
    K = 6
    chem = Chemostat(make_population(K))
    tmax = 3.0

    Chemostats.simulate!(chem, tmax, Chemostats.Forward(K), ensalg)
    snap = chem.snaps[end]

    @test snap.t == tmax
    @test snap.N == K
    @test length(chem.pop) == K
    @test isnan(snap.log_f)
    @test snap.nsim >= K

    ids = sort([ cell_params(cell).id for cell in chem.pop ])
    @test ids == 1:K
end

@testset "Strict population size (clone=$clone)" for clone in (true, false)
    K = 10
    chem = Chemostat(make_population(K))
    tmax = 3.0

    Chemostats.simulate!(chem, tmax, Chemostats.Strict(K, clone))
    snap = chem.snaps[end]

    @test snap.t == tmax
    @test snap.N == K
    @test length(chem.pop) == K

    @test all(cell -> cell_params(cell).id in 1:K, chem.pop)
end

@testset "Lax population size" for ensalg in ENSALGS
    L = 20
    chem = Chemostat(make_population(L))
    tmax = 5.0
    alg = Chemostats.Lax(L, 0.5)

    Chemostats.simulate!(chem, tmax, alg, ensalg; Nmax = 10_000)
    snap = chem.snaps[end]

    @test snap.t == tmax
    @test 0 < snap.N < 20*L
    @test length(chem.pop) == snap.N

    @test all(cell -> cell_params(cell).id in 1:L, chem.pop)
end

@testset "Thin(δ) bookkeeping" for ensalg in ENSALGS
    K = 30
    δ = 0.3
    tmax = 2.0
    chem = Chemostat(make_population(K); save_mode = :full)

    Chemostats.simulate!(chem, tmax, Chemostats.Thin(δ), ensalg)
    snap = chem.snaps[end]

    @test snap.t == tmax
    @test snap.N == length(chem.pop)

    @test snap.log_f ≈ δ*tmax

    tree = chem.tree
    N_dead = count(id -> !Chemostats.is_alive(tree[id]) && isempty(Chemostats.children(tree, id)), eachindex(tree))
    @test snap.nsim == 2*(snap.N + N_dead) - K
    @test N_dead > 0   # sanity: δ=0.3 over tmax=2 should cull at least a few
end

@testset "Lax constructor" begin
    @test_throws ArgumentError Chemostats.Lax(10, 1.0; tol = 0.5)
    @test_throws ArgumentError Chemostats.Lax(10, 1.0; smoothing = 1.5)
    @test_throws ArgumentError Chemostats.Lax(10, 1.0; smoothing = -0.1)
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

    Λ_direct = Chemostats.est_growth_rate(snaps[1], snaps[end])
    @test Λ_direct ≈ Λ_true atol = 1e-6

    for β in (0.0, 0.3, 0.5, 0.8, 1.0)
        alg = Chemostats.Lax(50, τ; smoothing = β)
        Λ_curr = Chemostats.est_Λ_curr(snaps, t, alg)

        @test Λ_curr ≈ Λ_true atol = 1e-6
        @test Λ_curr ≈ Λ_direct atol = 1e-6
    end
end

@testset "Lax and est_growth_rate count time from the first snapshot" begin
    Λ_true = 0.3
    τ = 0.5
    ts = 5:τ:15
    snaps = [ Chemostats.Snapshot(t, round(Int, 1.0e6 * exp(Λ_true * t)), 0, 0.0) for t in ts ]

    @test Chemostats.est_growth_rate(snaps) ≈ Λ_true atol = 1e-6
    @test Chemostats.est_Λ_curr(snaps, ts[end], Chemostats.Lax(50, τ)) ≈ Λ_true atol = 1e-6
    @test Chemostats.get_thinning_rate(snaps, ts[1], Chemostats.Lax(50, τ)) == 0          # less than τ after the start
end

@testset "Lax falls back to nearby snapshots or sets δ = 0" begin
    Λ_true = 0.3
    τ = 0.5
    N0 = 1.0e6
    alg = Chemostats.Lax(50, τ)
    mksnap(t) = Chemostats.Snapshot(t, round(Int, N0 * exp(Λ_true * t)), 0, 0.0)

    ts = 0:τ:10
    t = ts[end]

    # Snapshot at t - τ shifted within τ/2: estimate still matches
    snaps = [ mksnap(s) for s in ts ]
    snaps[end-1] = mksnap(t - τ + 0.2τ)
    @test Chemostats.est_Λ_curr(snaps, t, alg) ≈ Λ_true atol = 1e-3

    # No snapshot within τ/2 of t - τ: estimate undefined, δ = 0
    snaps = [ mksnap(s) for s in ts if s != t - τ ]
    @test isnan(Chemostats.est_Λ_curr(snaps, t, alg))
    @test Chemostats.get_thinning_rate(snaps, t, alg) == 0
end

@testset "Lax snapshot lookup matches a brute-force search" begin
    τ = 0.5
    # Irregular times, as with `saveat` between Lax's tstops
    ts = sort!(vcat(0:τ:10, 0.13:0.71:10))
    snaps = [ Chemostats.Snapshot(t, i, 0, 0.0) for (i, t) in enumerate(ts) ]

    function near_ref(t)
        dist, i = findmin(snap -> abs(snap.t - t), snaps)
        dist <= τ / 2 ? snaps[i] : nothing
    end

    # Decreasing targets, continuing from the previous index
    i = lastindex(snaps)
    for t in [ 11.0; 10.9:-0.37:0.1; -1.0 ]
        snap, i = Chemostats.get_snapshot_near(snaps, t, τ, i)
        @test snap == near_ref(t)
    end
end

@testset "step! rejects unsupported ensemble algorithms" begin
    chem = Chemostat(make_population(4))
    alg = Chemostats.Direct()
    int = Chemostats.PopIntegrator(chem, alg, EnsembleSerial())
    Chemostats.init!(alg, int)

    @test_throws "not supported" Chemostats.step!(int, 5.0, SciMLBase.EnsembleDistributed())
end

@testset "default_ensalg with is_parallel and nthreads" begin
    for alg in (Chemostats.Forward(10), Chemostats.Strict(10), Chemostats.Strict(10, false), Chemostats.Thin(0.1), Chemostats.Lax(10, 1.0))
        ensalg = Chemostats.default_ensalg(alg)
        expect_threaded = Chemostats.is_parallel(alg) && Threads.nthreads() > 1
        @test (ensalg isa EnsembleThreads) == expect_threaded
    end
end
