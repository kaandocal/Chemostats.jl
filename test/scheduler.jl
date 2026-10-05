using Test
using Chemostats

include("models/exponential.jl")

ENSALGS = (EnsembleSerial(), EnsembleThreads())

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

    @test_throws DimensionMismatch Chemostats.simulate!(chem, tmax, Chemostats.Forward(K + 1), ensalg)

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

@testset "Strict population size" begin
    K = 10
    chem = Chemostat(make_population(K))
    tmax = 3.0

    @test_throws "parallelisation" Chemostats.simulate!(chem, tmax, Chemostats.Strict(K), EnsembleThreads())

    Chemostats.simulate!(chem, tmax, Chemostats.Strict(K))
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
    chem = Chemostat(make_population(K); save_leaves = true)

    Chemostats.simulate!(chem, tmax, Chemostats.Thin(δ), ensalg)
    snap = chem.snaps[end]

    @test snap.t == tmax
    @test snap.N == length(chem.pop)

    @test snap.log_f ≈ δ*tmax

    N_dead = length(chem.tree.leaves)
    @test snap.nsim == 2*(snap.N + N_dead) - K
    @test N_dead > 0   # sanity: δ=0.3 over tmax=2 should cull at least a few
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
