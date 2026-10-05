using Test
using Chemostats

include("models/exponential.jl")

ENSALGS = (EnsembleSerial(), EnsembleThreads())

# NOTE: EnsembleThreads() exercises the parallel code path correctly even
# when Threads.nthreads() == 1 (the default for `julia test/runtests.jl`),
# but only genuinely stresses the threading/scheduling logic -- the thing
# `resilience.jl` mainly exists to protect -- when run with e.g. `julia -t auto`.

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

    # No cell ever dies in this model, so every starting lineage must still
    # have at least one descendant present.
    ids = Set(cell_params(cell).id for cell in chem.pop)
    @test ids == Set(1:K)
end

@testset "Forward(L) preserves lineage identity" for ensalg in ENSALGS
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

    # Forward replaces a dividing cell with one daughter, which inherits the
    # parent's parameters unchanged -- so each starting lineage's `id` must
    # still be present, exactly once, at the end.
    ids = sort([ cell_params(cell).id for cell in chem.pop ])
    @test ids == 1:K
end

@testset "Strict(L) keeps population size exactly constant" begin
    K = 10
    chem = Chemostat(make_population(K))
    tmax = 3.0

    @test_throws "parallelisation" Chemostats.simulate!(chem, tmax, Chemostats.Strict(K), EnsembleThreads())

    Chemostats.simulate!(chem, tmax, Chemostats.Strict(K))
    snap = chem.snaps[end]

    @test snap.t == tmax
    @test snap.N == K
    @test length(chem.pop) == K

    # Strict clones existing cells (inheriting their `id`) and culls others --
    # no cell ever introduces a new `id`, so every surviving `id` must still
    # be one of the K originals (though some originals may now be missing,
    # and others duplicated, as lineages drift under resampling).
    @test all(cell -> cell_params(cell).id in 1:K, chem.pop)
end

@testset "Lax(L, τ) keeps population roughly bounded" for ensalg in ENSALGS
    L = 20
    chem = Chemostat(make_population(L))
    tmax = 5.0
    alg = Chemostats.Lax(L, 0.5)

    Chemostats.simulate!(chem, tmax, alg, ensalg; Nmax = 10_000)
    snap = chem.snaps[end]

    @test snap.t == tmax
    # Loose sanity bound, not a tight one: Lax only resizes *between* rounds,
    # so a burst of growth within the final round can overshoot `L*tol`
    # before the next check would have caught it. This just guards against
    # gross regressions (unbounded blow-up or extinction), not exact tuning.
    @test 0 < snap.N < 20*L
    @test length(chem.pop) == snap.N

    # Same reasoning as Strict: cloning never introduces a new `id`.
    @test all(cell -> cell_params(cell).id in 1:L, chem.pop)
end

@testset "Thin(δ) culls cells but keeps bookkeeping exact" for ensalg in ENSALGS
    K = 30
    δ = 0.3
    tmax = 2.0
    chem = Chemostat(make_population(K); save_leaves = true)

    Chemostats.simulate!(chem, tmax, Chemostats.Thin(δ), ensalg)
    snap = chem.snaps[end]

    @test snap.t == tmax
    @test snap.N == length(chem.pop)

    # log_f accumulates `δ * (round duration)` every round regardless of what
    # happened to any individual cell, so it must equal δ*tmax exactly
    # (up to floating-point error) -- this is what est_N(snap) actually
    # relies on to correct for the culling, so it being right matters.
    @test snap.log_f ≈ δ*tmax

    # Extend the full-binary-forest identity from the Direct() test to
    # account for cells killed by δ-culling (recorded as leaves since
    # save_leaves=true): nsim = 2*(N_alive + N_dead) - K.
    N_dead = length(chem.tree.leaves)
    @test snap.nsim == 2*(snap.N + N_dead) - K
    @test N_dead > 0   # sanity: δ=0.3 over tmax=2 should cull at least a few
end

@testset "Simulation can be resumed across multiple simulate! calls" begin
    K = 4
    chem = Chemostat(make_population(K))

    Chemostats.simulate!(chem, 1.0, Chemostats.Direct())
    @test chem.snaps[end].t == 1.0

    Chemostats.simulate!(chem, 2.0, Chemostats.Direct())
    @test chem.snaps[end].t == 2.0
    @test length(chem.snaps) >= 3
end

@testset "saveat records intermediate snapshots" begin
    K = 4
    chem = Chemostat(make_population(K))
    tmax = 3.0

    Chemostats.simulate!(chem, tmax, Chemostats.Direct(); saveat = [1.0, 2.0, 3.0])
    ts = [ s.t for s in chem.snaps ]

    @test issorted(ts)
    @test all(t -> t in ts, (1.0, 2.0, 3.0))
end
