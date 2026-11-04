using Test
using Random
using StatsBase
using SciMLBase
using Chemostats

include("models/exponential.jl")
include("models/rngcell.jl")

ENSALGS = (EnsembleSerial(), EnsembleThreads())

node_of(tree, cell) = tree[tree.living[cell]]

# Living nodes reflect their cells' current state, and the tree is ready for analysis
function up_to_date(tree)
    issorted(tree) && all(tree.living) do (cell, id)
        tree[id].t_end == Chemostats.get_curr_t(cell) && tree[id].data === Chemostats.get_data(cell)
    end
end

# Most recent common ancestor of the living cells
function mrca(tree, cells)
    chain = collect(keys(Chemostats.ancestors(tree, first(cells))))
    others = [ Set(keys(Chemostats.ancestors(tree, c))) for c in cells[2:end] ]
    chain[findfirst(id -> all(s -> id in s, others), chain)]
end

@testset "No tree unless requested" begin
    chem = Chemostat(make_population(5))
    Chemostats.simulate!(chem, 1.0, Chemostats.Direct())
    @test isnothing(chem.tree)
end

@testset "Tree records the simulation ($ensalg)" for ensalg in ENSALGS
    K, tmax = 5, 3.0
    founders = make_population(K)
    chem = Chemostat(copy(founders); save_mode = :full)
    Chemostats.simulate!(chem, tmax, Chemostats.Direct(), ensalg)
    tree, snap = chem.tree, chem.snaps[end]

    @test up_to_date(tree)
    @test length(tree) == snap.nsim
    @test length(Chemostats.living(tree)) == length(chem.pop) == snap.N
    @test count(iszero, tree.nodes.parent) == K

    for cell in chem.pop
        node = node_of(tree, cell)
        @test node.t_end == tmax
        @test node.data === cell.int.sol

        lineage = collect(Chemostats.ancestors(tree, cell; self = true))
        @test first(lineage) == node
        @test length(lineage) == node.generation + 1
        @test iszero(last(lineage).parent)
    end

    # No deaths: every finished cell divided in two, at the end of its solution
    for id in eachindex(tree)
        node = tree[id]
        Chemostats.is_alive(node) && continue
        cs = Chemostats.children(tree, id)
        @test length(cs) == 2
        @test all(c -> c.t_birth == node.t_end, cs)
        @test last(node.data.t) == node.t_end
    end

    # `ancestors` with a cell needs a living cell of the tree
    @test_throws ArgumentError Chemostats.ancestors(tree, first(founders))
    @test_throws ArgumentError Chemostats.ancestors(tree, ExponentialCell((; id = 0)))
end

@testset "Killed cells are recorded ($ensalg)" for ensalg in ENSALGS
    tmax = 2.0
    chem = Chemostat(make_population(30); save_mode = :full)
    Chemostats.simulate!(chem, tmax, Chemostats.Thin(0.3), ensalg)
    tree = chem.tree

    killed = [ id for id in eachindex(tree) if !Chemostats.is_alive(tree[id]) && isempty(Chemostats.children(tree, id)) ]
    @test !isempty(killed)
    @test all(id -> tree[id].t_birth <= tree[id].t_end < tmax, killed)
    @test all(id -> !isnothing(tree[id].data), killed)
end

@testset "Forward records single lineages" begin
    K = 6
    chem = Chemostat(make_population(K); save_mode = :full)
    Chemostats.simulate!(chem, 3.0, Chemostats.Forward(K))
    tree = chem.tree

    @test up_to_date(tree)
    @test length(tree) == chem.snaps[end].nsim
    @test all(id -> length(Chemostats.children(tree, id)) <= 1, eachindex(tree))
    @test count(iszero, tree.nodes.parent) == length(Chemostats.living(tree)) == K
end

@testset "Continuing a simulation" begin
    chem = Chemostat(make_population(5); save_mode = :full)
    Chemostats.simulate!(chem, 1.5, Chemostats.Direct())
    Chemostats.simulate!(chem, 3.0, Chemostats.Direct())

    @test up_to_date(chem.tree)
    @test all(cell -> node_of(chem.tree, cell).t_end == 3.0, chem.pop)
    @test length(chem.tree) == chem.snaps[end].nsim
end

@testset "Copied populations are recorded" begin
    cells = make_population(5)
    chem = Chemostat(cells; copy = true, save_mode = :ancestors)
    @test all(c -> !haskey(chem.tree.living, c), cells)

    Chemostats.simulate!(chem, 2.0, Chemostats.Direct())
    @test up_to_date(chem.tree)
    @test all(c -> haskey(chem.tree.living, c), chem.pop)
end

@testset "Tree is up to date after early termination" begin
    @testset "Nmax" begin
        chem = Chemostat(make_population(5); save_mode = :ancestors)
        Chemostats.simulate!(chem, 10.0, Chemostats.Direct(); Nmax = 20)
        @test chem.status == SciMLBase.ReturnCode.MaxIters
        @test up_to_date(chem.tree)
    end

    @testset "Exception ($ensalg)" for ensalg in ENSALGS
        chem = Chemostat(make_throwing_population(5); save_mode = :ancestors)
        @test_throws Exception Chemostats.simulate!(chem, 10.0, Chemostats.Direct(), ensalg; throw_on_error = true)
        @test up_to_date(chem.tree)
        @test !isempty(chem.tree.living)
    end
end

@testset "Cells dropped after errors leave the tree ($ensalg)" for ensalg in ENSALGS
    chem = Chemostat(make_throwing_population(5); save_mode = :ancestors)
    redirect_stderr(devnull) do
        Chemostats.simulate!(chem, 10.0, Chemostats.Direct(), ensalg)
    end

    @test isempty(chem.pop)
    @test isempty(chem.tree.living)
    @test length(chem.tree) == 5
end

@testset "Extinct populations are not pruned" begin
    # Large enough for a lazy sweep, which would otherwise remove every node
    chem = Chemostat(make_dying_population(2000); save_mode = :ancestors)
    Chemostats.simulate!(chem, 20.0, Chemostats.Direct())

    @test isempty(chem.pop)
    @test length(chem.tree) == 2000
end

@testset "Strict cannot record clones" begin
    chem = Chemostat(make_dying_population(5); save_mode = :ancestors)
    @test_throws ArgumentError Chemostats.simulate!(chem, 10.0, Chemostats.Strict(5))
    @test up_to_date(chem.tree)

    # Fine without cloning
    chem = Chemostat(make_dying_population(5); save_mode = :ancestors)
    Chemostats.simulate!(chem, 10.0, Chemostats.Strict(5, false))
    @test up_to_date(chem.tree)
end

@testset "Pruning" begin
    L, tmax = 50, 60.0   # long enough for lazy pruning to sweep several times

    function run_lax(; save_mode)
        Random.seed!(20240521)
        chem = Chemostat(make_population(L); save_mode)
        Chemostats.simulate!(chem, tmax, Chemostats.Lax(L, 0.5), EnsembleSerial(); saveat)
        chem
    end

    saveat = [ 0., 15., 30., 45., 60. ]

    pruned, full = run_lax(; save_mode = :ancestors), run_lax(; save_mode = :full)
    snap = pruned.snaps[end]

    @test [ (s.N, s.nsim) for s in pruned.snaps ] == [ (s.N, s.nsim) for s in full.snaps ]
    @test length(full.tree) == snap.nsim
    @test length(pruned.tree) < snap.nsim ÷ 3
    @test up_to_date(pruned.tree)

    # Populations at `saveat` are marked, so they survive pruning
    for t in saveat
        @test length(Chemostats.alive_at(pruned.tree, t)) == Chemostats._get_snapshot(pruned.snaps, t).N
    end

    # Living lineages are unaffected
    lineage(chem, cell) = [ (node.t_birth, node.t_end) for node in Chemostats.ancestors(chem.tree, cell) ]
    @test all(((c1, c2),) -> lineage(pruned, c1) == lineage(full, c2), zip(pruned.pop, full.pop))
end

@testset "Tree is identical across ensemble algorithms" begin
    function run_thin(ensalg)
        chem = Chemostat([ RNGCell(1000 + i) for i in 1:20 ]; save_mode = :ancestors)
        Chemostats.simulate!(chem, 3.0, Chemostats.Thin(0.3), ensalg; saveat = 0.5:0.5:3.0)
        (; parent, generation, t_birth, t_end, alive) = chem.tree.nodes
        (parent, generation, t_birth, t_end, alive)
    end

    @test run_thin(EnsembleSerial()) == run_thin(EnsembleThreads())
end

@testset "Coalescence" begin
    Λ_gt = 1.0

    # No deaths, so Strict never clones
    @testset "Strict (clone=$clone)" for clone in (true, false)
        N = 100
        Random.seed!(42)
        chem = Chemostat(make_population(N); save_mode = :ancestors)
        Chemostats.simulate!(chem, 10 * N / Λ_gt, Chemostats.Strict(N, clone))

        @test length(chem.pop) == N
        roots = [ last(collect(keys(Chemostats.ancestors(chem.tree, c; self = true)))) for c in chem.pop ]
        @test allequal(roots)
    end

    @testset "Lax with $ensalg" for ensalg in ENSALGS
        N = 100
        Random.seed!(43)
        chem = Chemostat(make_population(N); save_mode = :ancestors)
        Chemostats.simulate!(chem, 10 * N / Λ_gt, Chemostats.Lax(N, 0.5), ensalg)

        @test !isempty(chem.pop)
        roots = [ last(collect(keys(Chemostats.ancestors(chem.tree, c; self = true)))) for c in chem.pop ]
        @test allequal(roots)
    end

    @testset "ancestral ID distribution (clone=$clone)" for clone in (true, false)
        N = 10
        tmax = 10 * N / Λ_gt
        nreps = 100

        winning_ids = Vector{Int}(undef, nreps)
        coalesced = Vector{Bool}(undef, nreps)

        Threads.@threads for i in 1:nreps
            Random.seed!(7000 + i)
            chem = Chemostat(make_population(N))
            Chemostats.simulate!(chem, tmax, Chemostats.Strict(N, clone))

            ids = [ cell_params(cell).id for cell in chem.pop ]
            coalesced[i] = allequal(ids)
            winning_ids[i] = ids[1]
        end

        @test all(coalesced)

        μ = (N + 1) / 2
        σ2 = (N^2 - 1) / 12
        se_mean = sqrt(σ2 / nreps)

        γ2 = -6 * (N^2 + 1) / (5 * (N^2 - 1))   # excess kurtosis, discrete uniform
        μ4 = σ2^2 * (3 + γ2)
        se_var = sqrt((μ4 - σ2^2) / nreps)

        @test abs(mean(winning_ids) - μ) < 5 * se_mean
        @test abs(var(winning_ids) - σ2) < 5 * se_var
    end

    @testset "ancestral ID distribution (Lax)" begin
        N = 20
        tmax = 10 * N / Λ_gt
        nreps = 100

        winning_ids = Vector{Int}(undef, nreps)
        coalesced = Vector{Bool}(undef, nreps)

        Threads.@threads for i in 1:nreps
            Random.seed!(8000 + i)
            chem = Chemostat(make_population(N))
            Chemostats.simulate!(chem, tmax, Chemostats.Lax(N, 1.0; L_min = 10), EnsembleSerial())

            ids = [ cell_params(cell).id for cell in chem.pop ]
            coalesced[i] = !isempty(ids) && allequal(ids)
            winning_ids[i] = isempty(ids) ? 0 : ids[1]
        end

        @test all(coalesced)

        μ = (N + 1) / 2
        σ2 = (N^2 - 1) / 12
        se_mean = sqrt(σ2 / nreps)

        γ2 = -6 * (N^2 + 1) / (5 * (N^2 - 1))   # excess kurtosis, discrete uniform
        μ4 = σ2^2 * (3 + γ2)
        se_var = sqrt((μ4 - σ2^2) / nreps)

        @test abs(mean(winning_ids) - μ) < 5 * se_mean
        @test abs(var(winning_ids) - σ2) < 5 * se_var
    end

    @testset "Time to MRCA (Strict, clone=$clone)" for clone in (true, false)
        N = 10
        tmax = 50 * N / Λ_gt
        nreps = 100

        depths = Vector{Float64}(undef, nreps)
        Threads.@threads for i in 1:nreps
            Random.seed!(13000 + i)
            chem = Chemostat(make_population(N); save_mode = :ancestors)
            Chemostats.simulate!(chem, tmax, Chemostats.Strict(N, clone))
            depths[i] = tmax - chem.tree[mrca(chem.tree, chem.pop)].t_end
        end

        μ = (N^2 - 1) / (N * Λ_gt)
        σ2 = ((N + 1) / Λ_gt)^2 * sum(1 / (k * (k - 1))^2 for k in 2:N)
        se_mean = sqrt(σ2 / nreps)

        @test abs(mean(depths) - μ) < 5 * se_mean
        @test var(depths) < 2 * σ2   # loose sanity check; no kurtosis derivation needed
    end
end
