using Test
using ArgCheck
using OrdinaryDiffEqTsit5
using Random
using StatsBase
using Chemostats

include("models/exponential.jl")

# Only used for identity
mk() = ExponentialCell((; id = 0))

# The oldest recorded ancestor of `cell` (the founder it descends from), or
# `cell` itself if it has no recorded ancestors.
find_root(tree, cell) = (anc = collect(Chemostats.ancestors(tree, cell)); isempty(anc) ? cell : last(anc))

@testset "PopTree construction" begin
    tree = Chemostats.PopTree{typeof(mk())}()

    @test isempty(tree.parents)
    @test isempty(tree.leaves)
    @test tree.save_ancestors == false
    @test tree.save_leaves == false
    @test tree.last_sweep_size == 0

    tree2 = Chemostats.PopTree{typeof(mk())}(; save_ancestors = true, save_leaves = true)
    @test tree2.save_ancestors
    @test tree2.save_leaves
end

@testset "parent / set_parent!" begin
    tree = Chemostats.PopTree{typeof(mk())}()
    a, b, c = mk(), mk(), mk()

    @test Chemostats.parent(tree, a) === missing

    Chemostats.set_parent!(tree, a, b)
    @test Chemostats.parent(tree, a) == b
    @test tree.parents[a] == b

    @test_throws ArgCheck.CheckError Chemostats.set_parent!(tree, a, c)
    @test tree.parents[a] == b   # unchanged after the rejected attempt
end

@testset "add_offspring! respects save_ancestors" begin
    p = mk()
    c1, c2 = mk(), mk()

    tree_on = Chemostats.PopTree{typeof(mk())}(; save_ancestors = true)
    Chemostats.add_offspring!(tree_on, p, (c1, c2))
    @test tree_on.parents[c1] == p
    @test tree_on.parents[c2] == p

    tree_off = Chemostats.PopTree{typeof(mk())}(; save_ancestors = false)
    Chemostats.add_offspring!(tree_off, p, (c1, c2))
    @test isempty(tree_off.parents)

    tree_single = Chemostats.PopTree{typeof(mk())}(; save_ancestors = true)
    c3 = mk()
    Chemostats.add_offspring!(tree_single, p, (c3,))
    @test tree_single.parents[c3] == p

    tree_empty = Chemostats.PopTree{typeof(mk())}(; save_ancestors = true)
    Chemostats.add_offspring!(tree_empty, p, ())
    @test isempty(tree_empty.parents)
end

@testset "add_clone! inherits parent" begin
    tree = Chemostats.PopTree{typeof(mk())}(; save_ancestors = true)
    root, source = mk(), mk()
    Chemostats.set_parent!(tree, source, root)

    clone = mk()
    Chemostats.add_clone!(tree, source, clone)
    @test Chemostats.parent(tree, clone) == root
    @test Chemostats.parent(tree, source) == root   # source itself is untouched

    # A clone of a root (no recorded parent) stays rootless too.
    tree_root = Chemostats.PopTree{typeof(mk())}(; save_ancestors = true)
    founder, clone_of_founder = mk(), mk()
    Chemostats.add_clone!(tree_root, founder, clone_of_founder)
    @test Chemostats.parent(tree_root, clone_of_founder) === missing

    tree_off = Chemostats.PopTree{typeof(mk())}(; save_ancestors = false)
    Chemostats.add_clone!(tree_off, source, mk())
    @test isempty(tree_off.parents)
end

@testset "add_leaf! respects save_leaves" begin
    cell = mk()

    tree_on = Chemostats.PopTree{typeof(mk())}(; save_leaves = true)
    Chemostats.add_leaf!(tree_on, cell)
    @test tree_on.leaves == [cell]

    tree_off = Chemostats.PopTree{typeof(mk())}(; save_leaves = false)
    Chemostats.add_leaf!(tree_off, cell)
    @test isempty(tree_off.leaves)
end

@testset "ancestors()" begin
    tree = Chemostats.PopTree{typeof(mk())}(; save_ancestors = true)
    root, a, b = mk(), mk(), mk()
    Chemostats.set_parent!(tree, a, root)
    Chemostats.set_parent!(tree, b, a)

    @test collect(Chemostats.ancestors(tree, b)) == [a, root]
    @test collect(Chemostats.ancestors(tree, a)) == [root]
    @test collect(Chemostats.ancestors(tree, root)) == [] 

    @test collect(Chemostats.ancestors(tree, mk())) == []
end

@testset "prune!" begin
    @testset "keeps exactly the ancestry samples need" begin
        tree = Chemostats.PopTree{typeof(mk())}(; save_ancestors = true)

        root = mk()
        a, b = mk(), mk()
        a1, a2 = mk(), mk()
        b1 = mk()           

        Chemostats.set_parent!(tree, a, root)
        Chemostats.set_parent!(tree, b, root)
        Chemostats.set_parent!(tree, a1, a)
        Chemostats.set_parent!(tree, a2, a)
        Chemostats.set_parent!(tree, b1, b)
        @test length(tree.parents) == 5

        Chemostats.prune!(tree, [a1, b1])

        @test Set(keys(tree.parents)) == Set([a1, b1, a, b])
        @test tree.parents[a1] == a
        @test tree.parents[a] == root
        @test tree.parents[b1] == b
        @test tree.parents[b] == root
        @test !haskey(tree.parents, a2)
        @test tree.last_sweep_size == length(tree.parents) == 4

        # ancestors() still walks all the way to the (unkeyed) root.
        @test collect(Chemostats.ancestors(tree, a1)) == [a, root]
    end

    @testset "no-op on root" begin
        tree = Chemostats.PopTree{typeof(mk())}(; save_ancestors = true)
        root, a = mk(), mk()
        Chemostats.set_parent!(tree, a, root)

        Chemostats.prune!(tree, [root])
        @test isempty(tree.parents)
    end

    @testset "empty sample set drops everything" begin
        tree = Chemostats.PopTree{typeof(mk())}(; save_ancestors = true)
        Chemostats.set_parent!(tree, mk(), mk())

        Chemostats.prune!(tree, typeof(mk())[])
        @test isempty(tree.parents)
        @test tree.last_sweep_size == 0
    end

    @testset "is idempotent" begin
        tree = Chemostats.PopTree{typeof(mk())}(; save_ancestors = true)
        root, a, b = mk(), mk(), mk()
        Chemostats.set_parent!(tree, a, root)
        Chemostats.set_parent!(tree, b, root)

        Chemostats.prune!(tree, [a])             # b's branch dropped
        snapshot = copy(tree.parents)
        Chemostats.prune!(tree, [a])             # nothing left to discard
        @test tree.parents == snapshot
    end

    @testset "shared ancestry is only walked once" begin
        # Not a timing test -- just confirms the early-stop-on-shared-path
        # logic doesn't change the *result* vs. walking each sample fully.
        tree = Chemostats.PopTree{typeof(mk())}(; save_ancestors = true)
        root = mk()
        mid = mk()
        Chemostats.set_parent!(tree, mid, root)
        leaves = [mk() for _ in 1:5]
        for l in leaves
            Chemostats.set_parent!(tree, l, mid)
        end

        Chemostats.prune!(tree, leaves)
        @test Set(keys(tree.parents)) == Set(vcat(leaves, [mid]))
        @test all(l -> tree.parents[l] == mid, leaves)
        @test tree.parents[mid] == root
    end
end

@testset "simplify!" begin
    @testset "only prunes once tree has doubled" begin
        tree = Chemostats.PopTree{typeof(mk())}(; save_ancestors = true)
        root, a = mk(), mk()
        Chemostats.set_parent!(tree, a, root)
        tree.last_sweep_size = 10   # pretend a prior prune settled at size 10

        Chemostats.simplify!(tree, [a])   # 1 <= 2*10, must not fire
        @test length(tree.parents) == 1    # untouched -- prune! would replace the dict
        @test tree.last_sweep_size == 10   # untouched -- only prune! updates this

        tree.last_sweep_size = 0   # now 1 > 2*0, must fire
        Chemostats.simplify!(tree, [a])
        @test tree.last_sweep_size == 1    # prune! ran and updated it
    end

    @testset "save_ancestors = false makes this a no-op" begin
        tree = Chemostats.PopTree{typeof(mk())}(; save_ancestors = false)

        never_iterate = Iterators.map(_ -> error("samples should not be touched"), 1:1)
        Chemostats.simplify!(tree, never_iterate)
    end
end

@testset "end-to-end (via simulate!)" begin
    @testset "pruning doesn't affect simulation" begin
        function run_sim(; disable_pruning)
            Random.seed!(20240521)
            chem = Chemostat(make_stochastic_population(10); save_ancestors = true)
            disable_pruning && (chem.tree.last_sweep_size = typemax(Int) ÷ 2)
            Chemostats.simulate!(chem, 5.0, Chemostats.Direct(), EnsembleSerial())
            chem
        end

        chem_pruned = run_sim(; disable_pruning = false)
        chem_full = run_sim(; disable_pruning = true)

        @test chem_pruned.snaps[end].N == chem_full.snaps[end].N
        @test chem_pruned.snaps[end].nsim == chem_full.snaps[end].nsim

        # Pruning must actually have discarded something...
        @test length(chem_pruned.tree.parents) < length(chem_full.tree.parents)

        # ...without breaking any live cell's ancestry: every surviving cell
        # must trace back exactly as far in both trees.
        for (c1, c2) in zip(chem_pruned.pop, chem_full.pop)
            d1 = length(collect(Chemostats.ancestors(chem_pruned.tree, c1)))
            d2 = length(collect(Chemostats.ancestors(chem_full.tree, c2)))
            @test d1 == d2
        end
    end

    @testset "pruning shrinks tree size" begin
        Random.seed!(7)
        L = 100
        chem = Chemostat(make_population(L); save_ancestors = true)
        Chemostats.simulate!(chem, 20.0, Chemostats.Strict(L), EnsembleSerial())
        snap = chem.snaps[end]

        @test snap.N == L
        @test snap.nsim > 20 * L   # make sure plenty of turnover actually happened

        # Without pruning, tree.parents would hold one entry per division
        # ever (≈ nsim); with it, it should stay close to the population size L
        @test length(chem.tree.parents) < snap.nsim ÷ 3
        @test length(chem.tree.parents) < 30 * L
    end
end

@testset "Coalescence" begin
    Λ_gt = 1.0
    N = 100
    tmax = 10 * N / Λ_gt

    @testset "Strict" begin
        Random.seed!(42)
        chem = Chemostat(make_population(N); save_ancestors = true)
        Chemostats.simulate!(chem, tmax, Chemostats.Strict(N))

        @test length(chem.pop) == N
        roots = [ find_root(chem.tree, cell) for cell in chem.pop ]
        @test all(r -> r === roots[1], roots)
    end

    @testset "Lax with $ensalg" for ensalg in (EnsembleSerial(), EnsembleThreads())
        Random.seed!(43)
        chem = Chemostat(make_population(N); save_ancestors = true)
        Chemostats.simulate!(chem, tmax, Chemostats.Lax(N, 0.5), ensalg)

        @test !isempty(chem.pop)   # didn't go extinct
        roots = [ find_root(chem.tree, cell) for cell in chem.pop ]
        @test all(r -> r === roots[1], roots)
    end

    @testset "ancestral ID distribution" begin
        N = 10
        tmax = 10 * N / Λ_gt
        nreps = 100

        winning_ids = Vector{Int}(undef, nreps)
        coalesced = Vector{Bool}(undef, nreps)

        Threads.@threads for i in 1:nreps
            Random.seed!(7000 + i)
            chem = Chemostat(make_population(N))
            Chemostats.simulate!(chem, tmax, Chemostats.Strict(N))

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
end
