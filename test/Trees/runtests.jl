using Trees
using Aqua
using Random
using StatsBase: mean, var, pweights
using Test

mutable struct Cell
    name::Symbol
end

# Hand-built tree, with the birth time of each cell as its data:
#   r1 (0-1) -> a (1-2) -> a1 (2-3, removed), a2 (2-, alive)
#            -> b (1-2.5) -> b1, b2 (2.5-, alive)
#   r2 (0-1.5, removed)
function example_tree()
    cells = Dict(name => Cell(name) for name in (:r1, :r2, :a, :b, :a1, :a2, :b1, :b2))
    tree = Tree{Cell,Float64}()

    cell_init!(tree, [ cells[:r1], cells[:r2] ], 0.)
    cell_divide!(tree, cells[:r1], [ cells[:a], cells[:b] ], 1., 0.)
    cell_die!(tree, cells[:r2], 1.5, 0.)
    cell_divide!(tree, cells[:a], [ cells[:a1], cells[:a2] ], 2., 1.)
    cell_divide!(tree, cells[:b], [ cells[:b1], cells[:b2] ], 2.5, 1.)
    cell_die!(tree, cells[:a1], 3., 2.)
    foreach(name -> cell_update!(tree, cells[name], 3.), (:a2, :b1, :b2))

    sort!(tree)
end

# Node id of each named cell in a sorted example tree, identified by birth and end time
function node_ids(tree)
    find(tb, te; alive=false) = only(i for i in eachindex(tree) if tree[i].t_birth == tb && tree[i].t_end == te &&
                                                                    tree[i].alive == alive)
    r1, r2 = find(0., 1.), find(0., 1.5)
    a, b = find(1., 2.), find(1., 2.5)
    a1, a2 = find(2., 3.), find(2., 3.; alive=true)
    b1, b2 = filter(i -> tree[i].t_birth == 2.5, eachindex(tree))
    (; r1, r2, a, b, a1, a2, b1, b2)
end

@testset "Trees" begin

@testset "Structure" begin
    tree = example_tree()
    n = node_ids(tree)

    @test length(tree) == 8

    @test roots(tree) == sort([ n.r1, n.r2 ])
    @test living(tree) == sort([ n.a2, n.b1, n.b2 ])
    @test tree[n.a2].alive && !tree[n.a1].alive
    @test is_alive(tree[n.a2]) && !is_alive(tree[n.a1])

    @test parent(tree, n.a1) == n.a && parent(tree, n.r1) == 0
    @test sort(collect(children(tree, n.r1))) == sort([ n.a, n.b ])
    @test sort(collect(children(tree, n.b))) == sort([ n.b1, n.b2 ])
    @test isempty(children(tree, n.r2)) && isempty(children(tree, n.a2))
    @test has_divided(tree, n.a) && !has_divided(tree, n.a1)

    @test tree[n.b1].generation == 2
    @test lifetime(tree[n.a2]) == 1. && lifetime(tree[n.r2]) == 1.5
    @test mean(lifetime, ancestors(tree, n.a1; self=true)) ≈ 1.

    @test alive_at(tree, 0.5) == sort([ n.r1, n.r2 ])
    @test alive_at(tree, 1.) == sort([ n.r2, n.a, n.b ])
    @test alive_at(tree, 2.7) == sort([ n.a1, n.a2, n.b1, n.b2 ])
    @test isempty(alive_at(tree, -1.))
end

@testset "Lineages and collections" begin
    tree = example_tree()
    n = node_ids(tree)

    @test collect(ancestors(tree, n.a1)) == [ n.a, n.r1 ]
    @test collect(ancestors(tree, n.a1; self=true)) == [ n.a1, n.a, n.r1 ]
    @test collect(ancestors(tree, n.a1; self=true, tspan=(1.5, 2.5))) == [ n.a1, n.a ]
    @test isempty(ancestors(tree, n.r1))

    @test collect(descendants(tree, n.a)) == sort([ n.a1, n.a2 ])
    @test length(collect(descendants(tree, n.r1))) == 6
    @test first(descendants(tree, n.r1; self=true)) == n.r1
    @test sort(collect(descendants(tree, n.r1; tspan=(0., 2.2)))) == sort([ n.a, n.b, n.a1, n.a2 ])
    @test sort(collect(descendants(tree, n.r1; tspan=(2.2, 3.)))) == sort([ n.b, n.a1, n.a2, n.b1, n.b2 ])      # skips a

    late = filter(node -> node.t_birth >= 2, tree)
    @test sort(late.ids) == sort([ n.a1, n.a2, n.b1, n.b2 ])
    @test alive_at(late, 2.2) == sort([ n.a1, n.a2 ])
    @test alive_at(ancestors(tree, n.b1; self=true), 1.2) == [ n.b ]

    @test ids(tree) == eachindex(tree)
    @test ids(view(tree, [ n.a, n.b ])) == [ n.a, n.b ]
    @test collect(ids(ancestors(tree, n.a1))) == collect(ancestors(tree, n.a1))
    @test collect(eachnode(tree)) == [ tree[i] for i in eachindex(tree) ]
    @test collect(eachnode(ancestors(tree, n.a1))) == [ tree[n.a], tree[n.r1] ]
end

@testset "Ancestors of a living cell" begin
    r1, r2, a, b, a1 = Cell(:r1), Cell(:r2), Cell(:a), Cell(:b), Cell(:a1)
    tree = Tree{Cell,Float64}()
    cell_init!(tree, [ r1, r2 ], 0.)
    cell_die!(tree, r2, 0.5, 0.)
    cell_divide!(tree, r1, [ a, b ], 1., 0.)
    cell_divide!(tree, a, [ a1 ], 2., 1.)

    by_id(cell; kwargs...) = collect(ancestors(tree, tree.living[cell]; kwargs...))
    @test collect(ancestors(tree, a1)) == by_id(a1)
    @test collect(ancestors(tree, a1; self=true, tspan=(1.5, 3.))) == by_id(a1; self=true, tspan=(1.5, 3.))
    @test length(collect(ancestors(tree, a1))) == 2

    @test_throws ArgumentError ancestors(tree, a)            # divided
    @test_throws ArgumentError ancestors(tree, r2)           # died
    @test_throws ArgumentError ancestors(tree, Cell(:x))     # never in the tree

    # Pruning renumbers nodes; cells still find their lineages
    prune!(tree)
    @test length(tree) == 4
    @test collect(ancestors(tree, a1)) == by_id(a1)
    @test [ tree[i].t_birth for i in ancestors(tree, a1; self=true) ] == [ 2., 1., 0. ]
    @test collect(ancestors(tree, b)) == roots(tree)
end

@testset "Statistics" begin
    tree = example_tree()
    n = node_ids(tree)
    tb(node) = node.t_birth

    vals = tree.nodes.t_birth
    @test mean(tb, tree) ≈ mean(vals)
    @test var(tb, tree) ≈ var(vals)
    @test all(mean_and_var(tb, tree) .≈ (mean(vals), var(vals)))

    # Time weights: lifetimes (living cells up to their last update), clipped to `tspan`
    w = [ tree[i].t_end - tree[i].t_birth for i in eachindex(tree) ]
    @test time_weight.(tree.nodes) ≈ w
    @test mean(tb, tree, time_weight) ≈ sum(w .* vals) / sum(w)
    @test all(mean_and_var(tb, tree, time_weight) .≈ (sum(w .* vals) / sum(w), var(vals, pweights(w); corrected=true)))
    @test time_weight(tree[n.r1]; tspan=(0.5, 1.2)) ≈ time_weight(; tspan=(0.5, 1.2))(tree[n.r1]) ≈ 0.5

    # Zero weights and `nothing` values are skipped
    @test mean(tb, tree, node -> node.t_end == 2.5 ? 1. : 0.) == 1.                  # only b
    @test mean(node -> node.data, tree) ≈ mean(filter(!isnothing, tree.nodes.data))
    @test isnan(mean(node -> nothing, tree))
    @test isnan(mean(tb, tree, node -> 0.))
    @test isnan(last(mean_and_var(tb, view(tree, [ n.r1 ]))))                         # a single value

    # Vector-valued summaries
    m, v = mean_and_var(node -> [ node.t_birth, node.t_end ], tree)
    @test m ≈ [ mean(vals), mean(tree.nodes.t_end) ] && v ≈ [ var(vals), var(tree.nodes.t_end) ]

    # Weighted means do not allocate per node
    chain = Tree{Cell,Float64}()
    cell = Cell(:x)
    cell_init!(chain, [ cell ], 0.)
    for k in 1:1000
        kid = Cell(:y)
        cell_divide!(chain, cell, [ kid ], Float64(k), 0.)
        cell = kid
    end
    sort!(chain)
    mean(tb, chain, time_weight)
    @test (@allocated mean(tb, chain, time_weight)) < 1000

    # Statistics over lineages and views
    @test mean(tb, ancestors(tree, n.b1; self=true)) ≈ mean([ 2.5, 1., 0. ])
    @test mean(tb, filter(node -> node.t_birth >= 2, tree)) ≈ 2.25

    @test Trees.quadgk_integral((node, a) -> a^2, nothing, 0., 2.) ≈ 8 / 3
    age(node, a) = a
    @test time_average(age, tree[n.b]) ≈ time_average(age)(tree[n.b]) ≈ 0.75
    @test mean(time_average(age), tree, time_weight) ≈ mean(node -> time_average(age, node), tree, time_weight)
    @test time_average(age, tree[n.b]; tspan=(1.5, 2.5)) ≈ 1.                       # ages 0.5 to 1.5
    @test isnothing(time_average(age, tree[n.b]; tspan=(3., 4.)))
    @test time_average(age, tree[n.b]; integrate=(f, node, a, b) -> 3. * (b - a)) ≈ 3.

    # Living cells: no data yet, then averaged up to their last update or the end of `tspan`
    @test isnothing(time_average(age, tree[n.a2]; tspan=(0., 3.)))
    cell_update!(tree, only(c for (c, i) in tree.living if i == n.a2), 3., 2.)
    @test time_average(age, tree[n.a2]; tspan=(0., 3.)) ≈ 0.5
    @test time_average(age, tree[n.a2]) ≈ 0.5                                        # up to the last update
end

@testset "Snapshots and records" begin
    tree = example_tree()
    n = node_ids(tree)
    abs_t(node, a) = node.data + a      # data is the birth time, so this is absolute time

    # Living cells have no data yet and are skipped
    snap = snapshot(abs_t, tree, 2.7)
    @test snap.id == [ n.a1 ] && snap.value ≈ [ 2.7 ]

    for (cell, id) in tree.living
        cell_update!(tree, cell, 3., tree[id].t_birth)
    end

    snap = snapshot(abs_t, tree, 2.7)
    @test snap.id == alive_at(tree, 2.7) && all(snap.value .≈ 2.7)

    snaps = snapshot(abs_t, tree, [ 0.5, 2.7 ])
    @test snaps.t == [ 0.5, 0.5, 2.7, 2.7, 2.7, 2.7 ] && snaps.value ≈ snaps.t
    @test isempty(snapshot(abs_t, tree, Float64[]).id)

    # Ages: living cells only count up to their age at their last update
    rec = at_age(abs_t, tree, 0.8)
    @test sort(rec.id) == sort([ n.r1, n.r2, n.a, n.b, n.a1, n.a2 ])
    @test all(rec.value .≈ getindex.(Ref(tree.nodes.t_birth), rec.id) .+ 0.8)
    @test length(at_birth(abs_t, tree).id) == 8

    pairs = at_division((node, τ) -> (abs_t(node, 0), τ), tree).value        # birth time and generation time
    @test sort(last.(pairs)) ≈ [ 1., 1., 1.5 ]
    @test sort(first.(pairs) .+ last.(pairs)) ≈ [ 1., 2., 2.5 ]
end

@testset "Updates, sorting and errors" begin
    tree = example_tree()
    @test_throws ArgumentError cell_die!(tree, Cell(:x), 4., 0.)
    @test_throws ArgumentError cell_init!(tree, [ first(keys(tree.living)) ], 4.)
    @test_throws ArgumentError cell_update!(tree, Cell(:x), 4.)
    @test_throws ArgumentError cell_divide!(tree, Cell(:x), [ Cell(:y) ], 4., 0.)

    # Invalid divisions and times leave the tree unchanged
    cell, other = collect(keys(tree.living))[1:2]
    nodes, n_living = deepcopy(tree.nodes), length(tree.living)
    d = Cell(:d)
    @test_throws ArgumentError cell_divide!(tree, cell, [ d, d ], 4., 0.)
    @test_throws ArgumentError cell_divide!(tree, cell, [ other ], 4., 0.)
    @test_throws ArgumentError cell_divide!(tree, cell, [ d ], 0.5, 0.)               # before its birth
    @test_throws ArgumentError cell_die!(tree, cell, 0.5, 0.)
    @test_throws ArgumentError cell_update!(tree, cell, 0.5)
    @test tree.nodes == nodes && length(tree.living) == n_living

    # Updates without data only change the time
    cell, id = first(tree.living)
    cell_update!(tree, cell, 3.5, 7.)
    cell_update!(tree, cell, 4.)
    @test tree[id].t_end == 4. && tree[id].data == 7. && tree[id].alive

    # Functions needing time order do not sort the tree themselves
    cell_divide!(tree, cell, [ Cell(:c1), Cell(:c2) ], 4.5, 3.)
    @test !issorted(tree)
    @test_throws ArgumentError children(tree, 1)
    @test_throws ArgumentError descendants(tree, 1)
    @test_throws ArgumentError alive_at(tree, 4.)
    @test_throws ArgumentError snapshot((node, a) -> a, tree, 4.)
    @test_throws ArgumentError at_division((node, a) -> a, tree)
    sort!(tree)
    @test issorted(tree) && length(alive_at(tree, 4.6)) == 4

    # Lazy pruning only sweeps large trees
    n = length(tree)
    prune!(tree; lazy=true)
    @test length(tree) == n
    prune!(tree)
    @test length(tree) < n

    # A mother can be reused as one of its daughters
    tree = Tree{Cell,Float64}()
    c, d = Cell(:c), Cell(:d)
    cell_init!(tree, [ c ], 0.)
    cell_divide!(tree, c, [ c, d ], 1., 0.)
    sort!(tree)
    @test living(tree) == [ 2, 3 ] && parent(tree, tree.living[c]) == 1

    # Divisions with zero lifetime
    tree = Tree{Cell,Float64}()
    c, d, e, f, g = (Cell(x) for x in (:c, :d, :e, :f, :g))
    cell_init!(tree, [ c, g ], 0.)
    cell_divide!(tree, c, [ d, e ], 0., 0.)
    cell_divide!(tree, d, [ f ], 1., 0.)
    sort!(tree)
    @test collect(children(tree, 1)) == [ 3, 4 ] && collect(descendants(tree, 1)) == [ 3, 5, 4 ]
    @test isnothing(time_average((node, a) -> a, tree[1])) && 1 ∉ alive_at(tree, 0.)

    # Data type defaults to `Any`
    @test Tree{Cell}() isa Tree{Cell,Any}
end

@testset "Cells removed at birth" begin
    # The tree is the same as if the removed cells had never been added
    c = [ Cell(Symbol(:c, i)) for i in 1:9 ]
    tree, ref = Tree{Cell,Float64}(), Tree{Cell,Float64}()
    cell_init!(tree, c[1:3], 0.)
    cell_init!(ref, c[2:3], 0.)
    foreach(tr -> cell_divide!(tr, c[2], c[4:5], 1., 0.), (tree, ref))

    cell_die!(tree, c[1], 0., 0.)                       # a founder, after later events
    cell_divide!(tree, c[3], c[6:7], 2., 0.)
    cell_die!(tree, c[6], 2., 0.)                       # not the last node
    cell_divide!(ref, c[3], c[7:7], 2., 0.)

    cell_divide!(tree, c[4], c[8:9], 3., 1.)
    foreach(x -> cell_die!(tree, x, 3., 0.), c[8:9])    # leaves the mother without children
    cell_die!(ref, c[4], 3., 1.)

    @test tree.nodes == ref.nodes && tree.living == ref.living
    @test sort!(tree).nodes == sort!(ref).nodes && tree.living == ref.living

    # A sorted tree stays sorted
    cell_die!(tree, c[7], 2., 0.)
    @test issorted(tree) && length(tree) == 4 && isempty(children(tree, 2))     # c3
    @test living(tree) == [ tree.living[c[5]] ]

    # Pruned trees, with frozen nodes
    function lineage(drop)
        tree = Tree{Cell,Float64}(; min_sweep=1)
        cell = Cell(:x)
        cell_init!(tree, [ cell ], 0.)
        for k in 1:20
            kids = drop ? [ Cell(:y), Cell(:z) ] : [ Cell(:z) ]
            cell_divide!(tree, cell, kids, Float64(k), 0.)
            drop && cell_die!(tree, kids[1], Float64(k), 0.)
            cell = last(kids)
            prune!(tree; lazy=true)
        end
        tree
    end

    dropped, ref = lineage(true), lineage(false)
    @test dropped.nodes == ref.nodes && Trees.n_frozen(dropped) == Trees.n_frozen(ref) > 0
end

@testset "Pruning edge cases" begin
    # An extinct population without marks is removed entirely
    tree = Tree{Cell,Float64}()
    c, d = Cell(:c), Cell(:d)
    cell_init!(tree, [ c ], 0.)
    cell_divide!(tree, c, [ d ], 1., 0.)
    cell_die!(tree, d, 2., 0.)
    @test length(prune!(tree)) == 0

    # Without a common ancestor of the living cells, nothing is frozen
    tree = Tree{Cell,Float64}()
    cells = [ Cell(:x), Cell(:y) ]
    cell_init!(tree, cells, 0.)
    for k in 1:5, i in 1:2
        kid = Cell(:z)
        cell_divide!(tree, cells[i], [ kid ], Float64(k), 0.)
        cells[i] = kid
    end
    prune!(tree)
    @test length(tree) == 12 && Trees.n_frozen(tree) == 0

    # A single forward lineage keeps its ancestors, which are frozen
    tree = Tree{Cell,Float64}(; min_sweep=1)
    cell = Cell(:x)
    cell_init!(tree, [ cell ], 0.)
    for k in 1:50
        kids = [ Cell(:y), Cell(:z) ]
        cell_divide!(tree, cell, kids, Float64(k), 0.)
        cell_die!(tree, kids[2], k + 0.5, 0.)
        cell = kids[1]
        prune!(tree; lazy=true)
    end
    prune!(tree)
    @test length(tree) == 51 && Trees.n_frozen(tree) == 50
    @test length(collect(ancestors(tree, only(living(tree))))) == 50
end

include("random_trees.jl")

@testset "Aqua" begin
    Aqua.test_all(Trees)
end

end
