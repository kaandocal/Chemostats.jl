using Trees
using Trees: cell_init!, cell_divide!, cell_die!, cell_update!
using Aqua
using Random
using StatsBase: mean, var, mean_and_var, pweights
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

    @test keys(roots(tree)) == sort([ n.r1, n.r2 ])
    @test keys(living(tree)) == sort([ n.a2, n.b1, n.b2 ])
    @test is_alive(tree[n.a2]) && !is_alive(tree[n.a1])

    @test parent(tree, n.a1) == n.a && parent(tree, n.r1) == 0
    @test sort(collect(keys(children(tree, n.r1)))) == sort([ n.a, n.b ])
    @test collect(children(tree, n.b)) == [ tree[n.b1], tree[n.b2] ] && first(living(tree)) == tree[n.a2]
    @test isempty(children(tree, n.r2)) && isempty(children(tree, n.a2))
    @test has_divided(tree, n.a) && !has_divided(tree, n.a1)

    @test tree[n.b1].generation == 2
    @test lifetime(tree[n.a2]) == 1. && lifetime(tree[n.r2]) == 1.5
    @test mean(lifetime, ancestors(tree, n.a1; self=true)) ≈ 1.

    @test keys(alive_at(tree, 1.)) == sort([ n.r2, n.a, n.b ])                        # born at 1, not ended at 1
    @test isempty(alive_at(tree, -1.))
end

@testset "Lineages and collections" begin
    tree = example_tree()
    n = node_ids(tree)

    @test collect(keys(ancestors(tree, n.a1))) == [ n.a, n.r1 ]
    @test collect(keys(ancestors(tree, n.a1; self=true))) == [ n.a1, n.a, n.r1 ]
    @test collect(keys(ancestors(tree, n.a1; self=true, tspan=(1.5, 2.5)))) == [ n.a1, n.a ]
    @test isempty(ancestors(tree, n.r1))

    @test collect(keys(descendants(tree, n.a))) == sort([ n.a1, n.a2 ])
    @test first(descendants(tree, n.r1; self=true)) == tree[n.r1]
    @test sort(collect(keys(descendants(tree, n.r1; tspan=(0., 2.2))))) == sort([ n.a, n.b, n.a1, n.a2 ])
    @test sort(collect(keys(descendants(tree, n.r1; tspan=(2.2, 3.))))) == sort([ n.b, n.a1, n.a2, n.b1, n.b2 ])      # skips a

    late = filter(node -> node.t_birth >= 2, tree)
    @test sort(keys(late)) == sort([ n.a1, n.a2, n.b1, n.b2 ])
    @test keys(alive_at(late, 2.2)) == sort([ n.a1, n.a2 ])
    @test keys(alive_at(ancestors(tree, n.b1; self=true), 1.2)) == [ n.b ]
    @test keys(filter(node -> node.t_birth >= 1, ancestors(tree, n.a1; self=true))) == [ n.a1, n.a ]    # keeps the order

    # Collections iterate over nodes, `keys` gives their ids; trees themselves are not iterable
    @test keys(tree) == eachindex(tree) == 1:8
    @test collect(values(tree)) == [ tree[i] for i in 1:8 ]
    @test collect(pairs(tree)) == [ i => tree[i] for i in 1:8 ]
    @test_throws MethodError iterate(tree)
    @test !isempty(tree) && isempty(Tree{Cell}())
    @test collect(view(tree, [ n.a, n.b ])) == [ tree[n.a], tree[n.b] ]
    @test collect(ancestors(tree, n.a1)) == [ tree[n.a], tree[n.r1] ]
    @test collect(pairs(ancestors(tree, n.a1))) == [ n.a => tree[n.a], n.r1 => tree[n.r1] ]
    @test eltype(collect(descendants(tree, n.a))) == Trees.Node{Float64}

    @test sprint(show, living(tree)) == "TreeView(3 nodes)"
    @test sprint(show, tree[n.a1]) == "Node{Float64}(2.0 → 3.0, data = 2.0)"
end

@testset "Ancestors of a living cell" begin
    r1, r2, a, b, a1 = Cell(:r1), Cell(:r2), Cell(:a), Cell(:b), Cell(:a1)
    tree = Tree{Cell,Float64}()
    cell_init!(tree, [ r1, r2 ], 0.)
    cell_die!(tree, r2, 0.5, 0.)
    cell_divide!(tree, r1, [ a, b ], 1., 0.)
    cell_divide!(tree, a, [ a1 ], 2., 1.)

    by_id(cell; kwargs...) = collect(keys(ancestors(tree, tree.living[cell]; kwargs...)))
    @test collect(keys(ancestors(tree, a1))) == by_id(a1)
    @test collect(keys(ancestors(tree, a1; self=true, tspan=(1.5, 3.)))) == by_id(a1; self=true, tspan=(1.5, 3.))

    @test_throws ArgumentError ancestors(tree, a)            # divided
    @test_throws ArgumentError ancestors(tree, r2)           # died
    @test_throws ArgumentError ancestors(tree, Cell(:x))     # never in the tree

    # Pruning renumbers nodes; cells still find their lineages
    prune!(tree)
    @test length(tree) == 4
    @test collect(keys(ancestors(tree, a1))) == by_id(a1)
    @test [ node.t_birth for node in ancestors(tree, a1; self=true) ] == [ 2., 1., 0. ]
    @test collect(keys(ancestors(tree, b))) == keys(roots(tree))
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
    @test time_weights.(tree.nodes) ≈ w
    @test mean(tb, tree, time_weights) ≈ sum(w .* vals) / sum(w)
    @test all(mean_and_var(tb, tree, time_weights) .≈ (sum(w .* vals) / sum(w), var(vals, pweights(w); corrected=true)))
    @test time_weights(tree[n.r1]; tspan=(0.5, 1.2)) ≈ time_weights(; tspan=(0.5, 1.2))(tree[n.r1]) ≈ 0.5

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
    mean(tb, chain, time_weights)
    @test (@allocated mean(tb, chain, time_weights)) < 1000

    # Statistics over lineages and views
    @test mean(tb, ancestors(tree, n.b1; self=true)) ≈ mean([ 2.5, 1., 0. ])
    @test mean(tb, filter(node -> node.t_birth >= 2, tree)) ≈ 2.25

    age(node, a) = a
    @test time_average(age, tree[n.b]) ≈ time_average(age)(tree[n.b]) ≈ 0.75
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
    @test snapshot(abs_t, tree, 2.7) ≈ [ 2.7 ]                                    # only a1

    for (cell, id) in tree.living
        cell_update!(tree, cell, 3., tree[id].t_birth)
    end

    snap = snapshot(abs_t, tree, 2.7)
    @test length(snap) == length(alive_at(tree, 2.7)) == 4 && all(snap .≈ 2.7)

    @test snapshot(abs_t, tree, [ 0.5, 2.7 ]) ≈ [ 0.5, 0.5, 2.7, 2.7, 2.7, 2.7 ]
    @test isempty(snapshot(abs_t, tree, Float64[]))

    # On a lineage, a snapshot is the single lineage cell alive at each time
    lineage = ancestors(tree, n.a2; self=true)                                    # a2, a, r1
    @test snapshot(abs_t, lineage, [ 0.5, 1.5, 2.5 ]) ≈ [ 0.5, 1.5, 2.5 ]
    @test length(at_division(abs_t, lineage)) == 2                                # a and r1

    # Ages: living cells only count up to their age at their last update
    rec = at_age(abs_t, tree, 0.8)
    @test sort(rec) ≈ [ 0., 0., 1., 1., 2., 2. ] .+ 0.8                        # r1, r2, a, b, a1, a2
    @test length(at_birth(abs_t, tree)) == 8

    divs = at_division((node, τ) -> (abs_t(node, 0), τ), tree)               # birth time and generation time
    @test sort(last.(divs)) ≈ [ 1., 1., 1.5 ]
    @test sort(first.(divs) .+ last.(divs)) ≈ [ 1., 2., 2.5 ]
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
    @test_throws ArgumentError cell_divide!(tree, cell, [ cell, d ], 4., 0.)           # its own daughter
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

    # Divisions with zero lifetime
    tree = Tree{Cell,Float64}()
    c, d, e, f, g = (Cell(x) for x in (:c, :d, :e, :f, :g))
    cell_init!(tree, [ c, g ], 0.)
    cell_divide!(tree, c, [ d, e ], 0., 0.)
    cell_divide!(tree, d, [ f ], 1., 0.)
    sort!(tree)
    @test collect(keys(children(tree, 1))) == [ 3, 4 ] && collect(keys(descendants(tree, 1))) == [ 3, 5, 4 ]
    @test isnothing(time_average((node, a) -> a, tree[1])) && 1 ∉ keys(alive_at(tree, 0.))

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
    @test keys(living(tree)) == [ tree.living[c[5]] ]

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
    @test length(collect(ancestors(tree, only(keys(living(tree)))))) == 50

    # Marked node collections are kept
    tree = example_tree()
    mark!(tree, alive_at(tree, 1.2))
    @test length(alive_at(prune!(tree), 1.2)) == 3                                 # keeps r2, removed at 1.5
end

@testset "Newick" begin
    tree = example_tree()
    n = node_ids(tree)

    @test newick(tree) == "(($(n.a1):1.0,$(n.a2):1.0)$(n.a):1.0,($(n.b1):0.5,$(n.b2):0.5)$(n.b):1.5)$(n.r1):1.0;\n" *
                          "$(n.r2):1.5;\n"
    @test sprint(newick, tree) == newick(tree)
    @test newick(descendants(tree, n.b; self=true)) == "($(n.b1):0.5,$(n.b2):0.5)$(n.b):1.5;\n"
    @test newick(ancestors(tree, n.a1; self=true)) == "(($(n.a1):1.0)$(n.a):1.0)$(n.r1):1.0;\n"
    @test newick(view(tree, [ n.r1, n.a1, n.b ])) == "($(n.b):1.5)$(n.r1):1.0;\n$(n.a1):1.0;\n"     # a1 without its parent
    @test newick(view(tree, Int[])) == ""

    # Polytomies, with a middle child missing from the collection
    poly = Tree{Cell,Float64}()
    cell, kids = Cell(:c), [ Cell(:k) for _ in 1:3 ]
    cell_init!(poly, [ cell ], 0.)
    cell_divide!(poly, cell, kids, 1., 0.)
    foreach(kid -> cell_update!(poly, kid, 2.), kids)
    sort!(poly)
    @test newick(poly) == "(2:1.0,3:1.0,4:1.0)1:1.0;\n"
    @test newick(view(poly, [ 1, 2, 4 ])) == "(2:1.0,4:1.0)1:1.0;\n"

    # Labels are quoted if needed, and annotations are written after labels
    r2 = view(tree, [ n.r2 ])
    @test newick(r2; label = node -> "it's (x)", annotate = node -> (t = node.t_birth, alive = node.alive, s = :x)) ==
          "'it''s (x)'[&t=0.0,alive=false,s=\"x\"]:1.5;\n"
    @test newick(r2; label = Returns(nothing), annotate = Returns(nothing)) == ":1.5;\n"
    @test newick(r2; label = node -> "x_1", annotate = Returns((;))) == "x_1:1.5;\n"

    # Deep lineages do not overflow the stack
    chain = Tree{Cell,Float64}()
    cell = Cell(:x)
    cell_init!(chain, [ cell ], 0.)
    for k in 1:100_000
        kid = Cell(:y)
        cell_divide!(chain, cell, [ kid ], Float64(k), 0.)
        cell = kid
    end
    @test_throws ArgumentError newick(chain)                                        # unsorted
    s = newick(sort!(chain))
    @test count(==('('), s) == 100_000 && endswith(s, ")1:1.0;\n")
end

include("random_trees.jl")

@testset "Aqua" begin
    Aqua.test_all(Trees)
end

end
