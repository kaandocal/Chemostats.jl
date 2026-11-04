# Properties of arbitrary random trees, checked against brute force. At each step, every living cell divides into
# 1-3 children, is removed, or survives; in half of the steps, all events happen at exactly the same time.
# Events within a step involve distinct cells, so they may arrive in any order. The data of a cell is its birth time.

function random_events(rng; n0=4, steps=30, nmax=40)
    cells = [ Cell(:x) for _ in 1:n0 ]
    birth = IdDict{Cell,Float64}(c => 0. for c in cells)
    parent = IdDict{Cell,Union{Cell,Nothing}}(c => nothing for c in cells)
    events = Any[ (:init, copy(cells), 0.) ]
    t_mark = rand(rng, 5:steps-5)

    for k in 1:steps
        batch = []
        ties = rand(rng, Bool)
        for c in copy(cells)
            u = rand(rng)
            t = ties ? Float64(k) : k - rand(rng) / 2
            if u < 0.3 && length(cells) < nmax
                kids = [ Cell(:x) for _ in 1:rand(rng, 1:3) ]
                push!(batch, (:divide, c, kids, t, birth[c]))
                filter!(x -> x !== c, cells)
                append!(cells, kids)
                foreach(kid -> (birth[kid] = t; parent[kid] = c), kids)
            elseif u < 0.45 && length(cells) > 2
                push!(batch, (:remove, c, t, birth[c]))
                filter!(x -> x !== c, cells)
            end
        end

        push!(events, batch, (:step, Float64(k)))
        k == t_mark && push!(events, (:mark,))
    end

    push!(events, (:update, copy(cells), Float64(steps)))

    events, (; cells, birth, parent, t_mark)
end

function apply!(tree, e)
    e[1] == :init && return cell_init!(tree, e[2], e[3])
    e[1] == :divide && return cell_divide!(tree, e[2:end]...)
    e[1] == :remove && return cell_die!(tree, e[2:end]...)
    e[1] == :step && return nothing
    e[1] == :update && return foreach(c -> cell_update!(tree, c, e[3]), e[2])
    mark!(tree)
end

# Replays `events`, shuffling events within a step if `rng` is given, and pruning lazily after each step if `prune=true`
function replay(events; rng=nothing, prune=false, kwargs...)
    tree = Tree{Cell,Float64}(; kwargs...)
    for e in events
        e isa Vector ? foreach(x -> apply!(tree, x), isnothing(rng) ? e : shuffle(rng, e)) : apply!(tree, e)
        prune && e isa Tuple && e[1] == :step && prune!(tree; lazy=true)
    end

    sort!(tree)
end

overlaps(node, tspan) = node.t_birth <= tspan[2] && (node.alive || node.t_end > tspan[1])

@testset "Random trees (seed $seed)" for seed in 1:parse(Int, get(ENV, "TREES_SEEDS", "10"))
    rng = Xoshiro(seed)
    events, ref = random_events(rng)
    tree = replay(events)
    pruned = replay(events; prune=true, min_sweep=10)
    (; parent, t_birth, t_end, alive) = tree.nodes
    ids = eachindex(tree)

    # Canonical order does not depend on the order of events within a step
    @test tree.nodes == replay(events; rng).nodes
    @test pruned.nodes == replay(events; rng, prune=true, min_sweep=10).nodes

    # Recording: ancestries of living cells match the events
    for c in ref.cells
        expected = Float64[]
        p = ref.parent[c]
        while !isnothing(p)
            push!(expected, ref.birth[p])
            p = ref.parent[p]
        end
        @test [ tree[i].t_birth for i in ancestors(tree, tree.living[c]) ] == expected
    end

    # Queries agree with brute force
    @test all(parent[i] < i for i in ids) && issorted(t_birth)
    @test roots(tree) == findall(iszero, parent)
    @test living(tree) == findall(alive)
    @test all(collect(children(tree, i)) == findall(==(i), parent) for i in ids)
    @test all(has_divided(tree, i) == (i in parent) for i in ids)

    for s in [ 0., 3., 7.5, 12., 29.9, rand(rng) * 30 ]
        @test alive_at(tree, s) == [ i for i in ids if t_birth[i] <= s && (alive[i] || s < t_end[i]) ]
    end

    tspan = (5 * rand(rng), 5 + 20 * rand(rng))
    for i in [ roots(tree); rand(rng, ids, 5) ]
        chain = Int[]
        p = parent[i]
        while p != 0
            push!(chain, p)
            p = parent[p]
        end
        @test collect(ancestors(tree, i)) == chain
        @test collect(ancestors(tree, i; tspan)) == filter(j -> overlaps(tree[j], tspan), chain)

        desc = [ j for j in ids if i in ancestors(tree, j) ]
        pre = collect(descendants(tree, i))
        @test sort(pre) == desc
        @test all(parent[j] == i || findfirst(==(parent[j]), pre) < k for (k, j) in enumerate(pre))
        @test sort(collect(descendants(tree, i; tspan))) == filter(j -> overlaps(tree[j], tspan), desc)
    end

    # Pruning keeps the ancestors of living cells and the marked population, and nothing else
    @test sort([ pruned[i].t_birth for i in alive_at(pruned, ref.t_mark) ]) ==
          sort([ tree[i].t_birth for i in alive_at(tree, ref.t_mark) ])
    for c in ref.cells
        @test [ pruned[i].t_birth for i in ancestors(pruned, pruned.living[c]) ] ==
              [ tree[i].t_birth for i in ancestors(tree, tree.living[c]) ]
    end

    # Sweeps are lazy; after an explicit sweep, all leaves are living or marked, and sweeping again changes nothing
    prune!(pruned)
    @test all(has_divided(pruned, i) || pruned[i].alive || pruned[i].marked for i in eachindex(pruned))
    nodes = copy(pruned.nodes)
    prune!(pruned)
    @test pruned.nodes == nodes

    # Statistics agree across node collections and with brute force
    tb(node) = node.t_birth
    wf(node) = 1 + sin(node.t_birth)^2
    @test mean(tb, tree) ≈ mean(tb, filter(Returns(true), tree)) ≈ mean(tb, view(tree, ids)) ≈ mean(t_birth)
    @test var(tb, tree, wf) ≈ var(t_birth, pweights(wf.(tree.nodes)); corrected=true)

    s = 10 * rand(rng) + 5
    snap = snapshot((node, a) -> node.t_birth + a, tree, s)
    @test snap.id == filter(i -> !isnothing(tree[i].data), alive_at(tree, s)) && all(snap.value .≈ s)

    rec = at_age((node, a) -> a, tree, 1.5)
    @test rec.id == [ i for i in ids if !isnothing(tree[i].data) && lifetime(tree[i]) >= 1.5 ]
end
