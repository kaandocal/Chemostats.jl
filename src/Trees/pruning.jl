"""
    mark!(tree)
    mark!(tree, ids)

Marks the nodes `ids` to be kept by pruning, together with their ancestors. `mark!(tree)` marks the current population, i.e. all living cells.
"""
function mark!(tree::Tree, ids = living(tree))
    for id in ids
        tree.nodes.marked[id] = true
    end

    tree
end

"""
    prune!(tree::Tree; lazy=false)

Removes all "dead" nodes, i.e. nodes that are not ancestors of living cells or of marked nodes (see [`mark!`](@ref)). If `lazy=true`, only prunes if the tree has doubled in size since the last pruning.
"""
function prune!(tree::Tree; lazy::Bool=false)
    (; parent, marked) = tree.nodes
    (; min_sweep, last_sweep_size) = tree.shears
    n, f = length(tree), n_frozen(tree)
    lazy && n - f <= max(2 * last_sweep_size, min_sweep) && return tree

    # Only non-frozen nodes are marked
    keep = falses(n - f)
    for id in values(tree.living)
        mark_ancestors!(keep, parent, id, f)
    end

    for id in f+1:n
        marked[id] && mark_ancestors!(keep, parent, id, f)
    end

    compact!(tree, keep)
    tree.shears.last_sweep_size = length(tree) - f
    freeze!(tree)
end

# Marks `id` and its ancestors, up to the first marked or frozen node
function mark_ancestors!(keep, parent, id, f)
    while id > f && !keep[id - f]
        keep[id - f] = true
        id = parent[id]
    end
end

# Nodes born before the most recent common ancestor (MRCA) of the living cells are ancestors of the MRCA
# or of marked nodes, so later sweeps never remove them. Requires canonical order.
function freeze!(tree::Tree)
    isempty(tree.living) && return tree
    (; parent, t_birth) = tree.nodes
    n, f = length(tree), n_frozen(tree)

    # Number of living descendants of non-frozen nodes; one backward pass suffices as children come after parents
    n_desc = zeros(Int, n - f)
    for id in values(tree.living)
        n_desc[id - f] += 1
    end

    mrca = 0
    for i in n:-1:f+1
        if n_desc[i - f] == length(tree.living)
            mrca = i
            break
        end

        p = parent[i]
        p > f && (n_desc[p - f] += n_desc[i - f])
    end

    mrca == 0 && return tree

    tree.shears.n_frozen = f + searchsortedfirst(view(t_birth, f+1:n), t_birth[mrca]) - 1
    tree.shears.last_sweep_size = n - tree.shears.n_frozen
    tree
end
