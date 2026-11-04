Base.getindex(tree::Tree, id) = tree.nodes[id]

"""
    parent(tree::Tree, id)

Returns the parent of node `id`, or `0` if it is a root.
"""
Base.parent(tree::Tree, id) = tree.nodes.parent[id]

"""
    children(tree::Tree, id)

Returns the children of node `id` as a node collection. Requires a sorted tree.
"""
children(tree::Tree, id) = view(tree, child_ids(tree, id))

# Ids of the children of `id`, which are contiguous in a sorted tree
function child_ids(tree::Tree, id)
    ensure_sorted(tree)
    (; parent, t_birth) = tree.nodes
    n = length(tree)
    t = tree[id].t_end
    tree.nodes.alive[id] && return 1:0

    # Galloping search for the first node born at or after `t`
    lo, step = id, 1
    while lo + step <= n && t_birth[lo + step] < t
        lo += step
        step *= 2
    end

    hi = min(lo + step, n)
    c = searchsortedfirst(view(t_birth, lo+1:hi), t) + lo

    # Nodes born at `t` include the children of `id`, which are contiguous
    while c <= n && t_birth[c] == t && parent[c] != id
        c += 1
    end

    c_end = c - 1
    while c_end < n && t_birth[c_end + 1] == t && parent[c_end + 1] == id
        c_end += 1
    end

    c:c_end
end

"""
    has_divided(tree::Tree, id)

Returns `true` if node `id` has children in the population tree. Note that pruning may remove unmarked children. Requires a sorted tree.
"""
has_divided(tree::Tree, id) = !isempty(child_ids(tree, id))

"""
    lifetime(node::Node)

Returns the lifetime of a node, which is `t_end - t_birth`. For living cells, returns their age at their last update.
"""
lifetime(node::Node) = node.t_end - node.t_birth

"""
    is_alive(node::Node)

Returns `true` if `node` is a living cell.
"""
is_alive(node::Node) = node.alive

"""
    roots(tree::Tree)

Returns all root nodes.
"""
roots(tree::Tree) = view(tree, findall(iszero, tree.nodes.parent))

"""
    living(tree::Tree)

Returns all living cells, as a node collection.
"""
living(tree::Tree) = view(tree, sort!(collect(values(tree.living))))

"""
    alive_at(tree::Tree, t)
    alive_at(nodes, t)

Returns the nodes in `nodes` alive at time `t`. Requires a sorted tree.
"""
function alive_at(tree::Tree, t)
    ensure_sorted(tree)
    (; t_birth, t_end, alive) = tree.nodes

    # Finished nodes alive at `t` were born less than `max_lifetime` before `t`
    lo = searchsortedfirst(t_birth, t - tree.max_lifetime)
    hi = searchsortedlast(t_birth, t)
    ids = [ i for i in lo:hi if !alive[i] && t_end[i] > t ]

    for i in values(tree.living)
        t_birth[i] <= t && push!(ids, i)
    end

    view(tree, sort!(ids))
end
