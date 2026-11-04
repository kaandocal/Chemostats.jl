struct AncestorIterator{T<:Tree}
    tree::T
    start::Int
    tspan::Tuple{Float64,Float64}
end

"""
    ancestors(tree::Tree, id; self=false, tspan = (-Inf, Inf))
    ancestors(tree::Tree, cell; self=false, tspan = (-Inf, Inf))

Iterates over the ancestors of node `id` (its backward lineage), starting with its parent. If `self=true`, start with the cell itself. Only cells alive at some time in `tspan` are included. `cell` must be a living cell of the tree.
"""
function ancestors(tree::Tree, id::Integer; self::Bool=false, tspan = (-Inf, Inf))
    AncestorIterator(tree, self ? Int(id) : tree[id].parent, Float64.(tspan))
end

# Ancestors of a living cell
ancestors(tree::Tree, cell; kwargs...) = ancestors(tree, living_id(tree, cell); kwargs...)

function Base.iterate(iter::AncestorIterator, id=iter.start)
    (; parent, t_birth, t_end, alive) = iter.tree.nodes
    tspan = iter.tspan

    # Skip ancestors born after `tspan`
    while id != 0 && t_birth[id] > tspan[2]
        id = parent[id]
    end

    # Ancestors are visited backwards in time, so stop at the first one that ended before `tspan`
    (id == 0 || (!alive[id] && t_end[id] <= tspan[1])) && return nothing
    id, parent[id]
end

Base.IteratorSize(::Type{<:AncestorIterator}) = Base.SizeUnknown()
Base.eltype(::Type{<:AncestorIterator}) = Int

struct DescendantIterator{T<:Tree}
    tree::T
    root::Int
    self::Bool
    tspan::Tuple{Float64,Float64}
end

"""
    descendants(tree::Tree, id; self=false, tspan = (-Inf, Inf))

Iterates over the descendants of node `id` (its forward lineages) in depth-first manner. If `self=true`, start with the cell itself. Only cells alive at some time in `tspan` are included. Requires a sorted tree.
"""
function descendants(tree::Tree, id; self::Bool=false, tspan = (-Inf, Inf))
    ensure_sorted(tree)
    DescendantIterator(tree, Int(id), self, Float64.(tspan))
end

# Next node after `i` in preorder within the subtree of `iter.root`, or 0. Subtrees born after `tspan` are skipped.
function next_descendant(iter::DescendantIterator, i)
    (; tree, root, tspan) = iter
    (; parent, t_birth) = tree.nodes

    cs = children(tree, i)
    !isempty(cs) && t_birth[first(cs)] <= tspan[2] && return first(cs)

    while i != root
        p = parent[i]
        i < length(tree) && parent[i + 1] == p && return i + 1      # next sibling
        i = p
    end

    0
end

function Base.iterate(iter::DescendantIterator, i = iter.self ? iter.root : next_descendant(iter, iter.root))
    (; t_birth, t_end, alive) = iter.tree.nodes
    tspan = iter.tspan

    while i != 0
        t_birth[i] <= tspan[2] && (alive[i] || t_end[i] > tspan[1]) && return i, next_descendant(iter, i)
        i = next_descendant(iter, i)
    end

    nothing
end

Base.IteratorSize(::Type{<:DescendantIterator}) = Base.SizeUnknown()
Base.eltype(::Type{<:DescendantIterator}) = Int
