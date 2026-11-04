# Node collections: all nodes of a tree, a view of a tree, or the ancestors or descendants of a node.
# `keys` returns their ids, `values` their nodes, and `pairs` returns `id => node`. Trees are not iterable themselves.

struct TreeView{T<:Tree,I}
    tree::T
    ids::I
end

"""
    view(tree::Tree, ids)

Returns the nodes `ids` of `tree` as a node collection, without copying.
"""
Base.view(tree::Tree, ids) = TreeView(tree, ids)

const NodeCollection = Union{Tree,TreeView,AncestorIterator,DescendantIterator}

get_tree(tree::Tree) = tree
get_tree(v::TreeView) = v.tree
get_tree(lin::Lineage) = lin.tree

Base.keys(tree::Tree) = eachindex(tree.nodes)     # also gives `eachindex(tree)`
Base.keys(v::TreeView) = v.ids
Base.keys(lin::Lineage) = LineageIds(lin)
Base.values(tree::Tree) = view(tree, keys(tree))

# Iterates over the nodes of `c`, following the iteration of its ids
function Base.iterate(c::Union{TreeView,Lineage}, state...)
    next = iterate(keys(c), state...)
    isnothing(next) ? nothing : (get_tree(c)[first(next)], last(next))
end

node_type(::Type{<:Tree{C,P}}) where {C,P} = Node{P}
Base.eltype(::Type{<:TreeView{T}}) where {T} = node_type(T)
Base.length(v::TreeView) = length(v.ids)
Base.show(io::IO, v::TreeView) = print(io, "TreeView(", length(v), " nodes)")

# Functions needing time order require a sorted tree. They do not sort it themselves, which would renumber nodes.
function ensure_sorted(nodes)
    @argcheck issorted(get_tree(nodes)) "Tree must be sorted, call `sort!(tree)` first"
end

"""
    filter(f, nodes)

Returns the nodes of `nodes` for which `f(node)` is true, as a node collection (a view of the tree).
"""
function Base.filter(f, nodes::NodeCollection)
    view(get_tree(nodes), [ id for (id, node) in pairs(nodes) if f(node) ])
end

function alive_at(nodes::NodeCollection, t)
    filter(node -> node.t_birth <= t && (node.alive || t < node.t_end), nodes)
end

# Marks the nodes of a node collection, e.g. `mark!(tree, alive_at(tree, t))`
mark!(tree::Tree, nodes::Union{TreeView,Lineage}) = mark!(tree, keys(nodes))
