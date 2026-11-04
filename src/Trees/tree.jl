"""
    struct Node{P}
        parent::Int
        generation::Int
        t_birth::Float64
        t_end::Float64
        alive::Bool
        marked::Bool
        data::Union{Nothing,P}
    end

A node in a [`Tree`](@ref), representing a single cell. Roots have `parent = 0`. For living cells (`alive = true`), `t_end` is the time up to which the cell was last updated (see [`cell_update!`](@ref)).
"""
struct Node{P}
    parent::Int
    generation::Int
    t_birth::Float64
    t_end::Float64              # for living cells: time of the last update
    alive::Bool
    marked::Bool                # kept by pruning
    data::Union{Nothing,P}
end

function Base.show(io::IO, node::Node{P}) where {P}
    print(io, "Node{", P, "}(", node.t_birth, " → ", node.t_end, ", data = ")
    show(IOContext(io, :compact => true, :limit => true), node.data)
    print(io, ")")
end

mutable struct Shears
    min_sweep::Int
    last_sweep_size::Int        # non-frozen nodes kept by the last sweep
    n_frozen::Int               # nodes 1:n_frozen are frozen
end

"""
    Tree{C,P=Any}(; min_sweep=1024)

Creates a population tree. Here `C` is the type of a cell, and `P` is the type of the data saved for each cell (at division, death or updates), `Any` if omitted. Each cell is stored in a [`Node`](@ref).

`min_sweep` is a parameter used for pruning (see [`prune!`](@ref)).

Every node in the tree has a unique ID, and nodes can be accessed as `tree[id]`. IDs are **not** stable, and can change when the tree is pruned or sorted.
"""
mutable struct Tree{C,P,V<:StructVector{Node{P}}}
    nodes::V
    living::IdDict{C,Int}       # living cell => node
    max_lifetime::Float64       # over finished nodes, used by `alive_at`
    sorted::Bool                # nodes are in canonical order
    shears::Shears
end

function Tree{C,P}(; min_sweep::Int=1024) where {C,P}
    nodes = StructVector{Node{P}}(undef, 0)
    Tree{C,P,typeof(nodes)}(nodes, IdDict{C,Int}(), 0., true, Shears(min_sweep, 0, 0))
end

Tree{C}(; kwargs...) where {C} = Tree{C,Any}(; kwargs...)

Base.length(tree::Tree) = length(tree.nodes)

n_frozen(tree::Tree) = tree.shears.n_frozen

function Base.show(io::IO, tree::Tree{C}) where {C}
    print(io, "Tree{$C}($(length(tree)) nodes, $(length(tree.living)) living")
    n_frozen(tree) > 0 && print(io, ", $(n_frozen(tree)) frozen")
    print(io, ")")
end

###

# Calls must be serialised. Events of different cells can arrive in any order; the events of a single cell arrive in order.

"""
    cell_init!(tree, cells, t)

Add `cells` to the population tree at time `t` (usually called at `t = 0`).
"""
function cell_init!(tree::Tree, cells, t)
    for cell in cells
        @argcheck !haskey(tree.living, cell) "Cell is already in the tree"
        tree.living[cell] = push_node!(tree, 0, 0, t)
    end

    nothing
end

"""
    cell_divide!(tree, parent, children, t, data)

Called when `parent` divides into `children` at time `t`. Here `data` is the data saved for `parent`.
"""
function cell_divide!(tree::Tree, parent, children, t, data)
    id = living_id(tree, parent)
    @argcheck allunique(objectid, children) "Children must be distinct"
    @argcheck !any(c -> haskey(tree.living, c), children) "Children must not be living cells of the tree"
    finish_node!(tree, id, t, data)
    delete!(tree.living, parent)

    for cell in children
        tree.living[cell] = push_node!(tree, id, tree.nodes.generation[id] + 1, t)
    end

    nothing
end

"""
    cell_die!(tree, cell, t, data)

Called when `cell` dies without dividing. Here `data` is the data saved for `cell`. A cell that dies at its time of birth is removed from the tree.
"""
function cell_die!(tree::Tree, cell, t, data)
    id = living_id(tree, cell)
    # A cell removed at birth was never alive, so it is dropped
    t == tree.nodes.t_birth[id] ? drop_node!(tree, id) : finish_node!(tree, id, t, data)
    delete!(tree.living, cell)
    nothing
end

"""
    cell_update!(tree, cell, t, data)
    cell_update!(tree, cell, t)

Update the cell's current time to `t`; if `data` is provided, updates the cell's saved data. Should be called for each living cell before doing tree operations (e.g. at the end of a simulation).
"""
function cell_update!(tree::Tree, cell, t, data)
    id = living_id(tree, cell)
    set_end!(tree, id, t)
    tree.nodes.data[id] = data
    nothing
end

function cell_update!(tree::Tree, cell, t)
    set_end!(tree, living_id(tree, cell), t)
    nothing
end

function push_node!(tree::Tree{C,P}, parent, generation, t) where {C,P}
    push!(tree.nodes, Node{P}(parent, generation, t, t, true, false, nothing))
    tree.sorted = false
    length(tree)
end

function living_id(tree::Tree, cell)
    id = get(tree.living, cell, 0)
    @argcheck id != 0 "Cell is not a living cell of the tree"
    id
end

function set_end!(tree::Tree, id, t)
    @argcheck t >= tree.nodes.t_birth[id] "Time is before the birth of the cell"
    tree.nodes.t_end[id] = t
end

# Removes the living leaf `id`, shifting later ids down. Keeps the order, so a sorted tree stays sorted.
function drop_node!(tree::Tree, id)
    (; nodes) = tree
    deleteat!(nodes, id)
    for i in id:length(nodes)
        nodes.parent[i] > id && (nodes.parent[i] -= 1)
    end

    for (cell, i) in collect(tree.living)
        i > id && (tree.living[cell] = i - 1)
    end
end

function finish_node!(tree::Tree, id, t, data)
    (; nodes) = tree
    set_end!(tree, id, t)
    nodes.alive[id] = false
    nodes.data[id] = data
    tree.max_lifetime = max(tree.max_lifetime, t - nodes.t_birth[id])
    nothing
end

###

"""
    sort!(tree::Tree)

Sorts nodes in a tree. Some functions (e.g. [`children`](@ref) and [`descendants`](@ref)) require the tree to be sorted. Does nothing if the tree is already sorted.
"""
Base.sort!(tree::Tree) = tree.sorted ? tree : compact!(tree, nothing)
Base.issorted(tree::Tree) = tree.sorted

# Keeps the non-frozen nodes `f + k` with `keep[k]` (all if `keep === nothing`) in canonical order, in place.
# Only allocates O(number of non-frozen nodes). Sorting by (birth time, generation) puts parents before
# children; ties are broken by new parent id, then old id (keeping siblings in order).
function compact!(tree::Tree, keep)
    (; parent, generation, t_birth) = tree.nodes
    n, f = length(tree), n_frozen(tree)

    birth_key(i) = (t_birth[i], generation[i])
    rows = isnothing(keep) ? collect(f+1:n) : [ i for i in f+1:n if keep[i - f] ]
    sort!(rows; by = birth_key)

    new_id = zeros(Int, n - f)
    new_parent(i) = (p = parent[i]; p <= f ? p : new_id[p - f])

    next, lo = f, 1
    while lo <= length(rows)
        hi = lo
        while hi < length(rows) && birth_key(rows[hi + 1]) == birth_key(rows[lo])
            hi += 1
        end

        # Parents have been renumbered in earlier groups
        group = sort!(view(rows, lo:hi); by = i -> (new_parent(i), i))
        for i in group
            next += 1
            new_id[i - f] = next
        end

        lo = hi + 1
    end

    # Copy the kept nodes, then write them back after the frozen nodes
    active = tree.nodes[rows]
    active.parent .= new_parent.(rows)
    resize!(tree.nodes, f + length(rows))
    foreach((dst, src) -> copyto!(dst, f + 1, src, 1, length(src)),
            StructArrays.components(tree.nodes), StructArrays.components(active))

    # Living cells are never frozen. Collect first, as assigning while iterating can rehash the `IdDict`.
    for (cell, id) in collect(tree.living)
        tree.living[cell] = new_id[id - f]
    end

    tree.sorted = true
    tree
end
