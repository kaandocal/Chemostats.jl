# Tree API

```@meta
CurrentModule = Chemostats.Trees
```

## Trees and nodes

A [`Tree`](@ref) stores one [`Node`](@ref) per cell as a [`StructVector`](@extref StructArrays :doc:`index`) under `tree.nodes`. Each node records the parent, generation, birth time and end time of a cell, and a user-defined field `data`. For living cells, the end time is the time of their last update.

Nodes are internally identified by their id and can be accessed as `tree[id]`.  Since `tree.nodes` is a `StructVector`, fields can also be accessed for all cells at once, e.g. `tree.nodes.t_birth`.

Cells are stored approximately in order, but some operations require the tree to be exactly sorted, achieved using [`sort!(tree)`](@ref Base.sort!(::Tree)).

**Note:** The id of a node is not stable and can change when calling [`sort!`](@ref Base.sort!(::Tree)) or [`prune!`](@ref).

```@docs
Tree
Node
Base.sort!(::Tree)
Base.parent(::Tree, ::Any)
has_divided
lifetime
is_alive
```

## Adding cells

The population tree is modified using the following functions:

| Function | Called when |
|:---|:---|
| [`cell_init!(tree, cells, t)`](@ref cell_init!) | add `cells` to the population (with no parent) |
| [`cell_divide!(tree, parent, children, t, data)`](@ref cell_divide!) | `parent` divides into `children` |
| [`cell_die!(tree, cell, t, data)`](@ref cell_die!) | `cell` dies |
| [`cell_update!(tree, cell, t[, data])`](@ref cell_update!) | cell is updated to time `t` |

```@docs
cell_init!
cell_divide!
cell_die!
cell_update!
```

## Pruning cells

For large trees, [`prune!(tree)`](@ref prune!) performs garbage collection by removing unused branches: only cells that are currently alive or marked using [`mark!`](@ref) are retained, together with their ancestors.

```@docs
mark!
unmark!
prune!
```

## Export

Trees can be exported in Newick format, which can be read by many tree-handling packages.

```@docs
newick
```
