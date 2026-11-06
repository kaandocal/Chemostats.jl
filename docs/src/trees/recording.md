# Recording trees

```@meta
CurrentModule = Trees
```

## Recording

The population tree is modified using the following functions:

| Function | Called when |
|---|---|
| [`cell_init!(tree, cells, t)`](@ref cell_init!) | add `cells` to the population (with no parent) |
| [`cell_divide!(tree, parent, children, t, data)`](@ref cell_divide!) | `parent` divides into `children` |
| [`cell_die!(tree, cell, t, data)`](@ref cell_die!) | `cell` dies |
| [`cell_update!(tree, cell, t[, data])`](@ref cell_update!) | cell is updated to time `t` |

[`Tree`](@ref) stores one [`Node`](@ref) per cell in a `StructVector` `tree.nodes`. Each node stores information about the parent, generation, and lifetime of each cell, and a user-defined field `data`.

Nodes are identified by their id and can be accessed as `tree[id]`. Cells are stored approximately in order, but some operations require the tree to be exactly sorted, achieved using [`sort!(tree)`](@ref Base.sort!(::Tree)). 

**Note:** The id of a node can change when calling [`sort!`][@ref] or [`prune!`](@ref). 

For large trees, [`prune!(tree)`](@ref prune!) performs garbage collection by removing all nodes that do not have a living or marked descendant (dead branches). You can call [`mark!`](@ref) to save cells that should not be garbage collected. 
