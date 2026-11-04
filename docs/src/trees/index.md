# Trees.jl

```@meta
CurrentModule = Trees
```

Trees.jl is a utility package for manipulating cell lineages and populations. The goal is to provide a simple, flexible API for basic lineage analysis tasks, while being lightweight enough for large-scale simulations with millions of cells. 

## Features
- Can store cells with any data format
- Utilities for iterating over ancestors and descendants
- Statistical functions to compute lineage and population averages
- Automatic garbage collection (pruning) of trees

## Example

The following example constructs a small tree by hand, where each cell stores its birth volume. A cell divides once it has doubled in size, and its daughters inherit fractions `f` and `1 - f` of its volume.

```@example index
using Trees
using Trees: cell_init!, cell_divide!, cell_die!, cell_update!
using Statistics

mutable struct Cell
    v_birth::Float64
end

# Divides `cell` at time `t` into daughters with fractions `f` and `1 - f` of its volume
function divide!(tree, cell, t, f)
    v_div = 2 * cell.v_birth
    daughters = [ Cell(f * v_div), Cell((1 - f) * v_div) ]
    cell_divide!(tree, cell, daughters, t, cell.v_birth)
    daughters
end

r = Cell(1.)
tree = Tree{Cell,Float64}()
cell_init!(tree, [ r ], 0.)
a, b = divide!(tree, r, 1., 0.25)                # r divides at t = 1: a and b have volumes 0.5 and 1.5
a1, a2 = divide!(tree, a, 2., 0.5)               # a divides at t = 2: a1 and a2 have volume 0.5
b1, b2 = divide!(tree, b, 2.5, 0.4)              # b divides at t = 2.5: b1 and b2 have volumes 1.2 and 1.8
cell_die!(tree, a2, 2.5, a2.v_birth)             # a2 is removed at t = 2.5
for cell in (a1, b1, b2)                         # update living cells before analysing the tree
    cell_update!(tree, cell, 3., cell.v_birth)
end
sort!(tree)
```

Nodes can be queried by id, which relates cells to each other in the tree:

```@example index
id = first(keys(living(tree)))      # 4 (a1)
tree[id]                            # Node{Float64}(2.0 → 3.0, data = 0.5)
collect(keys(ancestors(tree, id)))  # [2, 1] (a, r)
```

We can extract population statistics of the form `f(node, a)`, where `a` is the cell age. Here, cells grow exponentially and double their volume in one time unit, so the volume of a cell at age `a` is:

```@example index
volume(node, a) = node.data * 2^a
snapshot(volume, tree, 1.5)                      # [0.707, 2.121], i.e. a and b
```

We can also compute population statistics such as time or per-cell averages, over entire populations, descendants, or ancestral lineages.

```@example index
mean(time_average(volume), tree, time_weights)                        # 1.489
mean(time_average(volume), descendants(tree, 3; self=true), time_weights)  # 2.300 (b, b1, b2)
mean(node -> node.data, ancestors(tree, id; self=true))              # 0.667 (a1, a, r)
```

## See also

- [tskit](https://tskit.dev), whose tables and `simplify` algorithm inspired the pruning in Trees.
