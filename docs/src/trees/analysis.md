# Analysis

```@meta
CurrentModule = Chemostats.Trees
```

We use a small random tree in which the data of each cell is a random number.

```@example analysis
using Chemostats.Trees
using Chemostats.Trees: cell_init!, cell_divide!, cell_die!, cell_update!
using Random
using Statistics

mutable struct Cell end

Random.seed!(1)
tree = Tree{Cell,Float64}()
cells = [ Cell() for _ in 1:20 ]
cell_init!(tree, cells, 0.)

for t in 0.1:0.1:20
    i = rand(eachindex(cells))
    c = cells[i]
    kids = [ Cell(), Cell() ]
    cell_divide!(tree, c, kids, t, rand())
    cells[i] = kids[1]

    # Keep the population size constant by replacing a random cell
    j = rand(eachindex(cells))
    cell_die!(tree, cells[j], t, rand())
    cells[j] = kids[2]
end

foreach(c -> cell_update!(tree, c, 20., rand()), cells)
sort!(tree)
```

## Selecting cells

Most analysis functions take arbitrary node collections, which can be:
- a tree, using `values(tree)`
- a subtree `view(tree, ids)`
- cells satisfying certain conditions, using [`filter(f, nodes)`](@ref Base.filter(::Any, ::NodeCollection))
- ancestors of a cell using [`ancestors(tree, id)`](@ref ancestors)
- descendants of a cell using [`descendants(tree, id)`](@ref descendants)
- root cells, living cells or the children of a cell, using [`roots(tree)`](@ref roots), [`living(tree)`](@ref living) and [`children(tree, id)`](@ref children)

For these collections, `keys` returns the cell ids, `values` the nodes, and `pairs` returns `id => node`.

```@example analysis
id = first(keys(living(tree)))
lineage = ancestors(tree, id; self=true, tspan=(5., 20.))
late = filter(node -> node.t_birth >= 5, tree)
length(collect(lineage)), length(late)
```

```@docs
Base.view(::Tree, ::Any)
Base.filter(::Any, ::NodeCollection)
ancestors
descendants
roots
living
children
alive_at
alive_during
```

## Observables and summaries

Trees uses two kinds of user-defined functions:

- Observables `f(node, a)` of a cell at age `a` (time since birth). This could be a cell's position, size, protein content etc.
- Per-cell summaries `f(node)` of a [`Node`](@ref), e.g. its lifetime, or the time average of an observable using [`time_average(f, node)`](@ref time_average). 

```@docs
time_average
```

## Population and lineage statistics

[`mean`](@ref StatsBase.mean(::Any, ::NodeCollection, ::Any)), [`var`](@ref StatsBase.var(::Any, ::NodeCollection, ::Any)) and [`mean_and_var`](@ref StatsBase.mean_and_var(::Any, ::NodeCollection, ::Any)) are defined for trees and lineages. To compute time-averages, use [`time_weights`](@ref) to weight cells by their lifetime.

```@example analysis
age(node, a) = a
avg_age = time_average(age)

mean(avg_age, late),                                    # per cell
mean(avg_age, tree, time_weights),         # time-weighted, over the population
mean(avg_age, lineage, time_weights)    # time-weighted, along a lineage
```

```@docs
time_weights
StatsBase.mean(::Any, ::NodeCollection, ::Any)
StatsBase.var(::Any, ::NodeCollection, ::Any)
StatsBase.mean_and_var(::Any, ::NodeCollection, ::Any)
```

## Snapshots

[`snapshot(f, nodes, t)`](@ref snapshot) evaluates an observable over all cells alive at time `t`. `t` can be an array, in which case multiple snapshots are pooled.

```@example analysis
[ mean(snapshot(age, tree, t)) for t in (10., 15.) ]
```

Snapshots and age-conditioned records return vectors of values. Observables can return several values, e.g. to pair the birth and division times of each cell:

```@example analysis
div = at_division((node, τ) -> (node.t_birth, node.t_birth + τ), late)     # birth and division times
```

```@docs
snapshot
at_age
at_birth
at_division
```
