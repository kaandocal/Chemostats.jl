# Analysis

```@meta
CurrentModule = Trees
```

We use a small random tree in which the data of each cell is a random number.

```@example analysis
using Trees
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

## Node collections

Analysis functions take *node collections*:

- a tree (all of its nodes, iterated with `values(tree)`),
- `view(tree, ids)` or [`filter(f, nodes)`](@ref Base.filter(::Any, ::NodeCollection)) (e.g. cells born after a burn-in),
- [`ancestors(tree, id)`](@ref ancestors), the backward lineage of a node,
- [`descendants(tree, id)`](@ref descendants), its forward lineages.

Lineages can be restricted to cells alive at some time within a time span `tspan`. None of these copy any nodes.
Node collections iterate over nodes; `keys` returns their ids, and `pairs` returns `id => node`.

```@example analysis
id = first(keys(living(tree)))
lineage = ancestors(tree, id; self=true, tspan=(5., 20.))
late = filter(node -> node.t_birth >= 5, tree)
length(collect(lineage)), length(late)
```

## Observables and summaries

Trees uses two kinds of user-defined functions:

- **Observables** `f(node, a)` of a cell at age `a` (time since birth), usually evaluating its data `node.data`. Using
  the age makes observables independent of the time frame of the data, e.g. `sol(sol.t[1] + a)` for a solution object
  `sol = node.data`. The absolute time is `node.t_birth + a`.
- **Per-cell summaries** `f(node)` of a [`Node`](@ref), e.g. its lifetime, or the time average of an observable over
  a cell's life, [`time_average(f, node; tspan)`](@ref time_average). Time averages of living cells extend up to
  their last update.

Nodes without data, and nodes for which a function returns `nothing`, are skipped.

## Population and lineage statistics

[`mean`](@ref StatsBase.mean(::Any, ::NodeCollection, ::Any)), [`var`](@ref StatsBase.var(::Any, ::NodeCollection, ::Any))
and [`mean_and_var`](@ref StatsBase.mean_and_var(::Any, ::NodeCollection, ::Any)) average per-cell summaries over node
collections, optionally with weights. Weights are functions of nodes, and are evaluated in a single pass without
allocating per node. [`time_weights`](@ref) weights cells by the time they were alive
(within a time span), so that weighted means of per-cell time averages are time averages over the population.

```@example analysis
age(node, a) = a
tspan = (5., 20.)
avg_age = time_average(age; tspan)

mean(avg_age, late),                                    # per cell
mean(avg_age, tree, time_weights(; tspan)),         # time-weighted, over the population
mean(avg_age, lineage, time_weights(; tspan))    # time-weighted, along a lineage
```

Variances of per-cell summaries are not population variances: for the latter, use per-cell averages of `x` and `x²`.

## Snapshots and records

[`snapshot(f, nodes, t)`](@ref snapshot) evaluates an observable over all cells alive at time `t`, or over several times.
On a lineage, a snapshot is the single lineage cell alive at `t`. Population averages are unweighted, since random
removal keeps a uniform sample of the population.

```@example analysis
[ mean(snapshot(age, tree, t)) for t in (10., 15.) ]
```

Snapshots and age-conditioned records return vectors of values. Observables can return several values, e.g. to pair
the birth and division times of each cell:

```@example analysis
div = at_division((node, τ) -> (node.t_birth, node.t_birth + τ), late)     # birth and division times
mean(last.(div) .- first.(div)), [ mean(at_age(age, late, a)) for a in (0., 0.5) ]
```

[`at_age`](@ref) includes removed cells if they reached age `a`, and living cells up to their age at their last update.
