# API

```@meta
CurrentModule = Trees
```

```@docs
Trees
```

## Recording

```@docs
cell_init!
cell_divide!
cell_die!
cell_update!
```

## Population trees

```@docs
Tree
Node
Base.sort!(::Tree)
mark!
prune!
```

## Queries

```@docs
roots
living
Base.parent(::Tree, ::Any)
children
alive_at
ancestors
descendants
Base.view(::Tree, ::Any)
Base.filter(::Any, ::NodeCollection)
has_divided
lifetime
is_alive
```

## Statistics

```@docs
time_average
time_weights
StatsBase.mean(::Any, ::NodeCollection, ::Any)
StatsBase.var(::Any, ::NodeCollection, ::Any)
StatsBase.mean_and_var(::Any, ::NodeCollection, ::Any)
snapshot
```

## Records

```@docs
at_age
at_birth
at_division
```

## Export

```@docs
newick
```
