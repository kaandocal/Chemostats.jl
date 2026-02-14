# Usage 

Chemostats.jl simulates individual cells in a [`Chemostat`](@ref) object.

## Chemostats

```@docs
Chemostat
Chemostats.simulate!
Chemostats.get_snapshot
```

## Estimating population sizes

```@docs
Chemostats.Snapshot
Chemostats.est_Λ
Chemostats.est_N
Chemostats.est_logN
```

## The Cell API 

Cells in Chemostat.jl can be of an user-defined type that implements the following functions. For models defined using [`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl), the [`DECell`](@ref decell-page) class provides a direct wrapper.

(Mention `step!`)

```@docs
Chemostats.get_curr_t
Chemostats.get_state
Chemostats.divide!
Chemostats.die!
Chemostats.kill!
Chemostats.get_offspring
Chemostats.duplicate_cell
Chemostats.CellState
```
