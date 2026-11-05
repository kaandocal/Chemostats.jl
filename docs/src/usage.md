# Usage 

Chemostats.jl simulates individual cells in a [`Chemostat`](@ref) object.

## Chemostats

```@docs
Chemostat
Chemostats.simulate
Chemostats.simulate!
Chemostats.CellException
```

## Estimating population sizes

```@docs
Chemostats.Snapshot
Chemostats.est_Λ
Chemostats.est_N
Chemostats.est_logN
```

## The Cell API 

Cells in Chemostats.jl can be of any user-defined type that implements the following functions. For models defined using [DifferentialEquations.jl](@extref DifferentialEquations :doc:`index`), the [`DECell`](@ref decell-page) class provides a direct wrapper.

```@docs
Chemostats.get_curr_t
Chemostats.get_state
Chemostats.divide!
Chemostats.die!
Chemostats.kill!
Chemostats.get_children
Chemostats.clone_cell
Chemostats.get_cell_rng
Chemostats.get_data
Chemostats.CellState
```
