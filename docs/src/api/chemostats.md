# Chemostats

## Chemostat objects

Chemostats.jl simulates individual cells in a [`Chemostat`](@ref) object.

```@docs
Chemostat
Chemostats.simulate
Chemostats.simulate!
Chemostats.CellState
Chemostats.CellException
```

## Estimating population sizes

```@docs
Chemostats.Snapshot
Chemostats.est_growth_rate
Chemostats.est_N
Chemostats.est_logN
```
