# [Using Chemostats.jl with DifferentialEquations.jl](@id decell-page)

The class [`DECell`](@ref) provides a convenient way to use models defined using [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). This is includes models defined using the wider SciML ecosystem, e.g. using [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) and [Catalyst.jl](https://github.com/SciML/Catalyst.jl).

## The [`DECell`](@ref) interface
To use this functionality we require:
* An `AbstractDEProblem` describing how the cells behave
* A `DividingCallback` provided to the `AbstractDEProblem` (to be implemented!)
* A `divide` function

Chemostats.jl supports `AbstractDEProblem`s with arbitrary parameters `p` and state vectors `u`. This problem should include a `DividingCallback` to determine when a cell has reached the end of its life, either due to division or death. 

When a cell reaches the end of its life, the `divide` function is called to determine the parameters and initial state vectors of the daughter cells. This function takes a single argument `int`, the [integrator](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator) of the cell. It should an array of daughter cells, or `nothing` if there are no offspring. Each daughter cell is determined by a `NamedTuple` with fields `p` containing the parameters and `u0` containing the initial state of the daughter cell defined by the `AbstractDEProblem`.

```@docs
Chemostats.DECell
Chemostats.DividingCallback
```

## To do:
* Symbolic callbacks defined in ModelingToolkit.jl