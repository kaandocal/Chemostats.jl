# Cells

```@meta
CurrentModule = Chemostats
```

Cells in Chemostats.jl can be of any user-defined type that implements the [cell interface](@ref cell-interface). For models defined using [DifferentialEquations.jl](@extref DifferentialEquations :doc:`index`), the [`DECell`](@ref) class provides a direct wrapper.

## [DifferentialEquations.jl models](@id decell-page)

[`DECell`](@ref) also supports models defined using the wider SciML ecosystem, e.g. using [ModelingToolkit.jl](@extref ModelingToolkit :doc:`index`) and [Catalyst.jl](@extref Catalyst :doc:`index`).

To use this functionality we require:
* An `AbstractDEProblem` describing how the cells behave
* A `DivideCallback` provided to the `AbstractDEProblem`
* A `divide` function

Chemostats.jl supports `AbstractDEProblem`s with arbitrary parameters `p` and state vectors `u`. This problem should include a `DivideCallback` to determine when a cell has reached the end of its life, either due to division or death. When a cell reaches the end of its life, the `divide` function is called to determine the daughter cells. 

The `divide` argument should be a function of the form `int -> offspring` that takes a single argument `int`, the [integrator](@extref DifferentialEquations integrator) of the cell, and returns an [`Offspring`](@ref) object describing the daughter cells. This function can return `nothing` if there are no offspring. Alternatively, a fixed `Offspring` object can be passed directly.

```julia
# When the cell divides, the volume is halved, ribosomes are partitioned binomially,
# and each daughter draws a new cell cycle duration
Offspring(2, :V => ExactPartition(), :R => BinomialPartition(), :τ => Gamma(5, 1))

# Alternative: die with probability 1/10, otherwise divide asymmetrically
function divide(int)
    rand() < 0.1 && return nothing

    f = rand(Beta(20, 20))
    Offspring(:V => ExactPartition([f, 1 - f]), :R => BinomialPartition([f, 1 - f]))
end
```

**Deprecated:** `divide` can return a vector of `NamedTuple`s with fields `u0` and `p`, the initial state and parameters of each daughter cell. This is needed to change parameters without symbolic names, e.g. a `NamedTuple` of parameters.

```@docs
Chemostats.DECell
Chemostats.DivideCallback
Chemostats.Offspring
Chemostats.ExactPartition
Chemostats.BinomialPartition
Chemostats.AllOrNothingPartition
```

## [Custom cells](@id cell-interface)

Cells of a user-defined type must implement the following functions.

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
```
