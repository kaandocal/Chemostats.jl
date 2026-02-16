# ![](docs/src/assets/logo.svg) Chemostats.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://kaandocal.github.io/Chemostats.jl/dev)

Chemostats.jl is a package to efficiently simulate cell populations in Julia, featuring state-of-the-art algorithms that efficiently estimate population dynamics at sub-exponential cost. It is compatible with the [SciML](https://sciml.ai) ecosystem, including [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) and [Catalyst.jl](https://github.com/SciML/Catalyst.jl)

This package is experimental - any feedback is appreciated (either by email or by opening an issue on GitHub).

## Example

(This does not quite work yet...)

```julia
using Catalyst
using Chemostats
using JumpProcesses
using CairoMakie

# Define a reaction network in which ribosomes are produced in a volume-dependent manner
# Time is measured in units
rn = @reaction_network begin
    @species R(t) = 0.
    @variables V(t) = V_d / 2
    @parameters σ = 10. λ = 1. V_d = 2.

    σ * V(t), 0 --> R

    @equations begin
        D(V) ~ λ * R
    end
end 

# A cell divides once it doubles its initial volume. Implemented with a callback.
cb_div = Chemostats.DivideCallback(rn, rn.V ~ rn.V_d)

# This function determines what offspring a dividing cell produces.
function divide(int)
    rand() < 0.1 && return nothing      # Die with probability 1/10

    # Create two cells with half the volume and (roughly) half the ribosomes
    # The model parameters are inherited from the parent cell
    cell_1 = (u0 = (V = int[:V] / 2, R = rand(Binomial(Int(int[:R]), 0.5))), p=nothing)
    cell_2 = (u0 = (V = int[:V] / 2, R = int[:R] - cell_1.u0.R), p=nothing)
    
    cell_1, cell_2
end

prob = ODEProblem(rn, [], (0., 1.); callbacks=cb_div) # tspan here does not matter
cells = DECell[ DECell(prob, divide) ]
chem = Chemostat(cells)

# Simulate population for a day
Chemostats.simulate!(chem, 24 * 60, Chemostats.Strict(10); saveat=0:10:100.)

# The full population would require simulating a lot of cells!
est_Λ(chem)

plot(chem)
```

### See also
* [Agents.jl](https://juliadynamics.github.io/Agents.jl): A Julia framework for agent-based modelling. More general than Chemostats.jl, which is targeted towards efficient estimation of cell population dynamics.

### References 
* Giardinà, Kurchan, Peliti. "Direct evaluation of large-deviation functions," Phys. Rev. Lett. 96(12), 2006
* Lecomte & Tailleur. "A numerical approach to large deviations in continuous time," J. Stat. Mech.: Theory Exp. 2007(3), 2007

### Acknowledgments
* Thanks to [Chris Rackaukas](https://www.chrisrackauckas.com), [Aayush Sabharwal](https://github.com/AayushSabharwal), [Sam Isaacson](https://www.samisaacson.com), and the SciML community for their help and support (and their work on the SciML ecosystem!).
