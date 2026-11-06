# ![](docs/src/assets/logo.svg) Chemostats.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://kaandocal.github.io/Chemostats.jl/dev)
[![Run Tests](https://github.com/kaandocal/Chemostats.jl/actions/workflows/test.yml/badge.svg)](https://github.com/kaandocal/Chemostats.jl/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/kaandocal/Chemostats.jl/graph/badge.svg)](https://codecov.io/gh/kaandocal/Chemostats.jl)

Chemostats.jl is a package to efficiently simulate cell populations in Julia, featuring state-of-the-art algorithms that estimate population dynamics at sub-exponential cost. It is compatible with the [SciML](https://sciml.ai) ecosystem, including [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) and [Catalyst.jl](https://github.com/SciML/Catalyst.jl).

The theory behind this package is explained in our paper:

> Kaan Öcal, Lars Maas, Michael P.H. Stumpf. **"Probing bacterial fitness in virtual chemostats."** [*bioRxiv*](https://www.biorxiv.org/content/10.64898/2026.10.31.763286v1) (2026).

This package is somewhat experimental - any feedback is appreciated (either by email or by opening an issue).

## Example

A simple model of cell division, featuring volume growth coupled to ribosome production.

<!-- BEGIN demo/ribosomes.jl -->
```julia
using Catalyst, OrdinaryDiffEqTsit5, Chemostats

rn = @reaction_network begin
    @species R(t) = 0.
    @variables V(t) = 1.
    @parameters σ = 10. λ = 1. V_d = 2.

    σ * V, 0 --> R

    @equations begin
        D(V) ~ λ * R
    end
end;

cb_div = Chemostats.DivideCallback(rn.V ~ rn.V_d);

function divide(int)
    rand() < 0.1 && return nothing

    Chemostats.Offspring(2, :R => Chemostats.BinomialPartition(),
                            :V => Chemostats.ExactPartition())
end;

u0 = [ :R => 0., :V => 1. ]
p = [ :σ => 10., :λ => 1., :V_d => 2. ]
prob = HybridProblem(rn, u0, (0., Inf), p; callback = cb_div)

cells = DECell[ DECell(prob, Tsit5(), divide) ];

chem = Chemostats.simulate(cells, 2., Chemostats.Direct(); copy = true);

length(chem.pop)

[ cell.int[:V] for cell in first(chem.pop, 5) ]     # volumes of the first five cells

Chemostats.est_Λ(chem)

chem = Chemostats.simulate(cells, 20., Chemostats.Strict(20, false); copy = true)
length(chem.pop), Chemostats.est_Λ(chem)

chem = Chemostats.simulate(cells, 20., Chemostats.Strict(20, false); copy = true, save_mode = :ancestors);
```
<!-- END demo/ribosomes.jl -->

### See also

Chemostats.jl is designed for efficient estimation of growth rates in cell population models. If simple population simulations are required (corresponding to the `Direct` algorithm), the following packages may be useful:
* [Agents.jl](https://juliadynamics.github.io/Agents.jl/): A Julia framework for agent-based modelling.
* [AgentBasedModeling.jl](https://github.com/pihop/AgentBasedModeling.jl): A flexible framework built directly on top of Catalyst. 

### References 
* Giardinà, Kurchan, Peliti. "Direct evaluation of large-deviation functions," Phys. Rev. Lett. 96(12), 2006
* Lecomte & Tailleur. "A numerical approach to large deviations in continuous time," J. Stat. Mech.: Theory Exp. 2007(3), 2007

### Acknowledgments
* Thanks to [Chris Rackauckas](https://www.chrisrackauckas.com), [Aayush Sabharwal](https://github.com/AayushSabharwal), [Sam Isaacson](https://www.samisaacson.com), and the SciML community for their help and support (and their work on the SciML ecosystem!).
