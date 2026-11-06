# Chemostats.jl

```@meta
CurrentModule=Chemostats
```

## Introduction

Chemostats.jl is a package to simulate branching populations in Julia, with a focus on cells and bacteria. It implements state-of-the-art algorithms to efficiently estimate population dynamics at sub-exponential cost, based on the paper

> Kaan Öcal, Lars Maas, Michael P.H. Stumpf. **"Probing bacterial fitness in virtual chemostats."** [*bioRxiv*](https://www.biorxiv.org/content/10.64898/2026.10.31.763286v1) (2026).

This package is somewhat experimental - any feedback is appreciated (either by email or by opening an issue on GitHub).

To get started, have a look at the [Usage](@ref) section or browse the examples!

## Features

- Efficient parallel simulation of branching processes
- Directly compatible with the [SciML](https://sciml.ai) ecosystem, including [DifferentialEquations.jl](@extref DifferentialEquations :doc:`index`), [ModelingToolkit.jl](@extref ModelingToolkit :doc:`index`) and [Catalyst.jl](@extref Catalyst :doc:`index`)
- Support for basic cell lineage analysis

### See also

Chemostats.jl is designed for efficient estimation of growth rates in cell population models. If simple population simulations are required (corresponding to the `Direct` algorithm), the following packages may be relevant:
* [Agents.jl](https://juliadynamics.github.io/Agents.jl/): A Julia framework for agent-based modelling.
* [AgentBasedModeling.jl](https://github.com/pihop/AgentBasedModeling.jl): A flexible framework built directly on top of Catalyst. 

### References 
* Giardinà, Kurchan, Peliti. "Direct evaluation of large-deviation functions," Phys. Rev. Lett. 96(12), 2006
* Lecomte & Tailleur. "A numerical approach to large deviations in continuous time," J. Stat. Mech.: Theory Exp. 2007(3), 2007

### Acknowledgments
* Thanks to [Chris Rackauckas](https://www.chrisrackauckas.com), [Aayush Sabharwal](https://github.com/AayushSabharwal), [Sam Isaacson](https://www.samisaacson.com), and the SciML community for their help and support (and their work on the SciML ecosystem!).
