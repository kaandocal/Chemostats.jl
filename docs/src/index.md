# Chemostats.jl Documentation

```@meta
CurrentModule=Chemostats
```

## Introduction

Chemostats.jl is a package to simulate branching populations in Julia, with a focus on cells and bacteria. It implements state-of-the-art algorithms to efficiently estimate population dynamics at sub-exponential cost. 

## Features:

- Efficient parallel simulation of branching processes
- Directly compatible with the [SciML](https://sciml.ai) ecosystem, including [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) and [Catalyst.jl](https://github.com/SciML/Catalyst.jl)
- Implements state-of-the-art algorithms (thinning, strict and lax cloning)
- Support for lineage analysis

This package is experimental - any feedback is appreciated (either by email or by opening an issue on GitHub).

## Example

A simple model of cell division, featuring volume growth coupled to ribosome production (the full, runnable script lives at `demo/ribosomes.jl`):

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "demo", "ribosomes.jl"), String) * "\n```")
```

### See also

Chemostats.jl is designed for efficient estimation of growth rates in cell population models. If simple population simulations are required (corresponding to the `Direct` algorithm), the following packages may be relevant:
* [Agents.jl](https://juliadynamics.github.io/Agents.jl/): A Julia framework for agent-based modelling.
* [AgentBasedModeling.jl](https://github.com/pihop/AgentBasedModeling.jl): A flexible framework built directly on top of Catalyst. 

### References 
* Giardinà, Kurchan, Peliti. "Direct evaluation of large-deviation functions," Phys. Rev. Lett. 96(12), 2006
* Lecomte & Tailleur. "A numerical approach to large deviations in continuous time," J. Stat. Mech.: Theory Exp. 2007(3), 2007

### Acknowledgments
* Thanks to [Chris Rackaukas](https://www.chrisrackauckas.com), [Aayush Sabharwal](https://github.com/AayushSabharwal), [Sam Isaacson](https://www.samisaacson.com), and the SciML community for their help and support (and their work on the SciML ecosystem!).
