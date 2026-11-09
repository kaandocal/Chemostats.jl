# Usage

```@meta
CurrentModule = Chemostats
```

In this section, we show how to build a simple model of cell growth using `Catalyst.jl` and `Chemostats.jl`. In our
model, every cell has an exponentially growing volume `V`, and produces ribosomes `R` at random, following a jump
process. The rate of volume growth and the rate of ribosome production are coupled. The full script is in
`demo/ribosomes.jl`.

````julia
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
````

## Division

To enable division, we have to specify when a cell should divide. We can do this with a `DivideCallback`
that triggers once the cell has reached a fixed division volume `V_d`:

````julia
cb_div = Chemostats.DivideCallback(rn.V ~ rn.V_d);
````

We next have to tell `Chemostats.jl` what happens upon division. Our example is very simple: with probability 1/10,
the cell dies; otherwise, each daughter gets half of the volume, and the ribosomes are partitioned binomially:

````julia
function divide(int)
    rand() < 0.1 && return nothing

    Chemostats.Offspring(2, :R => Chemostats.BinomialPartition(),
                            :V => Chemostats.ExactPartition())
end;
````

## Simulation

Each cell is a `DECell`, which takes a `DifferentialEquations.jl` problem, a solver, and a `divide`
function. The division callback `cb_div` is part of the `HybridProblem`.

````julia
u0 = [ :R => 0., :V => 1. ]
p = [ :σ => 10., :λ => 1., :V_d => 2. ]
prob = HybridProblem(rn, u0, (0., Inf), p; callback = cb_div)

cells = [ DECell(prob, Tsit5(), divide) ];
````

We simulate the full population in a `Chemostat` until time 2, using the `Direct` algorithm. Since
`simulate` modifies the cells, we pass `copy = true` to reuse `cells` below:

````julia
chem = Chemostats.simulate(cells, 2., Chemostats.Direct(); copy = true);
````

We can then view the cells that are currently alive using `chem.pop`:

````julia
length(chem.pop)
````

````
108
````

````julia
[ cell.int[:V] for cell in first(chem.pop, 5) ]     # volumes of the first five cells
````

````
5-element Vector{Float64}:
 1.0423164542929497
 1.016818323022799
 1.0084091615113993
 1.0300651269626537
 1.1075921756578564
````

We can also estimate the growth rate of the population using `Chemostats.est_growth_rate`:

````julia
Chemostats.est_growth_rate(chem)
````

````
2.34106561356211
````

Note that our population grows exponentially, which makes long-term simulations difficult. To avoid this, we can
use the `Strict` or `Lax` simulation algorithms, which maintain a (roughly) constant population
size. For example, we can simulate a population consisting of at most 20 cells for ~100 generations:

````julia
chem = Chemostats.simulate(cells, 20., Chemostats.Strict(20, false); copy = true)
length(chem.pop), Chemostats.est_growth_rate(chem)
````

````
(20, 2.5697787564815346)
````

This allows us to compute much better estimates of population quantities like the growth rate. See
Simulation Algorithms for the available algorithms.

Chemostats can also record the lineage tree during a simulation, using `save_mode = :ancestors` or
`save_mode = :full`:

````julia
chem = Chemostats.simulate(cells, 20., Chemostats.Strict(20, false); copy = true, save_mode = :ancestors);
````

The resulting population tree is stored in `chem.tree` and can be analysed using the [Tree API](../trees/index.md).
See the [examples](growth_protein.md) for how to do this.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

