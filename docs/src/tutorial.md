# From a Single Cell to a Population

If you already use Catalyst or ModelingToolkit, you have something like this:

```julia
rn = @reaction_network begin
    @species R(t) = 0.
    @variables V(t) = 1.
    @parameters σ = 10. λ = 1. V_d = 2.

    σ * V, 0 --> R

    @equations begin
        D(V) ~ λ * R
    end
end
```

Ribosomes `R` are produced at a rate proportional to the cell's volume `V`, which in turn grows in proportion to the number of ribosomes. This is a perfectly ordinary reaction network - nothing about it is Chemostats-specific. Turning it into a population model takes three additions.

## 1. When does a cell divide?

A [`Chemostats.DivideCallback`](@ref) is a callback that stops the integrator once some condition holds. Here, a cell divides once it doubles its initial volume:

```julia
cb_div = Chemostats.DivideCallback(rn.V ~ rn.V_d)
```

It's passed to the problem like any other callback - `HybridProblem(rn, u0, tspan, p; callback = cb_div)`. By default the condition is treated as continuous (root-found as it crosses zero); for conditions that only make sense at discrete events, such as a counter hitting an integer, pass `discrete=true` and give a boolean expression instead of an equation.

## 2. What happens when it divides?

A plain function takes the integrator at the point of division and returns the daughter cells:

```julia
function divide(int)
    rand() < 0.1 && return nothing      # Die with probability 1/10

    R = round(Int, int[:R])
    V = int[:V]

    # Create two cells with half the volume and (roughly) half the ribosomes.
    R1 = rand(Binomial(R, 0.5))
    R2 = R - R1

    [ (u0 = [:R => R1, :V => V / 2], p = nothing),
      (u0 = [:R => R2, :V => V / 2], p = nothing) ]
end
```

Each daughter is a `NamedTuple` with `u0` and `p`; `p = nothing` means the daughter inherits the parent's parameters unchanged. Returning `nothing` instead of an array kills the cell instead - here with probability 1/10, standing in for whatever death condition your model needs.

## 3. Wrap it up

A [`Chemostats.DECell`](@ref) ties the problem, solver and `divide` function together; a [`Chemostat`](@ref) holds the population:

```julia
u0 = [ :R => 0., :V => 1. ]
p = [ :σ => 10., :λ => 1., :V_d => 2. ]
prob = HybridProblem(rn, u0, (0., Inf), p; callback = cb_div)

cells = DECell[ DECell(prob, Tsit5(), divide) ]
chem = Chemostat(cells)
```

Then simulate:

```julia
Chemostats.simulate!(chem, 10., Chemostats.Strict(10); saveat=0:0.1:10.)
println("Estimated growth rate: ", Chemostats.est_Λ(chem))
```

`Strict(10)` tracks a representative subpopulation of around 10 cells rather than the full, exponentially growing tree, while still giving an unbiased estimate of the growth rate. See [Simulation Algorithms](@ref) for what the other algorithms trade off against `Strict`.

The full script lives at `demo/ribosomes.jl`.
