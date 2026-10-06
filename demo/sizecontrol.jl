using Chemostats
using Distributions
using OrdinaryDiffEqTsit5

# We represent each cell by its size `u`, which grows exponentially in time:
f(u, p, t) = u

# We divide once a cell hits its division volume `V_d`, stored as a parameter
div_cond(u, t, int) = u - int.p.V_d
cb_div = Chemostats.DivideCallback(div_cond)

# The division threshold is the birth volume plus a Gamma-distributed amount;
# daughters inherit a random fraction of the parent's volume.
dist_Δ = Gamma(5, 0.2)
dist_div = Beta(1)

u0 = 1.
p = (V_d = 2.,)

prob = ODEProblem(f, u0, (0., Inf), p; callback=cb_div)

function divide(int)
    # Sample sizes of each daughter cell
    f = rand(dist_div)
    V1 = f * int.u
    V2 = (1 - f) * int.u

    # Return two cells with the given starting sizes and randomly sampled
    # division thresholds
    [ (u0 = V1, p = (V_d = V1 + rand(dist_Δ),)),
      (u0 = V2, p = (V_d = V2 + rand(dist_Δ),)) ]
end

cells = [ DECell(prob, Tsit5(), divide) ]
chem = Chemostat(cells)

# The true growth rate equals the physical growth rate of each cell, Λ = 1.
Chemostats.simulate!(chem, 1000., Chemostats.Strict(10))
println("Estimated growth rate: ", Chemostats.est_Λ(chem))
