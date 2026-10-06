using Chemostats
using Distributions
using OrdinaryDiffEqTsit5
using Random

# We represent each cell by a single variable `u`, the time to division.
# Then `u` follows a very simple ODE:
f(u, p, t) = -one(u)

# We divide when the time to division hits zero
div_cond(u, t, int) = u
cb_div = Chemostats.DivideCallback(div_cond)

# The time to division of every cell is exponentially distributed
u0 = randexp()

prob = ODEProblem(f, u0, (0., Inf); callback=cb_div)

function divide(int)
    # Return two cells with exponentially distributed times to division
    [ (u0 = randexp(), p = nothing), (u0 = randexp(), p = nothing) ]
end

cells = [ DECell(prob, Tsit5(), divide) ]
chem = Chemostat(cells)

# Estimate the growth rate from a single lineage. The true growth rate is
# exactly 1, but tracking only one lineage systematically underestimates it.
Chemostats.simulate!(chem, 1_000_000, Chemostats.Strict(1))
println("Estimated growth rate: ", Chemostats.est_Λ(chem))
