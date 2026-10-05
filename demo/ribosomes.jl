using Catalyst
using Chemostats
using Distributions
using OrdinaryDiffEqTsit5
using CairoMakie

# Define a reaction network in which ribosomes are produced in a volume-dependent manner.
# Time is measured in arbitrary units.
rn = @reaction_network begin
    @species R(t) = 0.
    @variables V(t) = 1.
    @parameters σ = 10. λ = 1. V_d = 2.

    σ * V, 0 --> R

    @equations begin
        D(V) ~ λ * R
    end
end

# A cell divides once it doubles its initial volume.
cb_div = Chemostats.DivideCallback(rn.V ~ rn.V_d)

# This function determines what offspring a dividing cell produces.
function divide(int)
    rand() < 0.1 && return nothing      # Die with probability 1/10

    R = round(Int, int[:R])
    V = int[:V]

    # Create two cells with half the volume and (roughly) half the ribosomes.
    # The model parameters are inherited from the parent cell.
    R1 = rand(Binomial(R, 0.5))
    R2 = R - R1

    [ (u0 = [:R => R1, :V => V / 2], p = nothing),
      (u0 = [:R => R2, :V => V / 2], p = nothing) ]
end

u0 = [ :R => 0., :V => 1. ]
p = [ :σ => 10., :λ => 1., :V_d => 2. ]
prob = HybridProblem(rn, u0, (0., Inf), p; callback = cb_div)

# Need a DECell array for type instability reasons (to be fixed)
cells = DECell[ DECell(prob, Tsit5(), divide) ] 
chem = Chemostat(cells)

# Simulate population for 10 time units
Chemostats.simulate!(chem, 10., Chemostats.Strict(10); saveat=0:0.1:10.)

# The full population would require simulating a lot of cells!
println("Estimated growth rate: ", Chemostats.est_Λ(chem))