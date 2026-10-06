using Random
using Distributions
using OrdinaryDiffEqTsit5
using Chemostats

# The "adder" model: size grows exponentially, a cell divides once it has
# added Δ ~ Gamma(5, 0.2) (mean 1) to its birth size, and daughters inherit a
# random fraction f ~ Beta(1,1) of the parent's volume. The true growth rate
# equals the physical growth rate of each cell, Λ = 1 (see docs/src/theory.md
# and demo/sizecontrol.jl, of which this is a test-suite simplification).
f_size(u, p, t) = u

cb_size = Chemostats.DivideCallback((u, t, int) -> u - int.p.V_d)

const dist_Δ_size = Gamma(5, 0.2)
const dist_frac_size = Beta(1)

function divide_size(int)
    f = rand(dist_frac_size)
    V1 = f * int.u
    V2 = (1 - f) * int.u

    [ (u0 = V1, p = (V_d = V1 + rand(dist_Δ_size),)),
      (u0 = V2, p = (V_d = V2 + rand(dist_Δ_size),)) ]
end

function SizeControlCell(alg = Tsit5(); u0 = 1., V_d = 2.)
    prob = ODEProblem(f_size, u0, (0., 0.), (; V_d); callback = cb_size)
    DECell(prob, alg, divide_size)
end

make_sizecontrol_population(K; alg = Tsit5()) = [ SizeControlCell(alg) for _ in 1:K ]
