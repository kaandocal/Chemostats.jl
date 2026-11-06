using OrdinaryDiffEqTsit5
using Chemostats

f_det(u, p, t) = u * p.Λ * log(2)

cb_det = Chemostats.DivideCallback((u, t, int) -> u - int.p.Vd; interp_points=0)

divide_det = Chemostats.Offspring(2, CartesianIndex() => Chemostats.ExactPartition()) 

function Deterministic(Λ = log(2); Vd = 2.)
    p = (; Vd, Λ)
    u0 = p.Vd / 2
    ODEProblem(f_det, u0, (0., 0.), p; callback=cb_det, alg=Tsit5(), abstol=1e-9, reltol=1e-6)
end 
