using OrdinaryDiffEq
using Chemostats

f_det(u, p, t) = u * p.Λ * log(2)

cb_det = Chemostats.DivideCallback((u, t, int) -> u - int.p.Vd; interp_points=0)

function divide_det(int)
    map(1:2) do _
        u0 = int.u / 2
        (; u0, p=nothing)
    end
end 

function Deterministic(Λ = log(2); Vd = 2.)
    p = (; Vd, Λ)
    u0 = p.Vd / 2
    ODEProblem(f_det, u0, (0., 0.), p; callback=cb_det)
end 
