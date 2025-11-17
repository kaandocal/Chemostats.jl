using OrdinaryDiffEq 

f_det(u, p, t) = u * p.λ
get_Λ_det(prob) = prob.p.λ

cb_det = ContinuousCallback((u, t, int) -> u - int.p.Vd, terminate!)

function divide_det(int)
    map(1:2) do _
        u0 = int.u / 2
        (; u0, p=int.p)
    end
end 

function Deterministic(; Vd = 2., λ = log(2))
    p = (; Vd, λ)
    u0 = p.Vd / 2
    ODEProblem(f, u0, (0., 0.), p; callback=cb, divide=divide_det)
end 