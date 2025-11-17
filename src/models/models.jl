module Models 

include("deterministic.jl")
include("multitypemarkov.jl")
#include("sizecontrol.jl")

function get_Λ(prob::ODEProblem)
    if hasproperty(prob.p, :s) && hasproperty(prob.p, :model)
        get_Λ_mtm(prob.p.model)
    elseif length(prob.p) == 1 && hasproperty(prob.p, :λ)
        get_Λ_det(prob.p.λ)
    end
end 

end