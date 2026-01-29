module Models 

include("deterministic.jl")
include("multitypemarkov.jl")
#include("sizecontrol.jl")

function sample_cell(prob::ODEProblem)
    if hasproperty(prob.p, :s) && hasproperty(prob.p, :model)
        sample_cell_mtm(prob)
    elseif length(prob.p) == 1 && hasproperty(prob.p, :λ)
        DECell(prob, divide_det)
    end
end 

function get_Λ(prob::ODEProblem)
    if hasproperty(prob.p, :s) && hasproperty(prob.p, :model)
        get_Λ_mtm(prob.p.model)
    elseif hasproperty(prob.p, :Vd) && hasproperty(prob.p, :λ)
        get_Λ_det(prob)
    end
end 

end