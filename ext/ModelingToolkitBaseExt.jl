module ModelingToolkitBaseExt

using Chemostats
using ModelingToolkitBase
using SciMLBase

MTK = ModelingToolkitBase

function affect_terminate!(x, obs, ctx, int)
    terminate!(int)
    x
end 

Chemostats.MTKDivideAffect() = MTK.ImperativeAffect(affect_terminate!)

function Chemostats.SymbolicDivideContinuous(cond; iv=MTK.t_nounits, name=:__chemostat_div__)
    event = MTK.SymbolicContinuousCallback(cond, Chemostats.MTKDivideAffect())
    MTK.System(MTK.Equation[], iv; name, continuous_events = [ event ])
end 

function Chemostats.SymbolicDivideDiscrete(cond; iv=MTK.t_nounits, name=:__chemostat_div__)
    event = MTK.SymbolicDiscreteCallback(cond, Chemostats.MTKDivideAffect())
    MTK.System(MTK.Equation[], iv; name, discrete_events = [ event ])
end

"""
    DivideCallback(eqs::MTK.Equation; kwargs...)
    DivideCallback(eqs::AbstractVector{<:MTK.Equation}; kwargs...)

Create a [`DivideCallback`](@ref) from one or more symbolic equalities of the form `lhs ~ rhs`. 
"""
function Chemostats.DivideCallback(eqs::Union{MTK.Equation, AbstractVector{<:MTK.Equation}}; kwargs...)
    eqs = eqs isa MTK.Equation ? [eqs] : eqs
    exprs = [eq.lhs - eq.rhs for eq in eqs]
    expr = length(exprs) == 1 ? only(exprs) : exprs

    getter = Ref{Any}(nothing)

    function condition(u, t, int)
        isnothing(getter[]) && (getter[] = MTK.getu(int, expr))
        getter[](MTK.ProblemState(; u, p = int.p, t))
    end

    Chemostats.DivideCallback(condition; kwargs...)
end

end