module ModelingToolkitBaseExt

using Chemostats
using ModelingToolkitBase
using SciMLBase

const MTK = ModelingToolkitBase

# MTK compilation of getu does not seem to be threadsafe, use global lock for now
const _GETU_LOCK = ReentrantLock()

# Accessor shared by all copies of a callback. Avoid deepcopy when callbacks 
# are copied upon division
mutable struct SharedGetter
    @atomic getter::Any
end

Base.deepcopy_internal(g::SharedGetter, ::IdDict) = g

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
    event = MTK.SymbolicDiscreteCallback(cond => Chemostats.MTKDivideAffect())
    MTK.System(MTK.Equation[], iv; name, discrete_events = [ event ])
end

"""
    DivideCallback(eqs::MTK.Equation; discrete=false, kwargs...)
    DivideCallback(eqs::AbstractVector{<:MTK.Equation}; discrete=false, kwargs...)

Create a [`DivideCallback`](@ref) from one or more symbolic equalities of the form `lhs ~ rhs`.
`discrete` determines whether to use a `DiscreteCallback` or a `ContinuousCallback` (defaults to `false`)
"""
function Chemostats.DivideCallback(eqs::Union{MTK.Equation, AbstractVector{<:MTK.Equation}}; discrete::Bool=false, kwargs...)
    eqs = eqs isa MTK.Equation ? [eqs] : eqs
    exprs = [eq.lhs - eq.rhs for eq in eqs]
    expr = length(exprs) == 1 ? only(exprs) : exprs

    getter = SharedGetter(nothing)

    function condition(u, t, int)
        # Lock only while building the accessor, not on every call
        g = @atomic getter.getter
        if isnothing(g)
            g = lock(_GETU_LOCK) do
                isnothing(@atomic getter.getter) && (@atomic getter.getter = MTK.getu(int, expr))
                @atomic getter.getter
            end
        end
        val = g(MTK.ProblemState(; u, p = int.p, t))
        discrete ? all(iszero, val) : val
    end

    Chemostats.DivideCallback(condition; discrete, kwargs...)
end

end