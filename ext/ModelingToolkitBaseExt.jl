module ModelingToolkitBaseExt

using Chemostats
using ModelingToolkitBase
using SciMLBase
using SymbolicIndexingInterface: variable_symbols, parameter_symbols

const MTK = ModelingToolkitBase

# MTK compilation of getu does not seem to be threadsafe, use global lock for now
const _GETU_LOCK = ReentrantLock()

# Accessor shared by all copies of a callback. Avoid deepcopy when callbacks
# are copied upon division
mutable struct SharedGetter
    @atomic getter::Any     # (accessor, layout hash) once built
end

Base.deepcopy_internal(g::SharedGetter, ::IdDict) = g

# Variable/parameter layout an accessor is compiled for; unchanged by deepcopy
layout_hash(int) = hash(variable_symbols(int), hash(parameter_symbols(int)))

# Function barriers keep the hot path in `condition` small and type-stable
@noinline function build_getter!(getter, expr, int)
    lock(_GETU_LOCK) do
        if isnothing(@atomic getter.getter)
            g = MTK.getu(int, expr)
            g(MTK.ProblemState(; u = int.u, p = int.p, t = int.t))    # compile under the lock
            @atomic getter.getter = (g, layout_hash(int))
        end
        @atomic getter.getter
    end
end

@noinline function check_layout!(checked, entry, int)
    last(entry) == layout_hash(int) || error("DivideCallback cannot be shared between different systems")
    checked[] = true
end

@noinline eval_getter(g, u, p, t) = g(MTK.ProblemState(; u, p, t))

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
    discrete && eqs isa AbstractVector && throw(ArgumentError("Discrete DivideCallback takes a single equation"))
    eqs = eqs isa MTK.Equation ? [eqs] : eqs
    exprs = [eq.lhs - eq.rhs for eq in eqs]
    expr = length(exprs) == 1 ? only(exprs) : exprs

    getter = SharedGetter(nothing)
    checked = Ref(false)    # Copied with the callback, so each copy checks once

    function condition(u, t, int)
        entry = @atomic getter.getter
        isnothing(entry) && (entry = build_getter!(getter, expr, int))
        checked[] || check_layout!(checked, entry, int)
        eval_getter(first(entry), u, int.p, t)
    end

    if discrete
        Chemostats.DivideCallback((u, t, int) -> iszero(condition(u, t, int)); discrete, kwargs...)
    elseif length(exprs) == 1
        Chemostats.DivideCallback(condition; kwargs...)
    else
        # Divide once any of the equations holds
        condition!(out, u, t, int) = (out .= condition(u, t, int))
        SciMLBase.VectorContinuousCallback(condition!, (int, _) -> terminate!(int), length(exprs); kwargs...)
    end
end

end