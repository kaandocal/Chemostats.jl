mutable struct DECell{I,DF,FT,P,Remake}
    int::Union{I,Nothing}
    prob::P
    sol
    divide::DF
    state::CellState.T
    tshift::FT
end

"""
    DECell(prob, alg, divide; t0 = 0., reset_t=false, remake::Bool=false)

Creates a new cell around a `SciMLBase.DEProblem`, simulated with algorithm `alg` starting at time `t0`.
The function `divide` takes a single integrator argument `int` and returns a list of offspring as a vector of `NamedTuples` with fields `u0` and `p`.

`reset_t` determines whether the integration time is started from `0` for every cell, which can avoid floating-point
errors for long simulations (many generations). If `reset_t` is set to `true`, each cell should keep track of its
starting time, e.g. via a parameter `t0` in `p`.

If `remake` is `true`, descendants will be created using `SciMLBase.remake`, which creates 
a copy of the underlying `AbstractProblem` (safest). Otherwise, descendants will be 
created using `SciMLBase.reinit`, which should be faster for `ODEProblem` and `SDEProblem` types.

**Note:** The problem should include a callback that calls `terminate!` when a cell reaches the end of its life.
"""
function DECell(prob, alg, divide; t0 = 0., reset_t=false, remake::Bool=false)
    tshift, int = if reset_t
        t0, copy_integrator(prob, alg, 0)
    else
        nothing, copy_integrator(prob, alg, t0)
    end

    # `reinit_remake` never reads `cell.prob`, so don't keep it alive for it.
    stored_prob = remake ? nothing : prob

    DECell{typeof(int), typeof(divide), typeof(tshift), typeof(stored_prob), remake}(
        int, stored_prob, missing, divide, CellState.Newborn, tshift
    )
end

# `DiscreteCallback` can hold mutable state
copy_callback(cb) = isempty(cb.discrete_callbacks) ? cb : deepcopy(cb)

# Copy to avoid cases where `init` aliases `p`
function copy_integrator(prob, alg, t0)
    init_kwargs = prob isa SciMLBase.AbstractJumpProblem ? (; alias_jump=false) : (;)
    int = init(prob, alg; tspan=(t0, Inf), init_kwargs...)
    int.p = deepcopy(int.p)
    int
end

process_int(int) = int.sol

function finalise!(cell::DECell)
    isnothing(cell.int) && return 

    cell.sol = process_int(cell.int)
    cell.int = nothing
end 

function reinit(old_int, prob; t0, p=nothing, u0=nothing)
    int = copy_integrator(prob, old_int.alg, t0)
    apply_reinit_overrides!(int, old_int; p, u0)
end

function reinit_remake(old_int; t0, p=nothing, u0=nothing)
    base_prob = SciMLBase.remake(
        old_int.sol.prob; u0 = copy(old_int.u), p = old_int.p, tspan = (t0, Inf)
    )

    int = init(base_prob, old_int.alg; callback = copy_callback(old_int.opts.callback))
    int.p = deepcopy(int.p)
    apply_reinit_overrides!(int, old_int; p, u0)
end

function apply_reinit_overrides!(int, old_int; p, u0)
    if !isnothing(p)
        if p isa AbstractVector{<:Pair}
            setp = SciMLBase.setp(int, first.(p))
            setp(int, last.(p))
        else
            int.p = deepcopy(p)
        end
    end

    u = copy(old_int.u)

    if !isnothing(u0) && !isempty(u0)
        if first(u0) isa Pair
            setu = SciMLBase.setu(int, first.(u0))
            setu(u, last.(u0))
        else
            u = copy(u0)
        end
    end

    SciMLBase.reinit!(int, u; reinit_dae=false)

    int
end

# Dispatches on a cell's `remake` type parameter; `prob` is ignored when `remake=true`.
reinit_dispatch(remake::Bool, old_int, prob; t0, p=nothing, u0=nothing) =
    remake ? reinit_remake(old_int; t0, p, u0) : reinit(old_int, prob; t0, p, u0)

function get_children(parent::DECell{I,DF,FT,P,Remake}, p=nothing) where {I,DF,FT,P,Remake}
    @assert !isnothing(parent.int)

    children = parent.divide(parent.int)
    isnothing(children) && return typeof(parent)[]

    t = get_curr_t(parent)

    map(children) do ch
        tshift, t0 = if isnothing(parent.tshift)
            nothing, t
        else
            t, zero(t)
        end

        int = reinit_dispatch(Remake, parent.int, parent.prob; t0, u0=ch.u0, p=ch.p)
        DECell{typeof(int), DF, typeof(tshift), P, Remake}(
            int, parent.prob, missing, parent.divide, CellState.Newborn, tshift
        )
    end
end

function get_curr_t(cell::DECell) 
    ret = if !isnothing(cell.int)
        cell.int.t 
    else 
        last(cell.sol.t)
    end 

    if !isnothing(cell.tshift)
        ret += cell.tshift 
    end 

    ret 
end 

get_state(cell::DECell) = cell.state

function set_state!(cell::DECell, state::CellState.T)
    savevalues!(cell)
    cell.state = state
end 

is_alive(cell::DECell) = get_state(cell) == CellState.Newborn || get_state(cell) == CellState.Alive

function init_cell!(cell::DECell)
    @argcheck cell.state == CellState.Newborn "Cannot reinitialise cell in state $(cell.state)"
    
    set_state!(cell, CellState.Alive)
end

function kill!(cell::DECell, t=get_curr_t(cell))
    @assert !isnothing(cell.int)

    cell.state = CellState.Killed

    SciMLBase.done(cell.int) || terminate!(cell.int)
    finalise!(cell)
end

function die!(cell::DECell)
    @assert !isnothing(cell.int)

    if get_state(cell) != CellState.EndOfLife
        @warn "Called `die!` on cell in state $(get_state(cell))"
        return 
    end 

    set_state!(cell, CellState.Dead)
    finalise!(cell)
end 

function divide!(cell::DECell)
    @assert !isnothing(cell.int)

    if get_state(cell) != CellState.EndOfLife
        @warn "Cell in state $(get_state(cell)) tried to divide"
        return 
    end 

    set_state!(cell, CellState.Divided)
    finalise!(cell)
end

savevalues!(cell::DECell) = savevalues!(cell.int)

function step!(cell::DECell, dt, p)
    @assert !isnothing(cell.int)

    isnothing(p) || @warn "Parameter arguments to DECell are ignored" maxlog=1

    if get_state(cell) != CellState.Alive
        @warn "Tried to simulate cell in state $(get_state(cell))"
        return 
    end 

    step!(cell.int, dt, true)
    savevalues!(cell)

    if SciMLBase.done(cell.int)
        if !SciMLBase.successful_retcode(cell.int.sol)
            @warn "Cell solver errored, removing cell..."
            set_state!(cell, CellState.Killed)
            finalise!(cell)
        else 
            set_state!(cell, CellState.EndOfLife)
        end
    end 
end

function clone_cell(cell::DECell{I,DF,FT,P,Remake}, t) where {I,DF,FT,P,Remake}
    @assert !isnothing(cell.int)

    tshift = isnothing(cell.tshift) ? zero(t) : cell.tshift
    @check cell.int.sol.t[1] <= t - tshift <= get_curr_t(cell)

    int = reinit_dispatch(Remake, cell.int, cell.prob; t0=t - tshift, u0=cell.int.sol(t - tshift), p=cell.int.p)
    DECell{typeof(int), DF, FT, P, Remake}(int, cell.prob, missing, cell.divide, CellState.Newborn, cell.tshift)
end

###

"""
    DivideCallback(condition; kwargs...)

Implements a differential equation callback to check whether a cell has reached the end of its lifetime. 
Internally, this returns a `ContinuousCallback` that calls `terminate!`. 
"""
function DivideCallback(condition; kwargs...)
    SciMLBase.ContinuousCallback(condition, SciMLBase.terminate!; kwargs...)
end
