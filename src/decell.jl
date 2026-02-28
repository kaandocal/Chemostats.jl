mutable struct DECell{I,DF,FT}
    int::Union{I,Nothing}
    sol
    divide::DF
    state::CellState.T
    tshift::FT
end

"""
    DECell([anc = missing, ]prob, divide; reset_t=false)

Creates a new cell around a `SciMLBase.DEProblem`. The function `divide` takes a single integrator argument 
`int` and returns a list of offspring as a vector of `NamedTuples` with fields `u0` and `p`. The argument 
`reset_t` determines whether the integration time is started from `0` for every cell, which can avoid floating-point 
errors for long simulations (many generations). If `reset_t` is set to `true`, each cell should keep track of its 
starting time, e.g. via a parameter `t0` in `p`.

**Note:** The problem should include a callback that calls `terminate!` when a cell reaches the end of its life.
"""
function DECell(prob, alg, divide; t0 = 0., reset_t=false)
    tshift, int = if reset_t
        t0, init(prob, alg; tspan=(0, Inf))
    else 
        nothing, init(prob, alg; tspan=(t0, Inf))
    end 

    DECell(int, missing, divide, CellState.Newborn, tshift)
end

process_int(int) = int.sol 

function finalise!(cell::DECell)
    isnothing(cell.int) && return 

    cell.sol = process_int(cell.int)
    cell.int = nothing
end 

function reinit(old_int; t0, p=nothing, u0=nothing)
    int = init(old_int.sol.prob, old_int.alg; tspan=(t0, Inf))

    if !isnothing(p)
        if p isa AbstractVector{<:Pair}
            setp = SciMLBase.setp(int, first.(p))
            setp(int, last.(p))
        else 
            int.p = p 
        end 
    end 

    u = old_int.u 

    if !isnothing(u0) && !isempty(u0)
        if first(u0) isa Pair 
            setu = SciMLBase.setu(int, first.(u0))
            setu(u, last.(u0))
        else 
            u = u0
        end
    end

    SciMLBase.reinit!(int, u; reinit_dae=false)

    int
end 

function get_children(parent, p=nothing)
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

        int = reinit(parent.int; t0, u0=ch.u0, p=ch.p)
        DECell(int, missing, parent.divide, CellState.Newborn, tshift)
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
    
    if t < get_curr_t(cell) || is_alive(cell)
        cell.state = CellState.Killed
    end 

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

    isnothing(p) || @warn "Parameter arguments to DECell are ignored"

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

function clone_cell(cell::DECell, t)
    @assert !isnothing(cell.int)

    tshift = isnothing(cell.tshift) ? zero(t) : cell.tshift
    @check cell.int.sol.t[1] <= t - tshift <= get_curr_t(cell)

    int = reinit(cell.int; t0=t - tshift, u0=cell.int.sol(t - tshift), p=cell.int.p)
    DECell(int, missing, cell.divide, CellState.Newborn, cell.tshift)
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
