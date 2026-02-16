mutable struct DECell{I,DF,FT}
    anc::Union{Missing,DECell}
    int::I
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
function DECell(anc::Union{Missing, DECell}, prob, divide; 
                t0=ismissing(anc) ? 0. : get_curr_t(anc), 
                reset_t=ismissing(anc) ? false : !isnothing(anc.tshift))
    tshift, int = if reset_t
        t0, init(prob; tspan=(0, Inf))
    else 
        nothing, init(prob; tspan=(t0, Inf))
    end 

    DECell(anc, int, divide, CellState.Newborn, tshift)
end

function DECell(prob, divide; kwargs...) 
    DECell(missing, prob, divide; kwargs...)
end

finalise_cell(cell::DECell) = cell

function reinit(old_int; t0, p=nothing, u0=nothing)
    int = init(old_int.sol.prob; tspan=(t0, Inf))

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

function get_offspring(int_, cell; save_lineages=false)
    anc = save_lineages ? finalise_cell(cell) : missing

    offspring = cell.divide(cell.int)
    isnothing(offspring) && return typeof(cell)[]

    t = get_curr_t(cell)

    map(offspring) do off 
        tshift, t0 = if isnothing(cell.tshift)
            nothing, t 
        else 
            t, zero(t)
        end 

        int = reinit(cell.int; t0, u0=off.u0, p=off.p)
        DECell(anc, int, cell.divide, CellState.Newborn, tshift)
    end
end

function get_curr_t(cell::DECell) 
    if isnothing(cell.tshift)
        cell.int.t
    else 
        cell.int.t + cell.tshift
    end 
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
    if t < get_curr_t(cell) || is_alive(cell)
        cell.state = CellState.Killed
    end 

    SciMLBase.done(cell.int) || terminate!(cell.int)
end 

function die!(cell::DECell)
    if get_state(cell) != CellState.EndOfLife
        @warn "Called `die!` on cell in state $(get_state(cell))"
        return 
    end 

    set_state!(cell, CellState.Dead)
end 

function divide!(cell::DECell)
    if get_state(cell) != CellState.EndOfLife
        @warn "Cell in state $(get_state(cell)) tried to divide"
        return 
    end 

    set_state!(cell, CellState.Divided)
end

savevalues!(cell::DECell) = savevalues!(cell.int)

function step!(cell::DECell, dt, p)
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
        else 
            set_state!(cell, CellState.EndOfLife)
        end
    end 
end

function duplicate_cell(cell::DECell, t)
    tshift = isnothing(cell.tshift) ? zero(t) : cell.tshift
    @check cell.int.sol.t[1] <= t - tshift <= get_curr_t(cell)

    int = reinit(cell.int; t0=t - tshift, u0=cell.int.sol(t - tshift), p=cell.int.p)
    DECell(cell.anc, int, cell.divide, CellState.Newborn, cell.tshift)
end 

###

"""
    DivideCallback(condition; kwargs...)

Implements a `DifferentialEquations.jl` callback to check whether a cell has reached the end of its lifetime. 
Internally, this returns a `ContinuousCallback` that calls `terminate!`. 
"""
function DivideCallback(condition; kwargs...)
    SciMLBase.ContinuousCallback(condition, SciMLBase.terminate!; kwargs...)
end
