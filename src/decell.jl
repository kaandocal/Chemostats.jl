@enumx CellState Newborn Alive Dead Divided Killed

mutable struct DECell{I}
    anc::Union{Missing,DECell{I}}
    int::I
    state::CellState.T      # Move to integrator?
end

function DECell(anc::Union{Missing, DECell}, prob; t0=ismissing(anc) ? 0. : get_curr_t(anc))
    int = init(prob; tspan=(t0, Inf))
    @check int.t == t0
    DECell(anc, int, CellState.Newborn)
end

DECell(prob; t0=0.) = DECell(missing, prob; t0)

get_curr_t(cell::DECell) = cell.int.t
get_state(cell::DECell) = cell.state

should_divide(cell::DECell) = SciMLBase.check_error(cell.int) == SciMLBase.ReturnCode.Terminated
function get_offspring end

function set_state!(cell::DECell, state::CellState.T)
    savevalues!(cell)
    cell.state = state
end 

is_alive(cell::DECell) = get_state(cell) == CellState.Newborn || get_state(cell) == CellState.Alive

function init_cell!(cell::DECell)
    @argcheck cell.state == CellState.Newborn "Cannot reinitialise cell in state $(cell.state)"
    
    set_state!(cell, CellState.Alive)
end

function kill_cell!(cell::DECell, t=get_curr_t(cell)) 
    if t < get_curr_t(cell) || is_alive(cell)
        cell.state = CellState.Killed
    end 

    if check_error(cell.int) == ReturnCode.Default 
        terminate!(cell.int)
    end
end 

function die!(cell::DECell)
    if get_state(cell) != CellState.Alive
        @warn "Called `die!` on cell in state $(get_state(cell))"
        return 
    end 

    set_state!(cell, CellState.Dead)
    terminate!(cell.int)
end 

function divide!(cell::DECell)
    if get_state(cell) != CellState.Alive
        @warn "Cell in state $(get_state(cell)) tried to divide"
        return 
    end 

    set_state!(cell, CellState.Divided)
    terminate!(cell.int)
end

savevalues!(cell::DECell) = savevalues!(cell.int)

function step!(cell::DECell, dt, model, env)
    if get_state(cell) != CellState.Alive
        @warn "Tried to simulate cell in state $(get_state(cell))"
        return 
    end 

    step!(cell.int, dt, true)
    savevalues!(cell)

    if should_divide(cell)
        divide!(cell)
    end
end

function DivideCallback(condition; kwargs...)
    SciMLBase.ContinuousCallback(condition, terminate!; kwargs...)
end

function duplicate_cell end

# function interpolate(cell::DECell, t)
#     @argcheck cell.t[1] <= t <= cell.t[end]
#     i = findlast(s -> s <= t, cell.t)
#     cell.u[i]
# end

# function copy_until(cell::DECell, t)
#     @argcheck cell.t[1] <= t <= cell.t[end]
#     #@check cell.status == CellState.Alive || cell.status == CellState.Newborn

#     idcs = findlast(s -> s <= t, cell.t)
#     tt = [ cell.t[idcs]; t ]
#     uu = [ cell.u[idcs]; interpolate(cell, t) ]

#     Cell(cell.anc, copy(tt), deepcopy(uu), deepcopy(cell.p), CellState.Alive)
# end
