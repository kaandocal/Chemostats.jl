@enumx CellState Newborn Alive Dead Divided Killed

mutable struct Cell{I}
    anc::Union{Missing,Cell{I}}
    int::I
    state::CellState.T      # Move to integrator?
end

function simulate_cell! end 
function get_offspring end 

### MUST RENAME, cf. init_cell! below
get_curr_t(cell::Cell) = cell.int.t
get_state(cell::Cell) = cell.state

function set_state!(cell::Cell, state::CellState.T)
    savevalues!(cell)
    cell.state = state
end 

is_alive(cell::Cell) = get_state(cell) == CellState.Newborn || get_state(cell) == CellState.Alive

function init_cell!(cell::Cell)
    @argcheck cell.state == CellState.Newborn "Cannot reinitialise cell in state $(cell.state)"
    
    set_state!(cell, CellState.Alive)
end

function kill_cell!(cell::Cell, t=get_curr_t(cell)) 
    if t < get_curr_t(cell) || is_alive(cell)
        cell.state = CellState.Killed
    end 

    if SciMLBase.check_error(cell.int) == SciMLBase.ReturnCode.Default 
        terminate!(cell.int)
    end
end 

function die!(cell::Cell)
    if get_state(cell) != CellState.Alive
        @warn "Called `die!` on cell in state $(get_state(cell))"
        return 
    end 

    set_state!(cell, CellState.Dead)
    terminate!(cell.int)
end 

function divide!(cell::Cell)
    if get_state(cell) != CellState.Alive
        @warn "Cell in state $(get_state(cell)) tried to divide"
        return 
    end 

    set_state!(cell, CellState.Divided)
    terminate!(cell.int)
end

savevalues!(cell::Cell) = savevalues!(cell.int)

function step!(cell::Cell, dt, model, env)
    if get_state(cell) != CellState.Alive
        @warn "Tried to simulate cell in state $(get_state(cell))"
        return 
    end 

    step!(cell.int, dt, true)
    savevalues!(cell)
end


# function duplicate_cell end

# function interpolate(cell::Cell, t)
#     @argcheck cell.t[1] <= t <= cell.t[end]
#     i = findlast(s -> s <= t, cell.t)
#     cell.u[i]
# end

# function copy_until(cell::Cell, t)
#     @argcheck cell.t[1] <= t <= cell.t[end]
#     #@check cell.status == CellState.Alive || cell.status == CellState.Newborn

#     idcs = findlast(s -> s <= t, cell.t)
#     tt = [ cell.t[idcs]; t ]
#     uu = [ cell.u[idcs]; interpolate(cell, t) ]

#     Cell(cell.anc, copy(tt), deepcopy(uu), deepcopy(cell.p), CellState.Alive)
# end
