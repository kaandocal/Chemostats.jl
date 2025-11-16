@enumx CellState Newborn Alive EndOfLife Dead Divided Killed

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

function DECell(prob; t0=0.) 
    if ismissing(get_divide_func(prob))
        @warn "DECell requires passing the `divide` kwarg to the AbstractDEProblem for cells to replicate"
    end 

    prob_ = remake(prob; kwargshandle=SciMLBase.KeywordArgSilent)
    DECell(missing, prob_; t0)
end

divide_nop(cell) = nothing 
function get_divide_func(prob)
    haskey(prob.kwargs, :divide) && return prob.kwargs[:divide]
    divide_nop
end

function get_offspring(int_, cell; save_lineages=false)
    anc = save_lineages ? cell : missing

    divide = get_divide_func(cell.int.sol.prob)
    offspring = divide(cell.int)
    isnothing(offspring) && return typeof(cell)[]

    map(offspring) do off
        @unpack u0, p = off
        prob = remake(cell.int.sol.prob; u0, p)
        Chemostats.DECell(anc, prob; t0=Chemostats.get_curr_t(cell))
    end
end

get_curr_t(cell::DECell) = cell.int.t
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

function kill_cell!(cell::DECell, t=get_curr_t(cell)) 
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
    @check cell.int.sol.t[1] <= t <= get_curr_t(cell)
    ret = deepcopy(cell)
    SciMLBase.change_t_via_interpolation!(ret.int, t, Val{true})
    ret
end 