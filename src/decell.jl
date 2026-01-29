@enumx CellState Newborn Alive EndOfLife Dead Divided Killed

mutable struct DECell{I,DF}
    anc::Union{Missing,DECell{I,DF}}
    int::I
    divide::DF
    state::CellState.T
end

function DECell(anc::Union{Missing, DECell}, prob, divide; t0=ismissing(anc) ? 0. : get_curr_t(anc))
    int = init(prob; tspan=(t0, Inf))
    DECell(anc, int, divide, CellState.Newborn)
end

function DECell(prob, divide) 
    DECell(missing, prob, divide)
end

finalise_cell(cell::DECell) = cell

function remake_cell(cell::DECell; anc=cell.anc, t=get_curr_t(cell), u0=nothing, p=nothing)
    int = init(cell.int.sol.prob; tspan=(t, Inf))

    if !isnothing(p)
        if p isa AbstractVector{<:Pair}
            setp = ModelingToolkit.setp(prob, first.(p))
            setp(int, last.(p))
        else 
            int.p = p 
        end 
    end 

    u = int.u 

    if !isnothing(u0) && !isempty(u0)
        if first(u0) isa Pair 
            setu = ModelingToolkit.setu(prob, first.(u0))
            setu(u, last.(u0))
        else 
            u = u0
        end
    end

    SciMLBase.reinit!(int, u; reinit_dae=false)

    DECell(anc, int, cell.divide, CellState.Newborn)
end 

function get_offspring(int_, cell; save_lineages=false)
    anc = save_lineages ? finalise_cell(cell) : missing

    offspring = cell.divide(cell.int)
    isnothing(offspring) && return typeof(cell)[]

    [ remake_cell(cell; anc, u0=off.u0, p=off.p) for off in offspring ]
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

    remake_cell(cell; t, u0=cell.int.sol(t))
end 