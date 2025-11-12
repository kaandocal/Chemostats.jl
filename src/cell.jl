@enumx CellState Newborn Alive Dead Divided Killed

mutable struct Cell{C,P}
    anc::Union{Missing,Cell{C,P}}
    t::Vector{Float64}
    u::Vector{C}
    p::P
    status::CellState.T      # Move to integrator?
end

Cell(t::Real, u, p=nothing) = Cell(missing, t, u, p)
function Cell(::Missing, t::Real, u, p = nothing)
    Cell(missing, [t], [u], p, CellState.Newborn)
end

function Cell(anc::Cell{C,P}, t::Real, u::C, p::P = nothing) where {C,P} 
    #@assert anc.status == CellState.Divided

    Cell(anc, [t], [u], p, CellState.Newborn)
end

function simulate_cell! end 
function get_offspring end 

function interpolate(cell::Cell, t)
    @argcheck cell.t[1] <= t <= cell.t[end]
    i = findlast(s -> s <= t, cell.t)
    cell.u[i]
end

function copy_until(cell::Cell, t)
    @argcheck cell.t[1] <= t <= cell.t[end]
    #@check cell.status == CellState.Alive || cell.status == CellState.Newborn

    idcs = findlast(s -> s <= t, cell.t)
    tt = [ cell.t[idcs]; t ]
    uu = [ cell.u[idcs]; interpolate(cell, t) ]

    Cell(cell.anc, copy(tt), deepcopy(uu), deepcopy(cell.p), CellState.Alive)
end

### 

mutable struct CellIntegrator{I,C}
    int::I
    sol::C
end

### MUST RENAME, cf. init_cell! below
function init_cell end
function init_offspring end 

get_t(int::CellIntegrator) = int.int.t
get_status(int::CellIntegrator) = int.sol.status
is_alive(int::CellIntegrator) = get_status(int) == CellState.Newborn || get_status(int) == CellState.Alive

function init_cell!(int::CellIntegrator)
    int.sol.status = CellState.Alive
    savevalues!(int)
end

function set_status!(int::CellIntegrator, status::CellState.T)
    savevalues!(int)
    int.sol.status = status
end 

function kill_cell!(int::CellIntegrator, t=get_t(int)) 
    if t < get_t(int) || is_alive(int)
        int.sol.status = CellState.Killed
    end 

    if SciMLBase.check_error(int.int) == SciMLBase.ReturnCode.Default 
        terminate!(int.int)
    end
end 

function die!(int::CellIntegrator)
    set_status!(int, CellState.Dead)
    terminate!(int.int)
end 

function divide!(int::CellIntegrator)
    set_status!(int, CellState.Divided)
    terminate!(int.int)
end

function savevalues!(int::CellIntegrator)
    t = get_t(int)
    #@argcheck get_status(int) == CellState.Alive

    push!(int.sol.t, t)
    push!(int.sol.u, int.int.sol(t))

    return (true, true)
end

function step!(int::CellIntegrator, dt, model, env)
    @assert get_status(int) != CellState.Newborn
    get_status(int) == CellState.Alive || return
    step!(int.int, dt, true)
    savevalues!(int)
end

function duplicate_cell end 