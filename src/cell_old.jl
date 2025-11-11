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

### 

mutable struct CellIntegrator{C,P}
    t::Float64
    u::C
    p::P
    sol::Cell{C,P}
end

CellIntegrator(cell::Cell) = CellIntegrator(cell.t[end], cell.u[end], cell.p, cell)

### MUST RENAME
init_cell(cell::Cell) = CellIntegrator(cell)

function init_offspring(int, cell::CellIntegrator; save_lineages=false)
    offspr = get_offspring(cell, int.chem.model, int.chem.env)

    anc = if save_lineages 
        cell.sol 
    else 
        missing 
    end 

    map(((u, p),) -> CellIntegrator(Cell(anc, get_t(cell), u, p)), offspr)
end 

get_t(int::CellIntegrator) = int.t
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
end 

die!(int::CellIntegrator) = set_status!(int, CellState.Dead)
divide!(int::CellIntegrator) = set_status!(int, CellState.Divided)

function savevalues!(int::CellIntegrator)
    if get_t(int) in int.sol.t
        return (true, true)
    end 

    @argcheck get_status(int) == CellState.Alive

    push!(int.sol.t, int.t)
    push!(int.sol.u, int.u)

    return (true, true)
end

function step!(int::CellIntegrator, dt, model, env)
    @assert get_status(int) != CellState.Newborn
    get_status(int) == CellState.Alive || return
    simulate_cell!(int, get_t(int) + dt, model, env)
    savevalues!(int)
end

# function get_u(cell::Cell, t)
#     i = findfirst(isequal(t), cell.t)
#     @check !isnothing(i) "Cell state not saved at time t=$t"

#     cell.u[i]
# end

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

function duplicate_cell(int::CellIntegrator, t::Float64)
    CellIntegrator(copy_until(int.sol, t))
end