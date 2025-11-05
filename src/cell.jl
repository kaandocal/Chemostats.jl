mutable struct Cell{C,P}
    anc::Union{Missing,Cell{C,P}}
    t::Vector{Float64}
    u::Vector{C}
    p::P
    alive::Bool      # Move to integrator?
end

Cell(t::Real, u, p=nothing) = Cell(missing, [t], [u], p, true)
Cell(anc::Union{Missing,Cell{C,P}}, t::Real, u::C, p::P = nothing) where {C,P} = Cell(anc, [t], [u], p, true)

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

function kill_cell!(int::CellIntegrator, t=int.t)
    savevalues!(int)
    int.sol.alive = false
end 

function savevalues!(int::CellIntegrator)
    if int.t in int.sol.t
        return (true, true)
    end 

    if !int.sol.alive
        @show int 
    end 

    @argcheck int.sol.alive

    push!(int.sol.t, int.t)
    push!(int.sol.u, int.u)

    return (true, true)
end

# function cell_alive(cell::Cell, t)
#     t >= cell.t[1] || return false
#     cell.alive && @check t <= cell.t[end]
#     return t <= cell.t[end]
# end

function step!(int::CellIntegrator, dt, model, env)
    int.sol.alive || return
    simulate_cell!(int, int.t + dt, model, env)
    savevalues!(int)
end

# function get_u(cell::Cell, t)
#     i = findfirst(isequal(t), cell.t)
#     @check !isnothing(i) "Cell state not saved at time t=$t"

#     cell.u[i]
# end


function interpolate(cell::Cell, t)
    @check cell.t[1] <= t <= cell.t[end]
    i = findlast(s -> s <= t, cell.t)

    cell.u[i]
end

function copy_until(cell::Cell, t)
    @check cell.t[1] <= t <= cell.t[end]

    idcs = findlast(s -> s <= t, cell.t)
    tt = [ cell.t[idcs]; t ]
    uu = [ cell.u[idcs]; interpolate(cell, t) ]

    Cell(cell.anc, deepcopy(tt), deepcopy(uu), deepcopy(cell.p), cell.alive)
end

function duplicate_cell(int::CellIntegrator, t::Float64)
    #@check cell_alive(cell, t)
    CellIntegrator(copy_until(int.sol, t))
end