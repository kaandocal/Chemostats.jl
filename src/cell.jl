
""" 
    @enumx CellState 

The current state of a cell as returned by [`Chemostats.get_state`](@ref).

* `Newborn`: Newborn cell that has  not been simulated yet.
* `Alive`: Cell that is currently being simulated, but has not reached its end of life.
* `EndOfLive`: Cell that has reached its end of life, before its offspring are determined.
* `Dead`: Cell that has died, leaving no offspring.
* `Divided`: Cell that has divided into daughter cells.
* `Killed`: Cell that has been killed (e.g. removed), by the [Simulation Algorithm](@ref "Simulation Algorithms").
"""
@enumx CellState Newborn Alive EndOfLife Dead Divided Killed

"""
    divide!(cell)

Called after a cell divides, as determined by [`Chemostats.get_offspring`](@ref). 
After this call, [`Chemostats.get_state`](@ref) should return `Divided`.
"""
function divide! end 

"""
    die!(cell)

Called after a cell dies, as determined by [`Chemostats.get_offspring`](@ref). 
After this call, [`Chemostats.get_state`](@ref) should return `Dead`.
"""
function die! end 

"""
    kill!(cell, t)

Called when a cell is killed by the simulation algorithm. The argument `t` is the time at which a cell is killed.
After this call, [`Chemostats.get_state`](@ref) should return `Killed`.
"""
function kill! end 


"""
    get_curr_t(cell)

Returns the current simulation time `t` of the cell.
"""
function get_curr_t end 

"""
    get_state(cell)

Returns the current cell state, see [`CellState`](@ref).
"""
function get_state end 

"""
    get_offspring(chem, cell; save_lineages=false)

Returns an array of cells representing the offspring of the given cell. Here `chem` is the 
surrounding chemostat. If `save_lineages` is `true`, the resulting offspring should keep track
of the parent cell.
"""
function get_offspring end 

"""
    clone_cell(cell, t)

Creates an independent of the copy of the cell at time `t`. Used by [`Chemostats.Strict`](@ref) to maintain
a constant population size.
"""
function clone_cell end 
