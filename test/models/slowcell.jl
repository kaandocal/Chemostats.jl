using Random
using Chemostats

# Custom cell type with random wall clock time to test race conditions
mutable struct SlowCell
    state::Chemostats.CellState.T
    t::Float64
    dt_next::Float64
end

SlowCell() = SlowCell(Chemostats.CellState.Newborn, 0.0, randexp())

Chemostats.get_curr_t(c::SlowCell) = c.t
Chemostats.get_state(c::SlowCell) = c.state
Chemostats.init_cell!(c::SlowCell) = (c.state = Chemostats.CellState.Alive)
Chemostats.die!(c::SlowCell) = (c.state = Chemostats.CellState.Dead)
Chemostats.divide!(c::SlowCell) = (c.state = Chemostats.CellState.Divided)
Chemostats.get_children(::SlowCell, p=nothing) = [ SlowCell(), SlowCell() ]

function Chemostats.step!(c::SlowCell, dt, p)
    sleep(0.001 + 0.01 * rand())
    rand() < 0.1 && error("SlowCell: random failure")

    if dt >= c.dt_next
        c.t += c.dt_next
        c.state = Chemostats.CellState.EndOfLife
    else
        c.t += dt
        c.dt_next -= dt
        c.state = Chemostats.CellState.Alive
    end
end
