module Chemostats

using Random
using DataStructures
using ArgCheck
using EnumX

using SciMLBase
using SciMLBase: ReturnCode
using SciMLBase: EnsembleAlgorithm, EnsembleSerial, EnsembleThreads
import SciMLBase: savevalues!, step!, add_tstop!
import SymbolicIndexingInterface

using Trees: Trees, Tree, ancestors, descendants, children, alive_at, living, lifetime,
             time_average, time_weight, snapshot, at_age, at_birth, at_division
import Trees: is_alive

include("queue.jl")
include("worker.jl")
include("cell.jl")
include("tree.jl")
include("chemostat.jl")
include("algorithms.jl")

include("integrate.jl")
include("decell.jl")

export Chemostat, DECell, est_Λ, CellException
public CellState, Snapshot
public Forward, Thin, Direct, Strict, Lax

public simulate!
public get_curr_t, get_state, get_rng, get_data, est_logN, est_N

function SymbolicDivideContinuous end
function SymbolicDivideDiscrete end
function MTKDivideAffect end

end 