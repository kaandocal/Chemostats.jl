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

include("Trees/Trees.jl")

using .Trees: Trees, Tree, ancestors, descendants, children, roots, has_divided, alive_at, living, lifetime,
             time_average, time_weights, snapshot, at_age, at_birth, at_division
import .Trees: is_alive

include("queue.jl")
include("worker.jl")
include("cell.jl")
include("recording.jl")
include("chemostat.jl")
include("algorithms.jl")

include("integrate.jl")
include("offspring.jl")
include("decell.jl")

export Chemostat, DECell, est_Λ, CellException
public CellState, Snapshot
public Forward, Thin, Direct, Strict, Lax

public simulate, simulate!, DivideCallback
public Offspring, ExactPartition, BinomialPartition, AllOrNothingPartition
public get_curr_t, get_state, get_cell_rng, get_data, est_logN, est_N


end 