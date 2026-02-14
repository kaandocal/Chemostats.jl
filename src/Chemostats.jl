module Chemostats

using Random
using DataStructures
using ArgCheck
using UnPack
using EnumX

using SciMLBase
using SciMLBase: ReturnCode
using SciMLBase: EnsembleAlgorithm, EnsembleSerial, EnsembleThreads
import SciMLBase: savevalues!, step!, add_tstop!

include("queue.jl")
include("cell.jl")
include("chemostat.jl")
include("algorithms.jl")

include("integrate.jl")
include("decell.jl")

export Chemostat, DECell, est_Λ
public CellState, Snapshot
public Forward, Thin, Direct, Strict, Lax

public simulate!
public get_curr_t, get_state, est_logN, est_N

include("models/models.jl")

end 