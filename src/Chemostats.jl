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
include("decell.jl")
include("chemostat.jl")
include("algorithms.jl")
include("integrate.jl")

export Chemostat
public CellState, Snapshot
public DECell, DivideCallback
public Forward, Thin, Direct, Strict, Lax

public simulate!
public get_curr_t, get_state, est_logN, est_N, est_Λ

end 