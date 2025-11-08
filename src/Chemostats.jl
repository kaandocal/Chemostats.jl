module Chemostats

using Random
using Distributions
using DataStructures
using ArgCheck
using UnPack
using EnumX

using SciMLBase
using SciMLBase: ReturnCode
using SciMLBase: EnsembleAlgorithm, EnsembleSerial, EnsembleThreads
import SciMLBase: savevalues!, step!, add_tstop!

include("cell.jl")
include("chemostat.jl")
include("algorithms.jl")
include("integrate.jl")

export Chemostat, Cell, simulate

end 