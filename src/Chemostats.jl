module Chemostats

using Random
using Distributions
using DataStructures
using ArgCheck
using UnPack

using SciMLBase
using SciMLBase: ReturnCode
import SciMLBase: savevalues!, step!, add_tstop!

include("cell.jl")
include("chemostat.jl")
include("algorithms.jl")
include("integrate.jl")

export Chemostat, Cell, simulate

end 