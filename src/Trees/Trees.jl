"""
    Trees

Functionality for representing cell population trees and lineages.
"""
module Trees

using ArgCheck
using StructArrays
import StatsBase
using StatsBase: mean, var, mean_and_var
using QuadGK: quadgk

# Recording
export Tree
public cell_init!, cell_divide!, cell_die!, cell_update!

# Pruning
export mark!, prune!

# Queries
export roots, living, children, alive_at, ancestors, descendants
export has_divided, lifetime, is_alive

# Statistics
export time_average, time_weights, snapshot
export at_age, at_birth, at_division

include("tree.jl")
include("pruning.jl")
include("iterators.jl")
include("queries.jl")
include("collections.jl")
include("observables.jl")

end
