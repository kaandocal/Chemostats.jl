if !@isdefined(_MODELS_RNGCELL_LOADED)
const _MODELS_RNGCELL_LOADED = true

include("exponential.jl")

using Distributions: Exponential

# Exponential lifetimes; dies with probability `p_die`, otherwise divides into two.
# All randomness comes from the cell's RNG.
function divide_rng(int)
    rng = Chemostats.get_cell_rng(int)
    rand(rng) < int.p.p_die && return nothing
    Chemostats.Offspring(CartesianIndex() => Exponential())
end

function RNGCell(seed; p_die = 0.2)
    rng = Xoshiro(seed)
    prob = ODEProblem(f_exp, randexp(rng), (0., 0.), (; p_die); callback = cb_exp)
    DECell(prob, Tsit5(), divide_rng; rng)
end

end
