# Generates the example pages in `src/examples` from the scripts in `demo/`, running the simulations.
#
# Usage: julia --project=docs docs/make_examples.jl

using Literate

for name in [ "growth_protein", "threshold_division" ]
    Literate.markdown(joinpath(@__DIR__, "..", "demo", "$name.jl"), joinpath(@__DIR__, "src", "examples");
                      execute = true, documenter = false)
end
