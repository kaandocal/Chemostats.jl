# Generates the example pages in `src/examples` from the scripts in `demo/`, running the simulations.
#
# Usage: julia --project=docs docs/make_examples.jl [names...]   (all pages if no names are given)

using Literate

names = isempty(ARGS) ? [ "growth_protein", "threshold_division" ] : ARGS

for name in names
    Literate.markdown(joinpath(@__DIR__, "..", "demo", "$name.jl"), joinpath(@__DIR__, "src", "examples");
                      execute = true, documenter = false)
end
