using Documenter
using DocumenterInterLinks
using Chemostats

links = InterLinks(
    "SciMLBase" => "https://docs.sciml.ai/SciMLBase/stable/",
    "ModelingToolkit" => "https://docs.sciml.ai/ModelingToolkit/stable/",
    "DifferentialEquations" => "https://docs.sciml.ai/DiffEqDocs/stable/",
    "Catalyst" => "https://docs.sciml.ai/Catalyst/stable/",
)

makedocs(
    sitename = "Chemostats.jl",
    modules = [ Chemostats ],
    format = Documenter.HTML(prettyurls = false),
    repo = Documenter.Remotes.GitHub("kaandocal", "Chemostats.jl"),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Usage" => "usage.md",
        "Theory" => "theory.md",
        "Algorithms" => "alg.md",
        "Using Chemostats.jl with DifferentialEquations.jl" => "decell.md",
    ],
    plugins = [ links ],
    linkcheck = true,
    checkdocs = :public,
)

deploydocs(
    repo = "github.com/kaandocal/Chemostats.jl.git",
    devbranch = "main",
    push_preview = true
)
