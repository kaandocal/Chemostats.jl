using Documenter
using Chemostats

makedocs(
    sitename = "Chemostats.jl",
    modules = [ Chemostats ],
    format = Documenter.HTML(prettyurls = false),
    repo = "..",
    pages = [
        "Home" => "index.md",
        "Usage" => "usage.md",
        "Algorithms" => "alg.md",
        "Using Chemostats.jl with DifferentialEquations.jl" => "decell.md",
    ],
    linkcheck = true,
)

deploydocs(
    repo = "github.com/kaandocal/Chemostats.jl.git",
    devbranch = "main",
    push_preview = true
)
