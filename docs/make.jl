using Documenter
using Chemostats

makedocs(
    sitename = "Chemostats.jl",
    modules = [ Chemostats ],
    format = Documenter.HTML(prettyurls = false),
    repo = Documenter.Remotes.GitHub("kaandocal", "Chemostats.jl"),
    pages = [
        "Home" => "index.md",
        "Usage" => "usage.md",
        "Theory" => "theory.md",
        "Algorithms" => "alg.md",
        "Using Chemostats.jl with DifferentialEquations.jl" => "decell.md",
    ],
    linkcheck = true,
    checkdocs = :public,
)

deploydocs(
    repo = "github.com/kaandocal/Chemostats.jl.git",
    devbranch = "main",
    push_preview = true
)
