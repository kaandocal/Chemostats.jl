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
        "Theory" => "theory.md",
        "Algorithms" => "alg.md",
        "Using Chemostats.jl with DifferentialEquations.jl" => "decell.md",
    ],
    linkcheck = true,
    # Several internal functions (e.g. `iter_unsafe`, `prune!`) carry
    # docstrings for the sake of in-code clarity on tricky implementation
    # details, without being part of the public API -- :exports scopes the
    # "every docstring must appear in the manual" check to exported/public
    # bindings only, instead of flagging every internal docstring too.
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/kaandocal/Chemostats.jl.git",
    devbranch = "main",
    push_preview = true
)
