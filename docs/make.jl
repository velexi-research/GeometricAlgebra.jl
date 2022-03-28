# --- Setup

# Imports
using Documenter
using GeometricAlgebra

# Set up DocMeta
DocMeta.setdocmeta!(
    GeometricAlgebra, :DocTestSetup, :(using GeometricAlgebra); recursive=true
)

# --- Generate documentation

makedocs(;
    modules=[GeometricAlgebra],
    authors="Kevin Chu <kevin@velexi.com> and contributors",
    repo="https://github.com/velexi-corporation/GeometricAlgebra.jl/blob/{commit}{path}#{line}",
    sitename="GeometricAlgebra",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://velexi-corporation.github.io/GeometricAlgebra.jl/stable",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        #"Installation" => "installation.md",
        "API" => "api.md",
        #"Acknowledgements" => "acknowledgements.md",
        "Index" => "docs-index.md",
    ],
)

# --- Deploy documentation

deploydocs(; repo="github.com/velexi-corporation/GeometricAlgebra.jl", devbranch="main")
