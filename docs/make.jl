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
    repo="https://github.com/velexi-research/GeometricAlgebra.jl/blob/{commit}{path}#{line}",
    sitename="GeometricAlgebra",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://velexi-research.github.io/GeometricAlgebra.jl/stable",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Types" => "types.md",
        "Operations" => "operations.md",
        "Implementation Details" => "implementation.md",
        "References" => "references.md",
        "Index" => "docs-index.md",
    ],
)

# --- Deploy documentation

deploydocs(; repo="github.com/velexi-research/GeometricAlgebra.jl", devbranch="main")
