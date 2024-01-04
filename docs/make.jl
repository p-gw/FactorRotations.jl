using FactorRotations
using Documenter

DocMeta.setdocmeta!(
    FactorRotations,
    :DocTestSetup,
    :(using FactorRotations);
    recursive = true,
)

makedocs(;
    checkdocs = :exported,
    doctestfilters = [r"(\d*)\.(\d{4})\d+"],
    modules = [FactorRotations],
    authors = "Philipp Gewessler",
    repo = "https://github.com/p-gw/FactorRotations.jl/blob/{commit}{path}#{line}",
    sitename = "FactorRotations.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://p-gw.github.io/FactorRotations.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Guides" => ["Rotate an existing loading matrix" => "guides/basic_usage.md"],
        "Rotation Methods" => "rotation_methods.md",
        "API" => "api.md",
    ],
)

deploydocs(; repo = "github.com/p-gw/FactorRotations.jl", devbranch = "main")
