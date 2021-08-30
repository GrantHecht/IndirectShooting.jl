using IndirectShooting
using Documenter

DocMeta.setdocmeta!(IndirectShooting, :DocTestSetup, :(using IndirectShooting); recursive=true)

makedocs(;
    modules=[IndirectShooting],
    authors="Grant Hecht",
    repo="https://github.com/GrantHecht/IndirectShooting.jl/blob/{commit}{path}#{line}",
    sitename="IndirectShooting.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://GrantHecht.github.io/IndirectShooting.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/GrantHecht/IndirectShooting.jl",
    devbranch="main",
)
