using TestParticle
using Documenter

makedocs(;
    modules=[TestParticle],
    authors="Hongyang Zhou <hyzhou@umich.edu> and contributors",
    repo="https://github.com/henry2004y/TestParticle.jl/blob/{commit}{path}#L{line}",
    sitename="TestParticle.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://henry2004y.github.io/TestParticle.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Example" => "example.md",
        "API" => "api.md",
        "Tutorial" => "tutorial.md"
    ],
)

deploydocs(;
    repo="github.com/henry2004y/TestParticle.jl",
)
