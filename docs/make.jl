using TestParticle
using TestParticleMakie
using Documenter, DemoCards

branch = "master"
# generate demo files
demos, postprocess_cb, demo_assets = makedemos("examples"; branch)
# if there are generated css assets, pass it to Documenter.HTML
assets = String[]
isnothing(demo_assets) || (push!(assets, demo_assets))

makedocs(;
    modules=[TestParticle],
    authors="Hongyang Zhou <hyzhou@umich.edu> and contributors",
    repo="https://github.com/henry2004y/TestParticle.jl/blob/{commit}{path}#L{line}",
    sitename="TestParticle.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://henry2004y.github.io/TestParticle.jl",
        assets=String[],
        size_threshold=35000000,
        size_threshold_warn=10000000,
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Examples" => demos,
        "API" => "api.md",
        "Plot Functions" => "plotfunctions.md"
    ],
)

deploydocs(;
    repo="github.com/henry2004y/TestParticle.jl",
)
