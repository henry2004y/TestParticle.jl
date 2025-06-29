using TestParticle
using Documenter, DemoCards
using DocumenterVitepress

branch = "master"
# generate demo files
demos, postprocess_cb, demo_assets = makedemos("examples"; branch)
# if there are generated css assets, pass it to Documenter.HTML
assets = String[]
isnothing(demo_assets) || (push!(assets, demo_assets))

makedocs(;
   modules = [TestParticle],
   authors = "Hongyang Zhou <hyzhou@umich.edu> and contributors",
   sitename = "TestParticle.jl",
   format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://henry2004y.github.io/TestParticle.jl",
    ),
   pages = [
      "Home" => "index.md",
      "Tutorial" => "tutorial.md",
      "Examples" => demos,
      "API" => "api.md",
      "Plot Functions" => "plotfunctions.md"
   ]
)

DocumenterVitepress.deploydocs(;
   repo = "github.com/henry2004y/TestParticle.jl",
   target = "build", # this is where Vitepress stores its output
   devbranch = "master",
   branch = "gh-pages",
   push_preview = true,
)
