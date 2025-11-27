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
      repo = "github.com/henry2004y/TestParticle.jl",
      devbranch = "master",
      devurl = "dev",
      assets = assets,
   ),
   pages = [
      "Home" => "index.md",
      "Tutorial" => "tutorial.md",
      "Examples" => demos,
      "API" => "api.md",
      "Plot Functions" => "plotfunctions.md"
   ]
)

deploydocs(;
   repo = "github.com/henry2004y/TestParticle.jl",
   target = "build",
   branch = "gh-pages",
   devbranch = "master",
   push_preview = true
)
