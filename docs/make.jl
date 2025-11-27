using Documenter, DocumenterVitepress

using TestParticle

makedocs(;
   modules = [TestParticle],
   authors = "Hongyang Zhou <hyzhou@umich.edu> and contributors",
   sitename = "TestParticle.jl",
   format = DocumenterVitepress.MarkdownVitepress(
      repo = "https://github.com/henry2004y/TestParticle.jl",
      devurl = "dev",
      deploy_url = "henry2004y.github.io/TestParticle.jl",
   ),
   pages = [
      "Home" => "index.md",
      "Tutorial" => "tutorial.md",
      "Examples" => example_pages,
      "API" => "api.md",
      "Plot Functions" => "plotfunctions.md"
   ],
   warnonly = true,
)

deploydocs(;
    repo = "github.com/henry2004y/TestParticle.jl",
    push_preview=true,
)
