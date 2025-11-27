using Documenter, DocumenterVitepress
using TestParticle

using Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "examples")
const OUTPUT_DIR = joinpath(@__DIR__, "generated")

# Create output directories
mkpath(OUTPUT_DIR)
mkpath(joinpath(OUTPUT_DIR, "basics"))
mkpath(joinpath(OUTPUT_DIR, "advanced"))

# Copy and clean index.md
index_src = joinpath(EXAMPLES_DIR, "index.md")
index_out = joinpath(OUTPUT_DIR, "index.md")
content = read(index_src, String)
content = replace(content, "{{{democards}}}" => "")
write(index_out, content)

# Define orders manually
basics_order = [
   "demo_energy_conservation.jl",
   "demo_boris.jl",
   "demo_Buniform_Ezero.jl",
   "demo_dimensionless.jl",
   "demo_dimensionless_periodic.jl",
   "demo_dimensionless_dimensional.jl",
   "demo_electron_proton.jl",
   "demo_multiple.jl",
   "demo_ExB_drift.jl",
   "demo_gravity_drift.jl",
   "demo_gradient_B.jl",
   "demo_curvature_B.jl",
   "demo_FLR.jl",
   "demo_polarization_drift.jl",
   "demo_array.md"
]

advanced_order = [
   "demo_boris_outofdomain.jl",
   "demo_cosmicray.jl",
   "demo_ensemble.jl",
   "demo_flux.jl",
   "demo_savingcallback.jl",
   "demo_output_func.jl",
   "demo_currentsheet.jl",
   "demo_magneticmirror.jl",
   "demo_magneticbottle.jl",
   "demo_proton_dipole.jl",
   "demo_analytic_magnetosphere.jl",
   "demo_shock.jl",
   "demo_fermi_foreshock.jl",
   "demo_spherical.jl",
   "demo_tokamak_coil.jl",
   "demo_tokamak_profile.jl",
   "demo_gc.jl",
   "demo_batsrus_3dstructured.md",
   "demo_radiation.jl",
   "demo_gpu.md"
]

function process_examples(subdir, order)
   pages = String[]
   for filename in order
      src_path = joinpath(EXAMPLES_DIR, subdir, filename)
      out_dir = joinpath(OUTPUT_DIR, subdir)

      if endswith(filename, ".jl")
         # Compile to markdown
         Literate.markdown(src_path, out_dir; documenter = true, credit = false)
         name = replace(filename, ".jl" => ".md")
         push!(pages, joinpath("generated", subdir, name))
      elseif endswith(filename, ".md")
         # Copy file
         cp(src_path, joinpath(out_dir, filename); force = true)
         push!(pages, joinpath("generated", subdir, filename))
      end
   end
   return pages
end

basics_pages = process_examples("basics", basics_order)
advanced_pages = process_examples("advanced", advanced_order)

example_pages = [
   "Overview" => "generated/index.md",
   "Basics" => basics_pages,
   "Advanced" => advanced_pages
]

makedocs(;
   modules = [TestParticle],
   authors = "Hongyang Zhou <hyzhou@umich.edu> and contributors",
   sitename = "TestParticle.jl",
   format = DocumenterVitepress.MarkdownVitepress(
      repo = "https://github.com/henry2004y/TestParticle.jl",
      devurl = "dev",
      deploy_url = "henry2004y.github.io/TestParticle.jl"
   ),
   pages = [
      "Home" => "index.md",
      "Tutorial" => "tutorial.md",
      "Examples" => example_pages,
      "API" => "api.md",
      "Plot Functions" => "plotfunctions.md"
   ],
   warnonly = true
)

deploydocs(;
   repo = "github.com/henry2004y/TestParticle.jl",
   push_preview = true
)
