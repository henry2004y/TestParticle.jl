using Documenter, DocumenterVitepress
using TestParticle

using Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "examples")
const OUTPUT_DIR = joinpath(@__DIR__, "src")

# Create output directories
mkpath(OUTPUT_DIR)
categories = ["drifts", "analytic", "applications", "features"]
for category in categories
   mkpath(joinpath(OUTPUT_DIR, category))
end

# Copy and clean index.md
index_src = joinpath(EXAMPLES_DIR, "overview.md")
index_out = joinpath(OUTPUT_DIR, "overview.md")
cp(index_src, index_out; force = true)

# Define orders manually
drifts_order = [
   "demo_ExB_drift.jl",
   "demo_gradient_B.jl",
   "demo_curvature_B.jl",
   "demo_gravity_drift.jl",
   "demo_polarization_drift.jl",
   "demo_FLR.jl"
]

analytic_order = [
   "demo_Buniform_Ezero.jl",
   "demo_electron_proton.jl",
   "demo_currentsheet.jl",
   "demo_magneticmirror.jl",
   "demo_magneticbottle.jl",
   "demo_adiabatic_periods.jl",
   "demo_proton_dipole.jl",
   "demo_analytic_magnetosphere.jl",
   "demo_tokamak_coil.jl",
   "demo_tokamak_profile.jl",
   "demo_spherical.jl",
   "demo_x_line.jl",
   "demo_zpinch.jl"
]

applications_order = [
   "demo_shock.jl",
   "demo_fermi_foreshock.jl",
   "demo_cosmicray.jl",
   "demo_radiation.jl",
   "demo_batsrus_3dstructured.md"
]

features_order = [
   "demo_boris.jl",
   "demo_energy_conservation.jl",
   "demo_phase_error.jl",
   "demo_performance.jl",
   "demo_interpolation.jl",
   "demo_dimensionless.jl",
   "demo_ensemble.jl",
   "demo_distributions.jl",
   "demo_savingcallback.jl",
   "demo_flux.jl",
   "demo_gc.jl",
   "demo_gpu.md",
   "demo_array.md",
   "demo_diskarray.jl",
   "demo_fieldlines.jl",
   "demo_current_loop.jl"
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
         push!(pages, joinpath(subdir, name))
      elseif endswith(filename, ".md")
         # Copy file
         cp(src_path, joinpath(out_dir, filename); force = true)
         push!(pages, joinpath(subdir, filename))
      end
   end
   return pages
end

drifts_pages = process_examples("drifts", drifts_order)
analytic_pages = process_examples("analytic", analytic_order)
applications_pages = process_examples("applications", applications_order)
features_pages = process_examples("features", features_order)

example_pages = [
   "Overview" => "overview.md",
   "Drifts" => drifts_pages,
   "Analytic Fields" => analytic_pages,
   "Applications" => applications_pages,
   "Features" => features_pages
]

makedocs(;
   modules = [TestParticle],
   authors = "Hongyang Zhou <hyzhou@umich.edu> and contributors",
   sitename = "TestParticle.jl",
   format = DocumenterVitepress.MarkdownVitepress(
      repo = "https://github.com/henry2004y/TestParticle.jl",
      devbranch = "master",
      devurl = "dev"
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

DocumenterVitepress.deploydocs(;
   repo = "github.com/henry2004y/TestParticle.jl",
   target = "build",
   devbranch = "master",
   branch = "gh-pages",
   push_preview = true
)
