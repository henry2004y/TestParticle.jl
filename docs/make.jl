using Documenter, DocumenterVitepress
using TestParticle
using SHA: sha256, bytes2hex

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
    "demo_curvature_gradient_B.jl",
    "demo_magnetic_drift.jl",
    "demo_drift_relation.jl",
    "demo_gravity_drift.jl",
    "demo_polarization_drift.jl",
    "demo_FLR.jl",
]

analytic_order = [
    "demo_Buniform_Ezero.jl",
    "demo_electron_proton.jl",
    "demo_currentsheet.jl",
    "demo_magneticmirror.jl",
    "demo_magneticbottle.jl",
    "demo_adiabatic_periods.jl",
    "demo_adiabatic_acceleration.jl",
    "demo_proton_dipole.jl",
    "demo_earth_gravity.jl",
    "demo_analytic_magnetosphere.jl",
    "demo_relativity_comparison.jl",
    "demo_tokamak_coil.jl",
    "demo_tokamak_profile.jl",
    "demo_spherical.jl",
    "demo_x_line.jl",
    "demo_zpinch.jl",
]

applications_order = [
    "demo_shock.jl",
    "demo_shock_drift.jl",
    "demo_shock_phase_space.jl",
    "demo_steady_phase_space.jl",
    "demo_fermi_foreshock.jl",
    "demo_cosmicray.jl",
    "demo_radiation.jl",
    "demo_gravitational_chaos.jl",
    "demo_batsrus_3dstructured.md",
]

features_order = [
    "demo_boris.jl",
    "demo_energy_conservation.jl",
    "demo_canonical.jl",
    "demo_accuracy.jl",
    "demo_adaptive_accuracy.jl",
    "demo_performance.jl",
    "demo_parallel_threads.md",
    "demo_parallel_distributed.md",
    "demo_interpolation.md",
    "demo_coevolution.jl",
    "demo_dimensionless.jl",
    "demo_unitful.md",
    "demo_ensemble.jl",
    "demo_distributions.jl",
    "demo_extrasaving.jl",
    "demo_flux.jl",
    "demo_gc.jl",
    "demo_hybrid.jl",
    "demo_gpu.md",
    "demo_array.md",
    "demo_diskarray.jl",
    "demo_fieldlines.jl",
    "demo_current_loop.jl",
    "demo_store_ensemble.jl",
]

const FORCE_REBUILD = get(ENV, "DOCS_FORCE_REBUILD", "false") == "true"

function process_examples(
        subdir::String, order::Vector{String};
        force_rebuild::Bool = FORCE_REBUILD
    )
    pages = String[]
    for filename in order
        src_path = joinpath(EXAMPLES_DIR, subdir, filename)
        out_dir = joinpath(OUTPUT_DIR, subdir)

        if endswith(filename, ".jl")
            name = replace(filename, ".jl" => ".md")
            md_path = joinpath(out_dir, name)
            hash_path = joinpath(out_dir, filename * ".sha256")
            current_hash = bytes2hex(sha256(read(src_path)))

            is_cached = !force_rebuild &&
                isfile(md_path) &&
                isfile(hash_path) &&
                read(hash_path, String) == current_hash

            if !is_cached
                @info "Generating $subdir/$filename"
                Literate.markdown(
                    src_path, out_dir;
                    flavor = Literate.DocumenterFlavor(),
                    execute = true,
                    credit = false
                )
                write(hash_path, current_hash)
            end
            push!(pages, joinpath(subdir, name))
        elseif endswith(filename, ".md")
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
    "Features" => features_pages,
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
        "Tutorial" => [
            "Introduction" => "tutorial/introduction.md",
            "Numerical Methods" => "tutorial/numerical_methods.md",
            "Boris Pusher" => "tutorial/boris.md",
            "Guiding Center" => "tutorial/guiding_center.md",
            "Features Walkthrough" => "tutorial/features.md",
        ],
        "Examples" => example_pages,
        "Python Wrapper" => "python_wrapper.md",
        "App" => "cli.md",
        "API" => "api.md",
        "Plot Functions" => "plotfunctions.md",
    ],
    warnonly = [:missing_docs, :linkcheck]
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/henry2004y/TestParticle.jl",
    target = "build",
    devbranch = "master",
    branch = "gh-pages",
    push_preview = true
)
