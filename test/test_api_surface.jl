module test_api_surface

using Test
import TestParticle as TP

module default_imports
    using TestParticle
end

@testset "Exported and public API" begin
    exported = (
        :prepare, :prepare_gc, :prepare_canonical, :prepare_erased,
        :PotentialField,
        :Proton, :Electron, :Ion,
        :TraceProblem, :TraceCanonicalProblem, :TraceGCProblem, :TraceHybridProblem,
        :solve, :remake,
        :trace, :trace!, :trace_relativistic, :trace_relativistic!,
        :trace_normalized, :trace_normalized!,
        :trace_relativistic_normalized, :trace_relativistic_normalized!,
        :trace_canonical, :trace_canonical!,
        :trace_canonical_relativistic, :trace_canonical_relativistic!,
        :trace_gc!, :trace_gc_drifts!, :trace_gc_flr!, :trace_gc_exb!,
        :AdaptiveHybrid, :GCRK4, :GCVern7, :GCRK45,
        :vdf_forward_trace, :vdf_backward_trace,
        :get_fields, :get_work_rates, :get_work, :orbit, :monitor,
        :save_trajectory, :build_interpolator,
        :EnsembleProblem, :EnsembleSerial, :EnsembleThreads, :EnsembleDistributed,
        :EnsembleSplitThreads, :DiscreteCallback, :TerminateOutside,
    )
    public_only = (
        :ReturnCode, :Field, :ErasedField, :ErasedFunction, :Species, :SpeciesDict,
        :ZeroField, :ZeroVector,
        :FullOrbit, :Canonical, :GuidingCenter,
        :CartesianGrid, :RectilinearGrid, :StructuredGrid,
        :qₑ, :mₑ, :qᵢ, :mᵢ, :c, :μ₀, :ϵ₀, :kB, :Rₑ, :BMoment_Earth, :eV,
        :get_gc, :get_gc_velocity,
        :full_to_gc, :gc_to_full, :velocity_to_canonical, :canonical_to_velocity,
        :canonical_hamiltonian,
        :get_gyrofrequency, :get_gyroperiod, :get_gyroradius, :get_energy,
        :energy2velocity, :get_curvature_radius, :get_adiabaticity,
        :adiabaticity_components,
        :get_particle_flux, :get_particle_crossings,
        :get_first_crossing,
        :Maxwellian, :BiMaxwellian, :Kappa, :BiKappa,
        :sample_unit_sphere, :sample_maxwellian, :generate_sphere,
        :jacobian, :derivative_t,
        :get_magnetic_properties, :get_work_channels,
        :LazyTimeInterpolator,
        :get_dx!,
        :get_dv!,
        :GPUGrid,
        :GPUUniformAxis,
        :GPUNonUniformAxis,
        :sph2cart,
        :cart2sph,
        :sph2cartvec,
        :cart2sphvec,
        :bin_centers,
        :bin_velocity_space,
        :project_vdf,
        :analytic_projection,
        :relative_l2,
        :velocity_moments,
        :vdf_grid_problem,
        :vdf_backward,
        :refine_vdf_window,
        :embed_vdf,
        :sample_velocity_ball,
        :PeriodicBC,
        :OnTheFly,
        :PreCompute,
        :StorePolicy,
        :FillExtrap, :ClampExtrap, :WrapExtrap,
    )

    @testset "Tracing workflow exports" begin
        for name in exported
            @test isdefined(TP, name)
            @test Base.isexported(TP, name)
            @test Base.ispublic(TP, name)
            @test isdefined(default_imports, name)
        end
        actual = Set(
            filter(
                name -> Base.isexported(TP, name) && name !== :TestParticle,
                names(TP),
            )
        )
        @test actual == Set(exported)
    end

    @testset "Qualified public APIs" begin
        for name in public_only
            @test isdefined(TP, name)
            @test Base.ispublic(TP, name)
            @test !Base.isexported(TP, name)
            @test !isdefined(default_imports, name)
        end
    end

    @test !Base.isexported(TP, :ZeroCurvBC)
end

end
