module test_api_surface

using Test
import TestParticle as TP

@testset "Exported and public API" begin
    public_only = (
        :trace_canonical!,
        :trace_canonical,
        :trace_canonical_relativistic!,
        :trace_canonical_relativistic,
        :get_dx!,
        :get_dv!,
        :get_B_parameters,
        :get_E_parameters,
        :get_gc_velocity,
        :GPUGrid1D,
        :GPUGrid2D,
        :GPUGrid3D,
        :GPUSphericalGrid,
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
    )

    @test all(name -> isdefined(TP, name), public_only)
    @test all(name -> Base.ispublic(TP, name), public_only)
    @test all(name -> !Base.isexported(TP, name), public_only)
    @test !Base.isexported(TP, :ZeroCurvBC)
    @test Base.isexported(TP, :vdf_forward_trace)
    @test Base.isexported(TP, :vdf_backward_trace)
    @test all(
        name -> Base.isexported(TP, name),
        (:FillExtrap, :ClampExtrap, :WrapExtrap),
    )
end

end
