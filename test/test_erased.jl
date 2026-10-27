module test_erased

using Test
using TestParticle
using StaticArrays
using BorisPushers
using KernelAbstractions
import TestParticle.ForwardDiff

const TP = TestParticle

@testset "ErasedField" begin
    @testset "Float64 and Float32 Dispatch" begin
        let
            B_func(r, t) = SA[r[1] + t, r[2], r[3]]
            field = TP.ErasedField(B_func)

            @test field isa TP.ErasedField{true}
            @test field isa TP.Field

            r64 = SA[1.0, 2.0, 3.0]
            val64 = field(r64, 0.5)
            @test val64 isa SVector{3, Float64}
            @test val64 ≈ SA[1.5, 2.0, 3.0]

            r32 = SA[1.0f0, 2.0f0, 3.0f0]
            val32 = field(r32, 0.5f0)
            @test val32 isa SVector{3, Float32}
            @test val32 ≈ SA[1.5f0, 2.0f0, 3.0f0]

            # Vector input conversion
            @test field([1.0, 2.0, 3.0], 0.5) ≈ SA[1.5, 2.0, 3.0]
            @test field(Float32[1.0, 2.0, 3.0], 0.5f0) ≈ SA[1.5f0, 2.0f0, 3.0f0]
        end
    end

    @testset "1-argument (Static) Fields" begin
        let
            B_static(r) = SA[r[2], r[1], 1.0]
            field = TP.ErasedField(B_static)

            @test field isa TP.ErasedField{false}

            r64 = SA[2.0, 3.0, 4.0]
            @test field(r64) ≈ SA[3.0, 2.0, 1.0]
            @test field(r64, 10.0) ≈ SA[3.0, 2.0, 1.0]

            r32 = SA[2.0f0, 3.0f0, 4.0f0]
            @test field(r32) ≈ SA[3.0f0, 2.0f0, 1.0f0]
            @test field(r32, 10.0f0) ≈ SA[3.0f0, 2.0f0, 1.0f0]
        end
    end

    @testset "ForwardDiff Autodiff Fallback" begin
        let
            B_func(r, t) = SA[r[1]^2, r[2], 0.0]
            field = TP.ErasedField(B_func)

            r0 = SA[2.0, 1.0, 0.0]
            jac = ForwardDiff.jacobian(r -> field(r, 0.0), r0)
            @test jac[1, 1] ≈ 4.0
            @test jac[2, 2] ≈ 1.0
        end
    end

    @testset "prepare integration and field_function unwrapping" begin
        let
            B_func(r) = [0.0, 0.0, 1.0]
            p = prepare(B_func)

            @test p.B isa TP.ErasedField{false}
            @test p.B.field_function === B_func

            # Idempotency
            field_erased = TP.ErasedField(p.B)
            @test field_erased === p.B
        end
    end

    @testset "GPU Adaptation" begin
        let
            B_func(r) = SA[0.0, 0.0, 1.0]
            field = TP.ErasedField(B_func)

            struct MockDevice end
            adapted = TP.Adapt.adapt(MockDevice(), field)
            @test !(adapted.field_function isa TP.ErasedFunction)
            @test adapted.field_function === B_func

            struct MockBackend <: KernelAbstractions.GPU end
            adapted_gpu = BorisPushers.adapt_field_to_gpu(field, MockBackend())
            @test !(adapted_gpu.field_function isa TP.ErasedFunction)
            @test adapted_gpu.field_function === B_func

            adapted_cpu = BorisPushers.adapt_field_to_gpu(field, BorisPushers.CPU())
            @test adapted_cpu === field
        end
    end

    @testset "Full-state (6-component) Evaluation & Type Stability" begin
        let
            # Velocity-dependent callback
            F_drag(xu) = SA[0.0, 0.0, -0.1 * xu[6]]
            p = prepare(; F = F_drag)
            y64 = SA[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
            val64 = @inferred p.F(y64, 0.0)
            @test val64 isa SVector{3, Float64}
            @test val64 ≈ SA[0.0, 0.0, -0.6]

            y32 = SA[1.0f0, 2.0f0, 3.0f0, 4.0f0, 5.0f0, 6.0f0]
            val32 = @inferred p.F(y32, 0.0f0)
            @test val32 isa SVector{3, Float32}
            @test val32 ≈ SA[0.0f0, 0.0f0, -0.6f0]

            # Vector input conversion for 6 components
            @test p.F([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 0.0) ≈ SA[0.0, 0.0, -0.6]

            # Strictly 3-component typed signature fallback
            B_typed(r::SVector{3, Float64}) = SA[r[1], r[2], r[3]]
            p_typed = prepare(; B = B_typed)
            val_typed = @inferred p_typed.B(y64, 0.0)
            @test val_typed ≈ SA[1.0, 2.0, 3.0]

            # ODE equations trace type stability
            B_func(r) = SA[0.0, 0.0, 1.0]
            p_trace = prepare(; B = B_func)
            trace_val = @inferred TestParticle.trace(y64, p_trace, 0.0)
            @test trace_val isa SVector{6, Float64}

            trace_rel_val = @inferred TestParticle.trace_relativistic(
                y64, p_trace, 0.0
            )
            @test trace_rel_val isa SVector{6, Float64}
        end
    end
end

end # module test_erased
