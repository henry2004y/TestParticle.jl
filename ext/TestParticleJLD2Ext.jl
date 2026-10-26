module TestParticleJLD2Ext

using TestParticle
import TestParticle: save_trajectory
using SciMLBase: AbstractODESolution, EnsembleSolution
import JLD2
using JLD2: jldopen, Group

function _write_single_trajectory!(
        group, t_vec::AbstractVector, u_vec::AbstractVector;
        states_matrix::Bool = true
    )
    n = length(t_vec)
    dim = isempty(u_vec) ? 0 : length(u_vec[1])
    group["t"] = collect(t_vec)
    if dim >= 3
        group["x"] = [Float64(u[1]) for u in u_vec]
        group["y"] = [Float64(u[2]) for u in u_vec]
        group["z"] = [Float64(u[3]) for u in u_vec]
    end
    if dim >= 6
        group["vx"] = [Float64(u[4]) for u in u_vec]
        group["vy"] = [Float64(u[5]) for u in u_vec]
        group["vz"] = [Float64(u[6]) for u in u_vec]
    end
    if states_matrix && dim > 0
        mat = Matrix{Float64}(undef, dim, n)
        for j in 1:n
            for i in 1:dim
                mat[i, j] = u_vec[j][i]
            end
        end
        group["states"] = mat
    end
    return
end

function TestParticle.save_trajectory(
        filepath::AbstractString, sol::AbstractODESolution;
        states_matrix::Bool = true
    )
    return save_trajectory(filepath, sol.t, sol.u; states_matrix = states_matrix)
end

function TestParticle.save_trajectory(
        filepath::AbstractString, t_vec::AbstractVector, u_vec::AbstractVector;
        states_matrix::Bool = true
    )
    jldopen(filepath, "w") do f
        _write_single_trajectory!(f, t_vec, u_vec; states_matrix = states_matrix)
    end
    return filepath
end

function TestParticle.save_trajectory(
        filepath::AbstractString, ensemble::EnsembleSolution;
        states_matrix::Bool = true
    )
    jldopen(filepath, "w") do f
        f["num_trajectories"] = length(ensemble.u)
        for (idx, sol) in enumerate(ensemble.u)
            g = JLD2.Group(f, string("trajectory_", idx))
            _write_single_trajectory!(g, sol.t, sol.u; states_matrix = states_matrix)
        end
    end
    return filepath
end

function TestParticle.save_trajectory(
        filepath::AbstractString, raw::NamedTuple{(:u, :t)};
        states_matrix::Bool = true
    )
    jldopen(filepath, "w") do f
        f["t"] = collect(raw.t)
        if raw.u isa AbstractArray{<:Any, 3}
            f["num_trajectories"] = size(raw.u, 1)
            if states_matrix
                f["states"] = collect(raw.u)
            end
        else
            _write_single_trajectory!(f, raw.t, raw.u; states_matrix = states_matrix)
        end
    end
    return filepath
end

TestParticle.save_trajectory(sol::AbstractODESolution, filepath::AbstractString; kwargs...) =
    save_trajectory(filepath, sol; kwargs...)

TestParticle.save_trajectory(
    ensemble::EnsembleSolution, filepath::AbstractString; kwargs...
) = save_trajectory(filepath, ensemble; kwargs...)

TestParticle.save_trajectory(
    t_vec::AbstractVector, u_vec::AbstractVector, filepath::AbstractString; kwargs...
) = save_trajectory(filepath, t_vec, u_vec; kwargs...)

TestParticle.save_trajectory(
    raw::NamedTuple{(:u, :t)}, filepath::AbstractString; kwargs...
) = save_trajectory(filepath, raw; kwargs...)

end # module TestParticleJLD2Ext
