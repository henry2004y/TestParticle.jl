module TestParticleCLI

using LinearAlgebra: norm
using StaticArrays: SVector, SA
import TOML
using JLD2: jldopen
using TestParticle
import TestParticle as TP
using BorisPushers

struct SimulationConfig
    tspan::Tuple{Float64, Float64}
    dt::Float64
    saveat::Float64
    output_path::String
end

struct ParticleConfig
    stateinit::Vector{Float64}
    species_symbol::Symbol
    q::Float64
    m::Float64
end

abstract type AbstractFieldConfig end

struct UniformFieldConfig <: AbstractFieldConfig
    E::SVector{3, Float64}
    B::SVector{3, Float64}
end

struct DipoleFieldConfig <: AbstractFieldConfig
    M::Float64
    E::SVector{3, Float64}
end

struct HarrisSheetConfig <: AbstractFieldConfig
    B0::Float64
    L::Float64
    Bg::Float64
    E::SVector{3, Float64}
end

function parse_svector3(raw::Vector)::SVector{3, Float64}
    return SA[Float64(raw[1]), Float64(raw[2]), Float64(raw[3])]
end

function parse_field_config(fdict::Dict{String, Any})::AbstractFieldConfig
    ftype = lowercase(get(fdict, "type", "uniform"))
    E = haskey(fdict, "E") ? parse_svector3(fdict["E"]) : SA[0.0, 0.0, 0.0]

    if ftype == "uniform"
        B = haskey(fdict, "B") ? parse_svector3(fdict["B"]) : SA[0.0, 0.0, 0.01]
        return UniformFieldConfig(E, B)
    elseif ftype == "dipole"
        M = Float64(get(fdict, "M", norm(TP.BMoment_Earth)))
        return DipoleFieldConfig(M, E)
    elseif ftype == "harris"
        B0 = Float64(get(fdict, "B0", 0.01))
        L = Float64(get(fdict, "L", 1.0e5))
        Bg = Float64(get(fdict, "Bg", 0.0))
        return HarrisSheetConfig(B0, L, Bg, E)
    else
        println(Core.stdout, "Warning: unknown field type '", ftype, "', using uniform.")
        B = haskey(fdict, "B") ? parse_svector3(fdict["B"]) : SA[0.0, 0.0, 0.01]
        return UniformFieldConfig(E, B)
    end
end

function make_field_functions(cfg::UniformFieldConfig)
    E_field = (x, t) -> cfg.E
    B_field = (x, t) -> cfg.B
    return E_field, B_field
end

function make_field_functions(cfg::DipoleFieldConfig)
    E_field = (x, t) -> cfg.E
    # Earth-centered dipole: B = (μ0 M / 4π r^5) * [3 x z, 3 y z, 2 z^2 - x^2 - y^2]
    # In SI: μ0/(4π) * M = 1e-7 * M
    const_fac = 1.0e-7 * cfg.M
    B_field = (x, t) -> begin
        r2 = x[1]^2 + x[2]^2 + x[3]^2
        r = sqrt(r2)
        if r < 1.0e-3
            return SA[0.0, 0.0, 0.0]
        end
        fac = const_fac / (r2 * r2 * r)
        Bx = fac * (3.0 * x[1] * x[3])
        By = fac * (3.0 * x[2] * x[3])
        Bz = fac * (2.0 * x[3]^2 - x[1]^2 - x[2]^2)
        return SA[Bx, By, Bz]
    end
    return E_field, B_field
end

function make_field_functions(cfg::HarrisSheetConfig)
    E_field = (x, t) -> cfg.E
    # B = [B0 * tanh(z/L), Bg, 0]
    B_field = (x, t) -> SA[cfg.B0 * tanh(x[3] / cfg.L), cfg.Bg, 0.0]
    return E_field, B_field
end

function parse_particle_config(pdict::Dict{String, Any})::ParticleConfig
    sname = lowercase(get(pdict, "species", "electron"))
    pos = haskey(pdict, "position") ? parse_svector3(pdict["position"]) :
        SA[0.0, 0.0, 0.0]

    vel = if haskey(pdict, "velocity")
        parse_svector3(pdict["velocity"])
    elseif haskey(pdict, "energy_eV")
        energy_eV = Float64(pdict["energy_eV"])
        alpha_deg = Float64(get(pdict, "pitch_angle_deg", 90.0))
        alpha = deg2rad(alpha_deg)
        m = sname == "proton" ? TP.mᵢ : TP.mₑ
        v_mag = sqrt(2.0 * energy_eV * TP.eV / m)
        SA[v_mag * sin(alpha), 0.0, v_mag * cos(alpha)]
    else
        SA[0.0, 1.0e5, 0.0]
    end

    stateinit = [pos[1], pos[2], pos[3], vel[1], vel[2], vel[3]]

    if sname == "proton"
        return ParticleConfig(stateinit, :proton, TP.qᵢ, TP.mᵢ)
    elseif sname == "electron"
        return ParticleConfig(stateinit, :electron, TP.qₑ, TP.mₑ)
    else
        q = Float64(get(pdict, "q", TP.qₑ))
        m = Float64(get(pdict, "m", TP.mₑ))
        return ParticleConfig(stateinit, :custom, q, m)
    end
end

function parse_simulation_config(
        sdict::Dict{String, Any}, overrides::Dict{String, String}
    )::SimulationConfig
    tspan_raw = get(sdict, "tspan", [0.0, 1.0e-7])
    t0 = Float64(tspan_raw[1])
    t1 = Float64(tspan_raw[2])
    dt = Float64(get(sdict, "dt", 1.0e-10))
    saveat = Float64(get(sdict, "saveat", 10.0 * dt))
    output_path = String(get(sdict, "output", ""))

    if haskey(overrides, "dt")
        dt = parse(Float64, overrides["dt"])
    end
    if haskey(overrides, "saveat")
        saveat = parse(Float64, overrides["saveat"])
    end
    if haskey(overrides, "output")
        output_path = overrides["output"]
    end
    if haskey(overrides, "tspan1")
        t1 = parse(Float64, overrides["tspan1"])
    end

    return SimulationConfig((t0, t1), dt, saveat, output_path)
end

function write_trajectory_csv(
        io::IO, t_vec::AbstractVector{Float64}, u_vec::AbstractVector
    )
    println(io, "t,x,y,z,vx,vy,vz")
    for i in eachindex(t_vec)
        t = t_vec[i]
        u = u_vec[i]
        println(io, t, ",", u[1], ",", u[2], ",", u[3], ",", u[4], ",", u[5], ",", u[6])
    end
    return
end

function write_trajectory_jld2(
        filepath::String, t_vec::AbstractVector{Float64}, u_vec::AbstractVector
    )
    n = length(t_vec)
    x = zeros(Float64, n)
    y = zeros(Float64, n)
    z = zeros(Float64, n)
    vx = zeros(Float64, n)
    vy = zeros(Float64, n)
    vz = zeros(Float64, n)
    for i in 1:n
        u = u_vec[i]
        x[i] = u[1]
        y[i] = u[2]
        z[i] = u[3]
        vx[i] = u[4]
        vy[i] = u[5]
        vz[i] = u[6]
    end
    jldopen(filepath, "w") do f
        f["t"] = collect(t_vec)
        f["x"] = x
        f["y"] = y
        f["z"] = z
        f["vx"] = vx
        f["vy"] = vy
        f["vz"] = vz
        states = Matrix{Float64}(undef, 6, n)
        states[1, :] .= x
        states[2, :] .= y
        states[3, :] .= z
        states[4, :] .= vx
        states[5, :] .= vy
        states[6, :] .= vz
        f["states"] = states
    end
    return
end

function save_trajectory(
        filepath::String, t_vec::AbstractVector{Float64}, u_vec::AbstractVector
    )
    ext = lowercase(splitext(filepath)[2])
    if ext in (".h5", ".hdf5", ".jld2")
        write_trajectory_jld2(filepath, t_vec, u_vec)
    else
        open(filepath, "w") do io
            write_trajectory_csv(io, t_vec, u_vec)
        end
    end
    return
end

function print_usage()
    println(Core.stdout, "TestParticleCLI - Command Line Particle Tracer")
    println(Core.stdout, "Usage: tp_cli <config.toml> [options]")
    println(Core.stdout, "Options:")
    println(Core.stdout, "  --dt <val>       Override time step dt")
    println(Core.stdout, "  --saveat <val>   Override output save interval")
    println(Core.stdout, "  --output <path>  Write trajectory to file (.csv, .h5, .jld2)")
    println(Core.stdout, "  --tspan1 <val>   Override end time")
    println(Core.stdout, "  -h, --help       Show this help message")
    return
end

function @main(args::Vector{String})::Cint
    if isempty(args) || args[1] == "-h" || args[1] == "--help"
        print_usage()
        return 0
    end

    config_file = args[1]
    if !isfile(config_file)
        println(Core.stdout, "Error: configuration file not found: ", config_file)
        return 1
    end

    # Parse CLI overrides
    overrides = Dict{String, String}()
    i = 2
    while i <= length(args)
        arg = args[i]
        if arg == "--dt" && i + 1 <= length(args)
            overrides["dt"] = args[i + 1]
            i += 2
        elseif arg == "--saveat" && i + 1 <= length(args)
            overrides["saveat"] = args[i + 1]
            i += 2
        elseif arg == "--output" && i + 1 <= length(args)
            overrides["output"] = args[i + 1]
            i += 2
        elseif arg == "--tspan1" && i + 1 <= length(args)
            overrides["tspan1"] = args[i + 1]
            i += 2
        else
            i += 1
        end
    end

    toml_data = TOML.parsefile(config_file)
    sim_dict = get(toml_data, "simulation", Dict{String, Any}())
    part_dict = get(toml_data, "particle", Dict{String, Any}())
    field_dict = get(toml_data, "field", Dict{String, Any}())

    sim_cfg = parse_simulation_config(sim_dict, overrides)
    part_cfg = parse_particle_config(part_dict)
    field_cfg = parse_field_config(field_dict)

    E_fn, B_fn = make_field_functions(field_cfg)
    param = if part_cfg.species_symbol == :proton
        prepare(E_fn, B_fn; species = Proton)
    elseif part_cfg.species_symbol == :electron
        prepare(E_fn, B_fn; species = Electron)
    else
        prepare(E_fn, B_fn; q = part_cfg.q, m = part_cfg.m)
    end

    prob = TraceProblem(part_cfg.stateinit, sim_cfg.tspan, param)
    println(Core.stdout, "Tracing particle...")
    sol = solve(prob, Boris(); dt = sim_cfg.dt, saveat = sim_cfg.saveat)
    if sol.retcode != TP.ReturnCode.Success
        println(Core.stderr, "Error: simulation failed with return code ", sol.retcode)
        return 1
    end

    n_steps = length(sol.u)
    println(Core.stdout, "Simulation complete: ", n_steps, " steps.")
    u_end = sol.u[end]
    println(Core.stdout, "Final Position: [", u_end[1], ", ", u_end[2], ", ", u_end[3], "]")
    println(Core.stdout, "Final Velocity: [", u_end[4], ", ", u_end[5], ", ", u_end[6], "]")

    if !isempty(sim_cfg.output_path)
        save_trajectory(sim_cfg.output_path, sol.t, sol.u)
        println(Core.stdout, "Trajectory written to: ", sim_cfg.output_path)
    end

    return 0
end

end
