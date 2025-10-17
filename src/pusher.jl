# Native particle pusher

struct TraceProblem{uType, tType, isinplace, P, F <: AbstractODEFunction, PF} <:
       AbstractODEProblem{uType, tType, isinplace}
   f::F
   "initial condition"
   u0::uType
   "time span"
   tspan::tType
   "(q2m, E, B)"
   p::P
   "function for setting initial conditions"
   prob_func::PF
end

get_EField(p::AbstractODEProblem) = get_EField(p.p)
get_BField(p::AbstractODEProblem) = get_BField(p.p)

struct TraceSolution{T, N, uType, uType2, DType, tType, rateType, P, A, IType, S, AC} <:
       AbstractODESolution{T, N, uType}
   "positions and velocities"
   u::uType
   u_analytic::uType2
   errors::DType
   "time stamps"
   t::tType
   k::rateType
   prob::P
   alg::A
   interp::IType
   dense::Bool
   tslocation::Int
   stats::S
   alg_choice::AC
   retcode::ReturnCode.T
end

function TraceSolution{T, N}(u, u_analytic, errors, t, k, prob, alg, interp, dense,
      tslocation, stats, alg_choice, retcode) where {T, N}
   return TraceSolution{T, N, typeof(u), typeof(u_analytic), typeof(errors), typeof(t),
      typeof(k), typeof(prob), typeof(alg), typeof(interp),
      typeof(stats),
      typeof(alg_choice)}(u, u_analytic, errors, t, k, prob, alg, interp,
      dense, tslocation, stats, alg_choice, retcode)
end

get_BField(sol::AbstractODESolution) = get_BField(sol.prob)
get_EField(sol::AbstractODESolution) = get_EField(sol.prob)

Base.length(ts::TraceSolution) = length(ts.t)

"""
Interpolate solution at time `x`. Forward tracing only.
"""
function (sol::TraceSolution)(
      t,
      ::Type{deriv} = Val{0};
      idxs = nothing,
      continuity = :left
) where {deriv}
   sol.interp(t, idxs, deriv, sol.prob.p, continuity)
end

DEFAULT_PROB_FUNC(prob, i, repeat) = prob

function TraceProblem(u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
   _f = ODEFunction{true, DEFAULT_SPECIALIZATION}(x -> nothing) # dummy func
   TraceProblem{typeof(u0), typeof(tspan), true, typeof(p), typeof(_f),
      typeof(prob_func)}(
      _f,
      u0,
      tspan,
      p,
      prob_func
   )
end
# For remake
function TraceProblem{iip}(; f, u0, tspan, p, prob_func) where {iip}
   TraceProblem{typeof(u0), typeof(tspan), iip, typeof(p), typeof(f),
      typeof(prob_func)}(
      f,
      u0,
      tspan,
      p,
      prob_func
   )
end

@inline ODE_DEFAULT_ISOUTOFDOMAIN(u, p, t) = false

"""
    update_velocity!(xv, param, dt, t)

Update velocity using the Boris method, Birdsall, Plasma Physics via Computer Simulation.
Reference: [DTIC](https://apps.dtic.mil/sti/citations/ADA023511)
"""
function update_velocity!(xv, p, dt, t)
   q2m, _, E, B = p
   n = size(xv, 2)
   E_field = Matrix{eltype(xv)}(undef, 3, n)
   B_field = Matrix{eltype(xv)}(undef, 3, n)
   # Get E and B fields for all particles
   for i in 1:n
      E_field[:, i] = E(xv[:, i], t)
      B_field[:, i] = B(xv[:, i], t)
   end

   @turbo for i in 1:n
      # half-step velocity update
      v_minus_x = xv[4, i] + q2m * E_field[1, i] * 0.5 * dt
      v_minus_y = xv[5, i] + q2m * E_field[2, i] * 0.5 * dt
      v_minus_z = xv[6, i] + q2m * E_field[3, i] * 0.5 * dt
      # rotation
      t_x = q2m * B_field[1, i] * 0.5 * dt
      t_y = q2m * B_field[2, i] * 0.5 * dt
      t_z = q2m * B_field[3, i] * 0.5 * dt

      t_mag2 = t_x * t_x + t_y * t_y + t_z * t_z

      s_x = 2 * t_x / (1 + t_mag2)
      s_y = 2 * t_y / (1 + t_mag2)
      s_z = 2 * t_z / (1 + t_mag2)
      # v' = v- + v- x t
      v_prime_x = v_minus_x + (v_minus_y * t_z - v_minus_z * t_y)
      v_prime_y = v_minus_y + (v_minus_z * t_x - v_minus_x * t_z)
      v_prime_z = v_minus_z + (v_minus_x * t_y - v_minus_y * t_x)
      # v+ = v- + v' x s
      v_plus_x = v_minus_x + (v_prime_y * s_z - v_prime_z * s_y)
      v_plus_y = v_minus_y + (v_prime_z * s_x - v_prime_x * s_z)
      v_plus_z = v_minus_z + (v_prime_x * s_y - v_prime_y * s_x)
      # half-step velocity update
      xv[4, i] = v_plus_x + q2m * E_field[1, i] * 0.5 * dt
      xv[5, i] = v_plus_y + q2m * E_field[2, i] * 0.5 * dt
      xv[6, i] = v_plus_z + q2m * E_field[3, i] * 0.5 * dt
   end
end


"""
Update location in one timestep `dt`.
"""
function update_location!(xv, dt)
   @turbo for i in 1:size(xv, 2)
      xv[1, i] += xv[4, i] * dt
      xv[2, i] += xv[5, i] * dt
      xv[3, i] += xv[6, i] * dt
   end
end

"""
     solve(prob::TraceProblem; trajectories::Int=1, dt::AbstractFloat,
    	 savestepinterval::Int=1, isoutofdomain::Function=ODE_DEFAULT_ISOUTOFDOMAIN)

Trace particles using the Boris method with specified `prob`.

# keywords

  - `trajectories::Int`: number of trajectories to trace.
  - `dt::AbstractFloat`: time step.
  - `savestepinterval::Int`: saving output interval.
  - `isoutofdomain::Function`: a function with input of position and velocity vector `xv` that determines whether to stop tracing.
"""
function solve(prob::TraceProblem, ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
      trajectories::Int = 1, savestepinterval::Int = 1, dt::AbstractFloat,
      isoutofdomain::Function = ODE_DEFAULT_ISOUTOFDOMAIN)
   sols = _solve(ensemblealg, prob, trajectories, dt, savestepinterval, isoutofdomain)
end

function _solve(::EnsembleSerial, prob, trajectories, dt, savestepinterval, isoutofdomain)
   sols, nt, nout = _prepare(prob, trajectories, dt, savestepinterval)
   irange = 1:trajectories
   _boris!(sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain)

   sols
end

function _solve(::EnsembleThreads, prob, trajectories, dt, savestepinterval, isoutofdomain)
   sols, nt, nout = _prepare(prob, trajectories, dt, savestepinterval)

   nchunks = Threads.nthreads()
   Threads.@threads for irange in index_chunks(1:trajectories; n = nchunks)
      _boris!(sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain)
   end

   sols
end

"""
Prepare for advancing.
"""
function _prepare(prob, trajectories, dt, savestepinterval)
   ttotal = prob.tspan[2] - prob.tspan[1]
   nt = round(Int, ttotal / dt) |> abs
   nout = nt ÷ savestepinterval + 1
   sols = Vector{TraceSolution}(undef, trajectories)

   sols, nt, nout
end

"""
Apply Boris method for particles with index in `irange`.
"""
function _boris!(sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain)
   if isoutofdomain !== ODE_DEFAULT_ISOUTOFDOMAIN
      @warn "The vectorized Boris pusher does not support `isoutofdomain`."
   end

   (; tspan, p, u0) = prob
   nparticles = length(irange)
   xv = Matrix{eltype(u0)}(undef, 6, nparticles)
   trajs = [Vector{SVector{6, eltype(u0)}}(undef, nout) for _ in 1:nparticles]

   # Set initial conditions
   for (j, i) in enumerate(irange)
      new_prob = prob.prob_func(prob, i, false)
      xv[:, j] .= new_prob.u0
      trajs[j][1] = SVector{6, eltype(u0)}(new_prob.u0)
   end

   # push velocity back in time by 1/2 dt
   update_velocity!(xv, p, -0.5*dt, -0.5*dt)

   iout = 1
   for it in 1:nt
      update_velocity!(xv, p, dt, (it-0.5)*dt)
      update_location!(xv, dt)
      if it % savestepinterval == 0
         iout += 1
         for j in 1:nparticles
            trajs[j][iout] = SVector{6, eltype(u0)}(view(xv, :, j))
         end
      end
   end

   t = range(tspan[1], step = dt*savestepinterval, length = nout) |> collect

   for (j, i) in enumerate(irange)
      traj_save = trajs[j]

      dense = false
      k = nothing
      alg = :boris
      alg_choice = nothing
      interp = LinearInterpolation(t, traj_save)
      retcode = ReturnCode.Default
      stats = nothing
      u_analytic = nothing
      errors = nothing
      tslocation = 0

      sols[i] = TraceSolution{Float64, 2}(traj_save, u_analytic, errors, t, k, prob, alg,
         interp, dense, tslocation, stats, alg_choice, retcode)
   end

   return
end
