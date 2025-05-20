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

struct BorisMethod{TV}
   # intermediate variables used in the solver
   v⁻::TV
   v′::TV
   v⁺::TV
   t_rotate::TV
   s_rotate::TV
   v⁻_cross_t::TV
   v′_cross_s::TV

   function BorisMethod()
      v⁻ = MVector{3, Float64}(undef)
      v′ = MVector{3, Float64}(undef)
      v⁺ = MVector{3, Float64}(undef)
      t_rotate = MVector{3, Float64}(undef)
      s_rotate = MVector{3, Float64}(undef)
      v⁻_cross_t = MVector{3, Float64}(undef)
      v′_cross_s = MVector{3, Float64}(undef)

      new{typeof(v⁻)}(v⁻, v′, v⁺, t_rotate, s_rotate, v⁻_cross_t, v′_cross_s)
   end
end

@inline ODE_DEFAULT_ISOUTOFDOMAIN(u, p, t) = false

"""
     update_velocity!(xv, paramBoris, param, dt, t)

Update velocity using the Boris method, Birdsall, Plasma Physics via Computer Simulation.
Reference: [DTIC](https://apps.dtic.mil/sti/citations/ADA023511)
"""
function update_velocity!(xv, paramBoris, param, dt, t)
   (; v⁻, v′, v⁺, t_rotate, s_rotate, v⁻_cross_t, v′_cross_s) = paramBoris
   q2m, _, Efunc, Bfunc = param
   E = Efunc(xv, t)
   B = Bfunc(xv, t)
   # t vector
   for dim in 1:3
      t_rotate[dim] = q2m*B[dim]*0.5*dt
   end
   t_mag2 = sum(abs2, t_rotate)
   # s vector
   for dim in 1:3
      s_rotate[dim] = 2*t_rotate[dim]/(1 + t_mag2)
   end
   # v-
   for dim in 1:3
      v⁻[dim] = xv[dim + 3] + q2m*E[dim]*0.5*dt
   end
   # v′
   cross!(v⁻, t_rotate, v⁻_cross_t)
   for dim in 1:3
      v′[dim] = v⁻[dim] + v⁻_cross_t[dim]
   end
   # v+
   cross!(v′, s_rotate, v′_cross_s)
   for dim in 1:3
      v⁺[dim] = v⁻[dim] + v′_cross_s[dim]
   end
   # v[n+1/2]
   for dim in 1:3
      xv[dim + 3] = v⁺[dim] + q2m*E[dim]*0.5*dt
   end

   return
end

"""
Update location in one timestep `dt`.
"""
function update_location!(xv, dt)
   xv[1] += xv[4]*dt
   xv[2] += xv[5]*dt
   xv[3] += xv[6]*dt

   return
end

"""
In-place cross product.
"""
function cross!(v1, v2, vout)
   vout[1] = v1[2]*v2[3] - v1[3]*v2[2]
   vout[2] = v1[3]*v2[1] - v1[1]*v2[3]
   vout[3] = v1[1]*v2[2] - v1[2]*v2[1]

   return
end

"""
     solve(prob::TraceProblem; trajectories::Int=1,
    	 savestepinterval::Int=1, isoutofdomain::Function=ODE_DEFAULT_ISOUTOFDOMAIN)

Trace particles using the Boris method with specified `prob`.

# keywords

  - `trajectories::Int`: number of trajectories to trace.
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
   (; tspan, p, u0) = prob
   paramBoris = BorisMethod()
   xv = MVector{6, eltype(u0)}(undef)
   traj = fill(MVector{6, eltype(u0)}(undef), nout)

   @fastmath @inbounds for i in irange
      # set initial conditions for each trajectory i
      iout = 1
      new_prob = prob.prob_func(prob, i, false)
      xv .= new_prob.u0
      traj[1] = copy(xv)

      # push velocity back in time by 1/2 dt
      update_velocity!(xv, paramBoris, p, -0.5*dt, -0.5*dt)

      for it in 1:nt
         update_velocity!(xv, paramBoris, p, dt, (it-0.5)*dt)
         update_location!(xv, dt)
         if it % savestepinterval == 0
            iout += 1
            traj[iout] = copy(xv)
         end
         isoutofdomain(xv, p, it*dt) && break
      end

      if iout == nout # regular termination
         traj_save = copy(traj)
         t = range(tspan[1], step = dt*savestepinterval, length = nout) |> collect
      else # early termination or savestepinterval != 1
         traj_save = traj[1:iout]
         t = tspan[1]:(dt * savestepinterval):(tspan[1] + dt * savestepinterval * (iout - 1)) |>
             collect
      end

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
