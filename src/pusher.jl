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
end

function BorisMethod(T::Type{<:AbstractFloat} = Float64)
   v⁻ = MVector{3, T}(undef)
   v′ = MVector{3, T}(undef)
   v⁺ = MVector{3, T}(undef)
   t_rotate = MVector{3, T}(undef)
   s_rotate = MVector{3, T}(undef)
   v⁻_cross_t = MVector{3, T}(undef)
   v′_cross_s = MVector{3, T}(undef)

   BorisMethod{typeof(v⁻)}(v⁻, v′, v⁺, t_rotate, s_rotate, v⁻_cross_t, v′_cross_s)
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
      t_rotate[dim] = q2m * B[dim] * 0.5 * dt
   end
   t_mag2 = sum(abs2, t_rotate)
   # s vector
   for dim in 1:3
      s_rotate[dim] = 2 * t_rotate[dim] / (1 + t_mag2)
   end
   # v-
   for dim in 1:3
      v⁻[dim] = xv[dim + 3] + q2m * E[dim] * 0.5 * dt
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
      xv[dim + 3] = v⁺[dim] + q2m * E[dim] * 0.5 * dt
   end

   return
end

"""
Update location in one timestep `dt`.
"""
function update_location!(xv, dt)
   xv[1] += xv[4] * dt
   xv[2] += xv[5] * dt
   xv[3] += xv[6] * dt

   return
end

"""
In-place cross product.
"""
function cross!(v1, v2, vout)
   vout[1] = v1[2] * v2[3] - v1[3] * v2[2]
   vout[2] = v1[3] * v2[1] - v1[1] * v2[3]
   vout[3] = v1[1] * v2[2] - v1[2] * v2[1]

   return
end

"""
    solve(prob::TraceProblem; trajectories::Int=1, dt::AbstractFloat,
    savestepinterval::Int=1, isoutofdomain::Function=ODE_DEFAULT_ISOUTOFDOMAIN,
        n::Int=1)

Trace particles using the Boris method with specified `prob`.

# keywords

  - `trajectories::Int`: number of trajectories to trace.
  - `dt::AbstractFloat`: time step.
  - `savestepinterval::Int`: saving output interval.
  - `isoutofdomain::Function`: a function with input of position and velocity vector `xv` that determines whether to stop tracing.
  - `n::Int`: number of substeps for the Multistep Boris method. Default is 1 (standard Boris).
  - `save_start::Bool`: save the initial condition. Default is `true`.
  - `save_end::Bool`: save the final condition. Default is `true`.
  - `save_everystep::Bool`: save the state at every `savestepinterval`. Default is `true`.
"""
function solve(prob::TraceProblem, ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
      trajectories::Int = 1, savestepinterval::Int = 1, dt::AbstractFloat,
      isoutofdomain::Function = ODE_DEFAULT_ISOUTOFDOMAIN, n::Int = 1,
      save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true)
   sols = _solve(ensemblealg, prob, trajectories, dt, savestepinterval, isoutofdomain, n,
      save_start, save_end, save_everystep)
end

function _dispatch_boris!(
      sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain, n,
      save_start, save_end, save_everystep)

   is_td = is_time_dependent(get_EField(prob)) || is_time_dependent(get_BField(prob))
   if n == 1
      _boris!(sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain,
         save_start, save_end, save_everystep, Val(is_td))
   else
      _multistep_boris!(
         sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain, n,
         save_start, save_end, save_everystep, Val(is_td))
   end
end

function _solve(
      ::EnsembleSerial, prob, trajectories, dt, savestepinterval, isoutofdomain, n,
      save_start, save_end, save_everystep)
   sols, nt,
   nout = _prepare(prob, trajectories, dt, savestepinterval,
      save_start, save_end, save_everystep)
   irange = 1:trajectories
   _dispatch_boris!(sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain, n,
      save_start, save_end, save_everystep)

   sols
end

function _solve(
      ::EnsembleThreads, prob, trajectories, dt, savestepinterval, isoutofdomain, n,
      save_start, save_end, save_everystep)
   sols, nt,
   nout = _prepare(prob, trajectories, dt, savestepinterval,
      save_start, save_end, save_everystep)

   nchunks = Threads.nthreads()
   Threads.@threads for irange in index_chunks(1:trajectories; n = nchunks)
      _dispatch_boris!(
         sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain, n,
         save_start, save_end, save_everystep)
   end

   sols
end

"""
Prepare for advancing.
"""
function _prepare(prob::TraceProblem, trajectories, dt, savestepinterval,
      save_start, save_end, save_everystep)
   ttotal = prob.tspan[2] - prob.tspan[1]
   nt = round(Int, ttotal / dt) |> abs

   nout = 0
   if save_start
      nout += 1
   end

   if save_everystep
      steps = nt ÷ savestepinterval
      last_is_step = (nt > 0) && (nt % savestepinterval == 0)
      nout += steps
      if !save_end && last_is_step
         nout -= 1
      end
      if save_end && !last_is_step
         nout += 1
      end
   elseif save_end
      nout += 1
   end

   sols = Vector{TraceSolution}(undef, trajectories)

   sols, nt, nout
end

"""
Apply Boris method for particles with index in `irange`.
"""
function _boris!(sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain,
      save_start, save_end, save_everystep, ::Val{ITD}) where ITD
   (; tspan, p, u0) = prob
   paramBoris = BorisMethod(eltype(u0))
   xv = MVector{6, eltype(u0)}(undef)
   v_old = MVector{3, eltype(u0)}(undef)
   # Safe initialization avoiding shared mutable elements
   traj = Vector{MVector{6, eltype(u0)}}(undef, nout)
   tsave = Vector{typeof(tspan[1] + dt)}(undef, nout)

   @fastmath @inbounds for i in irange
      # set initial conditions for each trajectory i
      iout = 0
      new_prob = prob.prob_func(prob, i, false)
      xv .= new_prob.u0

      if save_start
         iout += 1
         traj[iout] = copy(xv)
         tsave[iout] = tspan[1]
      end

      # push velocity back in time by 1/2 dt
      update_velocity!(xv, paramBoris, p, -0.5 * dt, tspan[1])

      it = 1
      while it <= nt
         v_old .= @view xv[4:6]
         t = ITD ? (it - 0.5) * dt : zero(dt)
         update_velocity!(xv, paramBoris, p, dt, t)

         if save_everystep && (it - 1) > 0 && (it - 1) % savestepinterval == 0
            iout += 1
            if iout <= nout
               traj[iout] = copy(xv)
               traj[iout][4:6] .= v_old
               t_current = tspan[1] + (it - 1) * dt
               update_velocity!(traj[iout], paramBoris, p, 0.5 * dt,
                  ITD ? t_current : zero(dt))
               tsave[iout] = t_current
            end
         end

         update_location!(xv, dt)
         if isoutofdomain(xv, p, it * dt)
            break
         end
         it += 1
      end

      final_step = min(it, nt)
      should_save_final = false
      if save_end
         should_save_final = true
      elseif save_everystep && (final_step > 0) && (final_step % savestepinterval == 0)
         should_save_final = true
      end

      if iout < nout && should_save_final
         iout += 1
         t_final = final_step == nt ? tspan[2] : tspan[1] + final_step * dt
         dt_final = t_final - (tspan[1] + (final_step - 0.5) * dt)
         update_velocity!(xv, paramBoris, p, dt_final, ITD ? t_final : zero(dt))
         traj[iout] = copy(xv)
         tsave[iout] = t_final
      end

      if iout == nout # regular termination
         traj_save = copy(traj)
         t = tsave
      else # early termination or savestepinterval != 1
         traj_save = traj[1:iout]
         t = tsave[1:iout]
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

      sols[i] = TraceSolution{eltype(u0), 2}(
         traj_save, u_analytic, errors, t, k, prob, alg,
         interp, dense, tslocation, stats, alg_choice, retcode)
   end

   return
end
