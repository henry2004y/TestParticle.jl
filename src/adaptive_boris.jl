# Adaptive Boris method

struct AdaptiveBoris{T}
   dtmin::T
   dtmax::T
   safety::T
end

"""
    AdaptiveBoris(; dtmin, dtmax, safety=0.1)

Adaptive Boris method with adaptive time stepping based on local gyrofrequency.
The time step is determined by `dt = safety / |qB/m|`, clamped by `dtmin` and `dtmax`.
"""
function AdaptiveBoris(; dtmin, dtmax, safety = 0.1)
   T = promote_type(typeof(dtmin), typeof(dtmax), typeof(safety))
   AdaptiveBoris{T}(T(dtmin), T(dtmax), T(safety))
end

function solve(prob::TraceProblem, alg::AdaptiveBoris,
      ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
      trajectories::Int = 1, savestepinterval::Int = 1,
      isoutofdomain::Function = ODE_DEFAULT_ISOUTOFDOMAIN,
      save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true)
   _solve(ensemblealg, prob, trajectories, alg, savestepinterval, isoutofdomain,
      save_start, save_end, save_everystep)
end

function _solve(
      ::EnsembleSerial, prob, trajectories, alg::AdaptiveBoris, savestepinterval,
      isoutofdomain, save_start, save_end, save_everystep)
   # We cannot precalculate nt for adaptive steps
   sols = Vector{AbstractODESolution}(undef, trajectories)
   irange = 1:trajectories

   _adaptive_boris!(sols, prob, irange, alg, savestepinterval, isoutofdomain,
      save_start, save_end, save_everystep)

   sols
end

function _solve(
      ::EnsembleThreads, prob, trajectories, alg::AdaptiveBoris, savestepinterval,
      isoutofdomain, save_start, save_end, save_everystep)
   sols = Vector{AbstractODESolution}(undef, trajectories)

   nchunks = Threads.nthreads()
   Threads.@threads for irange in index_chunks(1:trajectories; n = nchunks)
      _adaptive_boris!(
         sols, prob, irange, alg, savestepinterval, isoutofdomain,
         save_start, save_end, save_everystep)
   end

   sols
end

function _adaptive_boris!(sols, prob, irange, alg, savestepinterval, isoutofdomain,
      save_start, save_end, save_everystep)
   (; tspan, p, u0) = prob
   q2m, _, Efunc, Bfunc = p
   T = eltype(u0)
   paramBoris = BorisMethod(T)
   xv = MVector{6, T}(undef)
   v_old = MVector{3, T}(undef)
   xv_save = MVector{6, T}(undef)

   # Determine if fields are time dependent
   is_td = is_time_dependent(get_EField(prob)) || is_time_dependent(get_BField(prob))

   @fastmath @inbounds for i in irange
      # Initialize solution containers
      initial_capacity = 1000
      traj = Vector{SVector{6, T}}(undef, 0)
      tsave = Vector{typeof(tspan[1])}(undef, 0)
      sizehint!(traj, initial_capacity)
      sizehint!(tsave, initial_capacity)

      iout = 0
      new_prob = prob.prob_func(prob, i, false)
      xv .= new_prob.u0
      t = tspan[1]

      if save_start
         iout += 1
         push!(traj, SVector{6, T}(xv))
         push!(tsave, t)
      end

      # Initial dt calculation
      B = Bfunc(xv, t)
      B_mag = norm(B)
      omega = abs(q2m) * B_mag
      dt = alg.safety / omega
      dt = clamp(dt, alg.dtmin, alg.dtmax)

      # Backstep velocity: v(0) -> v(-1/2) using dt
      update_velocity!(xv, paramBoris, p, -0.5 * dt, t)

      it = 1
      while t < tspan[2]
         # Check if next step exceeds tspan[2]
         if t + dt > tspan[2]
            dt_step = tspan[2] - t
            # Resync v from `t - 0.5*dt` to `t - 0.5*dt_step`
            t_sync = is_td ? t : zero(T)
            update_velocity!(xv, paramBoris, p, 0.5 * dt, t_sync)
            update_velocity!(xv, paramBoris, p, -0.5 * dt_step, t_sync)
            dt = dt_step
         end

         v_old .= @view xv[4:6] # v_{n-1/2} relative to current dt

         # Saving logic (start of step)
         if save_everystep && (it - 1) > 0 && (it - 1) % savestepinterval == 0
            iout += 1
            xv_save .= xv
            xv_save[4:6] .= v_old
            t_save = is_td ? t : zero(T)
            # Advance to t to get v_n
            update_velocity!(xv_save, paramBoris, p, 0.5 * dt, t_save)
            push!(traj, SVector{6, T}(xv_save))
            push!(tsave, t)
         end

         # Update velocity to v_{n+1/2}
         t_mid = is_td ? t + 0.5 * dt : zero(T)
         update_velocity!(xv, paramBoris, p, dt, t_mid)

         # Update location x_{n} -> x_{n+1}
         update_location!(xv, dt)
         t += dt

         if isoutofdomain(xv, p, t)
            break
         end
         it += 1

         # Calculate new dt for next step
         B = Bfunc(xv, t)
         B_mag = norm(B)
         omega = abs(q2m) * B_mag
         if omega > 0
            dt_new = alg.safety / omega
            dt_new = clamp(dt_new, alg.dtmin, alg.dtmax)
         else
            dt_new = alg.dtmax
         end

         # Resync for next step
         # v is at t_{new} - 0.5 * dt_old (relative to t_{new})
         # i.e. it is v_{n+1/2} from step we just took.
         t_sync = is_td ? t : zero(T)
         update_velocity!(xv, paramBoris, p, 0.5 * dt, t_sync)
         update_velocity!(xv, paramBoris, p, -0.5 * dt_new, t_sync)

         dt = dt_new
      end

      # Final save
      should_save_final = false
      if save_end
         should_save_final = true
      elseif save_everystep && (it - 1) > 0 && (it - 1) % savestepinterval == 0
         should_save_final = true
      end

      if should_save_final
         # xv currently has x at t (which is >= tspan[2] or boundary)
         # xv[4:6] has v at t - 0.5*dt (start of next step)
         # We want v at t.
         # So we just need to advance by 0.5 * dt
         t_final = is_td ? t : zero(T)
         update_velocity!(xv, paramBoris, p, 0.5 * dt, t_final)
         push!(traj, SVector{6, T}(xv))
         push!(tsave, t)
      end

      # Construct solution
      sol_alg = :adaptive_boris
      interp = LinearInterpolation(tsave, traj)
      retcode = ReturnCode.Default
      stats = nothing

      sols[i] = build_solution(
         prob, sol_alg, tsave, traj; interp = interp, retcode = retcode, stats = stats)
   end

   return
end
