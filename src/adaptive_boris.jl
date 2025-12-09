abstract type AbstractParticleSolver end

struct AdaptiveBoris{T <: AbstractFloat} <: AbstractParticleSolver
   dtmax::T
   dtmin::T
   safety::T
end

function AdaptiveBoris(; dtmax::AbstractFloat = Inf, dtmin::AbstractFloat = 0.0,
      safety::AbstractFloat = 0.9)
   AdaptiveBoris(dtmax, dtmin, safety)
end

function solve(prob::TraceProblem, alg::AdaptiveBoris,
      ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
      trajectories::Int = 1, savestepinterval::Int = 1,
      isoutofdomain::Function = ODE_DEFAULT_ISOUTOFDOMAIN,
      save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true)
   sols = _solve_adaptive(ensemblealg, prob, trajectories, savestepinterval,
      isoutofdomain, save_start, save_end, save_everystep, alg)
end

function _solve_adaptive(::EnsembleSerial, prob, trajectories, savestepinterval,
      isoutofdomain, save_start, save_end, save_everystep, alg)
   sols = Vector{TraceSolution}(undef, trajectories)
   irange = 1:trajectories
   _dispatch_boris_adaptive!(sols, prob, irange, savestepinterval, isoutofdomain,
      save_start, save_end, save_everystep, alg)
   sols
end

function _solve_adaptive(::EnsembleThreads, prob, trajectories, savestepinterval,
      isoutofdomain, save_start, save_end, save_everystep, alg)
   sols = Vector{TraceSolution}(undef, trajectories)
   nchunks = Threads.nthreads()
   Threads.@threads for irange in index_chunks(1:trajectories; n = nchunks)
      _dispatch_boris_adaptive!(sols, prob, irange, savestepinterval, isoutofdomain,
         save_start, save_end, save_everystep, alg)
   end
   sols
end

function _dispatch_boris_adaptive!(
      sols, prob, irange, savestepinterval, isoutofdomain,
      save_start, save_end, save_everystep, alg)
   is_td = is_time_dependent(get_EField(prob)) || is_time_dependent(get_BField(prob))
   _boris_adaptive!(sols, prob, irange, savestepinterval, isoutofdomain,
      save_start, save_end, save_everystep, Val(is_td), alg)
end

function _boris_adaptive!(sols, prob, irange, savestepinterval, isoutofdomain,
      save_start, save_end, save_everystep, ::Val{ITD}, alg) where {ITD}
   (; tspan, p, u0) = prob
   paramBoris = BorisMethod(eltype(u0))
   xv = MVector{6, eltype(u0)}(undef)
   v_old = MVector{3, eltype(u0)}(undef)

   # Use accessors if available, or robust unpacking
   # TestParticle standard p structure is (q2m, m, E, B, F)
   # update_velocity! unpacks: q2m, _, Efunc, Bfunc = param
   # So q2m is p[1], Bfunc is p[4] (or get_BField(prob))
   q2m = p[1]
   Bfunc = get_BField(prob)
   (; dtmax, dtmin, safety) = alg

   @fastmath @inbounds for i in irange
      new_prob = prob.prob_func(prob, i, false)
      xv .= new_prob.u0

      traj_dynamic = Vector{MVector{6, eltype(u0)}}()
      tsave_dynamic = Vector{typeof(tspan[1])}()

      if save_start
         push!(traj_dynamic, copy(xv))
         push!(tsave_dynamic, tspan[1])
      end

      # Calculate initial adaptive dt for the backward push
      B_local = Bfunc(xv, tspan[1])
      B_mag = norm(B_local)
      if B_mag > 0
         dt_ideal = safety / (abs(q2m) * B_mag)
      else
         dt_ideal = dtmax
      end
      dt_next = clamp(dt_ideal, dtmin, dtmax)

      # Initial backward half-step
      update_velocity!(xv, paramBoris, p, -0.5 * dt_next, tspan[1])

      t_current = tspan[1]
      it = 1

      while t_current < tspan[2]
         # Calculate adaptive dt for the NEXT step
         # We are at x_n. We want to move to x_{n+1}.
         # B at current position x_n
         B_local = Bfunc(xv, t_current)
         B_mag = norm(B_local)
         if B_mag > 0
            dt_ideal = safety / (abs(q2m) * B_mag)
         else
            dt_ideal = dtmax
         end

         dt_next = clamp(dt_ideal, dtmin, dtmax)

         # Don't overshoot
         if t_current + dt_next > tspan[2]
            dt_next = tspan[2] - t_current
         end

         v_old .= @view xv[4:6]
         t_field = ITD ? t_current + 0.5 * dt_next : zero(t_current)

         update_velocity!(xv, paramBoris, p, dt_next, t_field)

         if save_everystep && (it - 1) > 0 && (it - 1) % savestepinterval == 0
             xv_save = copy(xv) # x_n
             xv_save[4:6] .= v_old # v_{n-1/2}
             update_velocity!(xv_save, paramBoris, p, 0.5 * dt_next, ITD ? t_current : zero(t_current))

             push!(traj_dynamic, xv_save)
             push!(tsave_dynamic, t_current)
         end

         update_location!(xv, dt_next)

         if isoutofdomain(xv, p, t_current + dt_next)
             break
         end

         t_current += dt_next
         it += 1
      end

      # Save end
      if save_end
         xv_final = copy(xv)
         xv_final[4:6] .= v_old
         update_velocity!(xv_final, paramBoris, p, 0.5 * dt_next, ITD ? t_current : zero(t_current))
         push!(traj_dynamic, xv_final)
         push!(tsave_dynamic, t_current)
      end

      # Construct TraceSolution
      t = tsave_dynamic
      traj_save = traj_dynamic

      # Dummy fields
      dense = false
      k = nothing
      alg = :boris_adaptive
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
end
