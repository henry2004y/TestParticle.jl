using LinearAlgebra: dot

struct MultistepBorisMethod{TV}
   t_n::TV
   e_n::TV
   v_cross_t::TV
   e_cross_t::TV
end

function MultistepBorisMethod(T::Type{<:AbstractFloat}=Float64)
   t_n = MVector{3, T}(undef)
   e_n = MVector{3, T}(undef)
   v_cross_t = MVector{3, T}(undef)
   e_cross_t = MVector{3, T}(undef)

   MultistepBorisMethod{typeof(t_n)}(t_n, e_n, v_cross_t, e_cross_t)
end

"""
    update_velocity_multistep!(xv, paramBoris, param, dt, t, n)

Update velocity using the Multistep Boris method.
Reference: [Zenitani & Kato 2025](https://arxiv.org/abs/2505.02270)
"""
function update_velocity_multistep!(xv, paramBoris, param, dt, t, n::Int)
   (; t_n, e_n, v_cross_t, e_cross_t) = paramBoris
   q2m, _, Efunc, Bfunc = param
   E = Efunc(xv, t)
   B = Bfunc(xv, t)

   # t_n and e_n vectors
   # Note: t_n in paper is (q/m * dt/(2n)) * B
   # e_n in paper is (q/m * dt/(2n)) * E
   # q2m is q/m
   factor = q2m * dt / (2 * n)

   @inbounds for dim in 1:3
      t_n[dim] = factor * B[dim]
      e_n[dim] = factor * E[dim]
   end

   t_n_mag2 = sum(abs2, t_n)
   t_n_mag = sqrt(t_n_mag2)

   # Calculate coefficients
   # c_n1 = cos(2n * alpha_n)
   # c_n2 = sin(2n * alpha_n) / t_n_mag
   # c_n3 = 2 * sin^2(n * alpha_n) / t_n_mag2
   # c_n6 = (2n - c_n2) / t_n_mag2

   # Check for small t_n to avoid division by zero or precision loss
   if t_n_mag < 1e-4
      # Taylor expansion limits as t_n -> 0
      # c_n1 ≈ 1 - 2*n^2*t_n^2
      # c_n2 ≈ 2n - 4/3*n^3*t_n^2
      # c_n3 ≈ 2*n^2
      # c_n6 ≈ 4/3*n^3

      c_n1 = 1.0 - 2*n^2*t_n_mag2
      c_n2 = 2*n - (4.0/3.0)*n^3*t_n_mag2
      c_n3 = 2.0*n^2
      c_n6 = (4.0/3.0)*n^3
   else
      alpha_n = atan(t_n_mag)
      n_alpha_n = n * alpha_n
      sin_n_alpha, cos_n_alpha = sincos(n_alpha_n)
      sin_2n_alpha = 2 * sin_n_alpha * cos_n_alpha
      cos_2n_alpha = cos_n_alpha^2 - sin_n_alpha^2

      c_n1 = cos_2n_alpha
      c_n2 = sin_2n_alpha / t_n_mag
      c_n3 = 2 * sin_n_alpha^2 / t_n_mag2
      c_n6 = (2*n - c_n2) / t_n_mag2
   end

   c_n4 = c_n2
   c_n5 = c_n3

   # Extract velocity
   v = @view xv[4:6]

   # Dot products
   v_dot_t = dot(v, t_n)
   e_dot_t = dot(e_n, t_n)

   # Cross products
   # v x t_n
   cross!(v, t_n, v_cross_t)

   # e_n x t_n
   cross!(e_n, t_n, e_cross_t)

   # Update velocity
   # Equation 39:
   # v_new = c_n1*v + c_n2*(v x t_n) + c_n3*(v . t_n)t_n + c_n4*e_n + c_n5*(e_n x t_n) + c_n6*(e_n . t_n)t_n

   @inbounds for i in 1:3
      xv[i+3] = c_n1 * xv[i+3] +
                c_n2 * v_cross_t[i] +
                c_n3 * v_dot_t * t_n[i] +
                c_n4 * e_n[i] +
                c_n5 * e_cross_t[i] +
                c_n6 * e_dot_t * t_n[i]
   end

   return
end

function _multistep_boris!(sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain, n_steps::Int)
   (; tspan, p, u0) = prob
   paramBoris = MultistepBorisMethod(eltype(u0))
   xv = MVector{6, eltype(u0)}(undef)
   # Safe initialization avoiding shared mutable elements
   traj = Vector{MVector{6, eltype(u0)}}(undef, nout)

   @fastmath @inbounds for i in irange
      # set initial conditions for each trajectory i
      iout = 1
      new_prob = prob.prob_func(prob, i, false)
      xv .= new_prob.u0
      traj[1] = copy(xv)

      # push velocity back in time by 1/2 dt using n=1 for initialization
      update_velocity_multistep!(xv, paramBoris, p, -0.5*dt, -0.5*dt, 1)

      for it in 1:nt
         update_velocity_multistep!(xv, paramBoris, p, dt, (it-0.5)*dt, n_steps)
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
      alg = :multistep_boris
      alg_choice = nothing
      interp = LinearInterpolation(t, traj_save)
      retcode = ReturnCode.Default
      stats = nothing
      u_analytic = nothing
      errors = nothing
      tslocation = 0

      sols[i] = TraceSolution{eltype(u0), 2}(traj_save, u_analytic, errors, t, k, prob, alg,
         interp, dense, tslocation, stats, alg_choice, retcode)
   end

   return
end
