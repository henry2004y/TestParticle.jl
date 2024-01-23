# Native particle pusher

struct TraceProblem{TS, T<:Real, TP, PF}
   u0::TS
   tspan::Tuple{T, T}
   dt::T
   param::TP
   prob_func::PF
end

DEFAULT_PROB_FUNC(prob, i, repeat) = prob

function TraceProblem(u0, tspan, dt, param; prob_func=DEFAULT_PROB_FUNC)
   TraceProblem(u0, tspan, dt, param, prob_func)
end

struct BorisMethod{T, TV}
   "(q2m, E, B)"
   param::T
   # intermediate variables used in the solver
   v⁻::TV
   v′::TV
   v⁺::TV
   t_rotate::TV
   s_rotate::TV
   v⁻_cross_t::TV
   v′_cross_s::TV

   function BorisMethod(param)
      # intermediate variables
      v⁻ = Vector{Float64}(undef, 3)
      v′ = Vector{Float64}(undef, 3)
      v⁺ = Vector{Float64}(undef, 3)
      t_rotate = Vector{Float64}(undef, 3)
      s_rotate = Vector{Float64}(undef, 3)
      v⁻_cross_t = Vector{Float64}(undef, 3)
      v′_cross_s = Vector{Float64}(undef, 3)

      new{typeof(param), typeof(v⁻)}(param,
         v⁻, v′, v⁺, t_rotate, s_rotate, v⁻_cross_t, v′_cross_s)
   end
end

@inline ODE_DEFAULT_ISOUTOFDOMAIN(u) = false

"""
    update_velocity!(xv, paramBoris, dt)

Updates velocity using the Boris method, Birdsall, Plasma Physics via Computer Simulation,
p.62. Reference: https://apps.dtic.mil/sti/citations/ADA023511
"""
function update_velocity!(xv, paramBoris, dt)
	(; param, v⁻, v′, v⁺, t_rotate, s_rotate, v⁻_cross_t, v′_cross_s) = paramBoris
   q2m = param[1]
   E = param[2](xv, 0.0)
   B = param[3](xv, 0.0)
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
	   v⁻[dim] = xv[dim+3] + q2m*E[dim]*0.5*dt
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
	   xv[dim+3] = v⁺[dim] + q2m*E[dim]*0.5*dt
   end

   return
end

function update_location!(xv, dt)
   xv[1] += xv[4]*dt
   xv[2] += xv[5]*dt
   xv[3] += xv[6]*dt

   return
end

"In-place cross product."
function cross!(v1, v2, vout)
   vout[1] = v1[2]*v2[3] - v1[3]*v2[2]
   vout[2] = -v1[1]*v2[3] + v1[3]*v2[1]
   vout[3] = v1[1]*v2[2] - v1[2]*v2[1]

   return
end

function trace_trajectory(prob::TraceProblem; trajectories::Int=1, 
   savestepinterval::Int=1, isoutofdomain::Function=ODE_DEFAULT_ISOUTOFDOMAIN)

   trajs = Vector{Array{eltype(prob.u0),2}}(undef, trajectories)

   for i in 1:trajectories
      new_prob = prob.prob_func(prob, i, false)
      (; u0, tspan, dt, param) = new_prob
      xv = copy(u0)
      # prepare advancing
      ttotal = tspan[2] - tspan[1]
      nt = Int(ttotal ÷ dt)
      iout, nout = 1, nt ÷ savestepinterval + 1
      traj = zeros(eltype(u0), 6, nout)

      traj[:,1] = xv

      # push velocity back in time by 1/2 dt
      update_velocity!(xv, param, -0.5*dt)

      for it in 1:nt
         update_velocity!(xv, param, dt)
         update_location!(xv, dt)
         if it % savestepinterval == 0
         	iout += 1
         	traj[:,iout] .= xv
         end
         isoutofdomain(xv) && break
      end

      if iout == nout # regular termination
         # final step if needed
         dtfinal = ttotal - nt*dt
         if dtfinal > 1e-3
         	update_velocity!(xv, param, dtfinal)
         	update_location!(xv, dtfinal)
         	traj = hcat(traj, xv)
         end
      else # early termination
         traj = traj[:, 1:iout]
      end
      trajs[i] = traj
   end

   trajs
end
