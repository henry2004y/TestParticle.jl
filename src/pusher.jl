# Native particle pusher

struct TraceProblem{TS, T<:Real, TP, PF}
   "initial condition"
   u0::TS
   "time span"
   tspan::Tuple{T, T}
   "time step"
   dt::T
   "(q2m, E, B)"
   p::TP
   "function for setting initial conditions"
   prob_func::PF
end

struct TraceSolution{TU<:Array, T<:AbstractVector}
   "positions and velocities"
   u::TU
   "time stamps"
   t::T
end

DEFAULT_PROB_FUNC(prob, i, repeat) = prob

function TraceProblem(u0, tspan, dt, p; prob_func=DEFAULT_PROB_FUNC)
   TraceProblem(u0, tspan, dt, p, prob_func)
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
      v⁻ = Vector{Float64}(undef, 3)
      v′ = Vector{Float64}(undef, 3)
      v⁺ = Vector{Float64}(undef, 3)
      t_rotate = Vector{Float64}(undef, 3)
      s_rotate = Vector{Float64}(undef, 3)
      v⁻_cross_t = Vector{Float64}(undef, 3)
      v′_cross_s = Vector{Float64}(undef, 3)

      new{typeof(v⁻)}(v⁻, v′, v⁺, t_rotate, s_rotate, v⁻_cross_t, v′_cross_s)
   end
end

@inline ODE_DEFAULT_ISOUTOFDOMAIN(u, p, t) = false

"""
    update_velocity!(xv, paramBoris, param, dt)

Updates velocity using the Boris method, Birdsall, Plasma Physics via Computer Simulation,
p.62. Reference: https://apps.dtic.mil/sti/citations/ADA023511
"""
function update_velocity!(xv, paramBoris, param,  dt)
	(; v⁻, v′, v⁺, t_rotate, s_rotate, v⁻_cross_t, v′_cross_s) = paramBoris
   q2m, E, B = param[1], param[2](xv, 0.0), param[3](xv, 0.0)
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

"Update location in one timestep `dt`."
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

"""
    solve(prob::TraceProblem; trajectories::Int=1,
       savestepinterval::Int=1, isoutofdomain::Function=ODE_DEFAULT_ISOUTOFDOMAIN)

Trace particles using the Boris method with specified `prob`.

# keywords
- `trajectories::Int`: number of trajectories to trace.
- `savestepinterval::Int`: saving output interval.
- `isoutofdomain::Function`: a function with input of position and velocity vector `xv` that determines whether to stop tracing.
"""
function solve(prob::TraceProblem, ensemblealg::BasicEnsembleAlgorithm=EnsembleSerial();
   trajectories::Int=1, savestepinterval::Int=1,
   isoutofdomain::Function=ODE_DEFAULT_ISOUTOFDOMAIN)

   sols = _solve(ensemblealg, prob, trajectories, savestepinterval, isoutofdomain)
end

function _solve(::EnsembleSerial, prob, trajectories, savestepinterval, isoutofdomain)
   sols, ttotal, nt, nout = _prepare(prob, trajectories, savestepinterval)

   _boris!(sols, prob, 1:trajectories, savestepinterval, ttotal, nt, nout, isoutofdomain)

   sols
end

function _solve(::EnsembleThreads, prob, trajectories, savestepinterval, isoutofdomain)
   sols, ttotal, nt, nout = _prepare(prob, trajectories, savestepinterval)

   nchunks = Threads.nthreads()
   Threads.@threads for (irange, ichunk) in chunks(1:trajectories, nchunks)
      _boris!(sols, prob, irange, savestepinterval, ttotal, nt, nout, isoutofdomain)
   end

   sols
end

"Prepare for advancing."
function _prepare(prob, trajectories, savestepinterval)
   (; tspan, dt) = prob
   ttotal = tspan[2] - tspan[1]
   nt = round(Int, ttotal / dt)
   nout = nt ÷ savestepinterval + 1
   sols = Vector{TraceSolution}(undef, trajectories)

   sols, ttotal, nt, nout
end

"Apply Boris method for particles with index in `irange`."
function _boris!(sols, prob, irange, savestepinterval, ttotal, nt, nout, isoutofdomain)
   (; tspan, dt, p, u0) = prob
   paramBoris = BorisMethod()
   xv = similar(u0)
   traj = zeros(eltype(u0), 6, nout)

   for i in irange
      # set initial conditions for each trajectory i
      iout = 1
      new_prob = prob.prob_func(prob, i, false)
      xv .= new_prob.u0
      traj[:,1] = xv

      # push velocity back in time by 1/2 dt
      update_velocity!(xv, paramBoris, p, -0.5*dt)

      for it in 1:nt
         update_velocity!(xv, paramBoris, p, dt)
         update_location!(xv, dt)
         if it % savestepinterval == 0
            iout += 1
            traj[:,iout] .= xv
         end
         isoutofdomain(xv, p, it*dt) && break
      end

      if iout == nout # regular termination
         dtfinal = ttotal - nt*dt
         if dtfinal > 1e-3 # final step if needed
            update_velocity!(xv, paramBoris, p, dtfinal)
            update_location!(xv, dtfinal)
            traj_save = hcat(traj, xv)
            t = [collect(tspan[1]:dt*savestepinterval:tspan[2])..., tspan[2]]
         else
            traj_save = copy(traj)
            t = collect(tspan[1]:dt*savestepinterval:tspan[2])
         end
      else # early termination or savestepinterval != 1
         traj_save = traj[:, 1:iout]
         t = collect(tspan[1]:dt*savestepinterval:tspan[1]+dt*savestepinterval*iout)
      end
      sols[i] = TraceSolution(traj_save, t)
   end

   return
end