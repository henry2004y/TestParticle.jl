# Native pusher

struct TraceProblem{T<:Real, TP}
   stateinit::AbstractVector{T}
   tspan::Tuple{T, T}
   param::TP
end

"updates velocity using the Boris method, Birdsall, Plasma Physics via Computer Simulation, p.62"
function update_velocity!(xv, param, dt, v_minus, v_prime, v_plus, t, s, v_minus_cross_t,
   v_prime_cross_s)
   q, m = param[1], param[2]
   E = param[3](xv, 0.0)
   B = param[4](xv, 0.0)
	t_mag2 = 0.0
	# t vector
	for dim in 1:3
	   t[dim] = q/m*B[dim]*0.5*dt
	end
	# magnitude of t, squared
	t_mag2 = sum(abs2, t)
	# s vector
	for dim in 1:3
	   s[dim] = 2*t[dim]/(1 + t_mag2)
	end
	# v-
	for dim in 1:3
	   v_minus[dim] = xv[dim+3] + q/m*E[dim]*0.5*dt
   end
	# v′
   cross!(v_minus, t, v_minus_cross_t)
	for dim in 1:3
	   v_prime[dim] = v_minus[dim] + v_minus_cross_t[dim]
   end
	# v+
   cross!(v_prime, s, v_prime_cross_s)
	for dim in 1:3
	   v_plus[dim] = v_minus[dim] + v_prime_cross_s[dim]
   end
	# v[n+1/2]
	for dim in 1:3
	   xv[dim+3] = v_plus[dim] + q/m*E[dim]*0.5*dt
   end

   return
end

function update_location!(xv, dt)
	xv[1] += xv[4]*dt
	xv[2] += xv[5]*dt
	xv[3] += xv[6]*dt

   return
end

function cross!(v1, v2, vout)
	vout[1] = v1[2]*v2[3] - v1[3]*v2[2]
	vout[2] = -v1[1]*v2[3] + v1[3]*v2[1]
	vout[3] = v1[1]*v2[2] - v1[2]*v2[1]

   return
end

function trace_trajectory(param; dt, stateinit, tspan)
   xv = copy(stateinit)
   # intermediate variables
   v_minus = Vector{Float64}(undef, 3)
	v_prime = Vector{Float64}(undef, 3)
	v_plus = Vector{Float64}(undef, 3)
	t = Vector{Float64}(undef, 3)
	s = Vector{Float64}(undef, 3)
   v_minus_cross_t = Vector{Float64}(undef, 3)
   v_prime_cross_s = Vector{Float64}(undef, 3)
	# push velocity back in time by 1/2 dt
	update_velocity!(xv, param, -0.5*dt, v_minus, v_prime, v_plus, t, s, v_minus_cross_t,
      v_prime_cross_s)
   #TODO: get this right!
   nt = 1000
   output = zeros(6, nt)
   for it in 1:nt
      update_velocity!(xv, param, -0.5*dt, v_minus, v_prime, v_plus, t, s, v_minus_cross_t,
         v_prime_cross_s)
      update_location!(xv, dt)

      output[:,it] .= xv
   end

   output
end

