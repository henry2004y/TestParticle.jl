using StaticArrays

"The multistep Boris method."
abstract type AbstractBorisAlgorithm end

"A struct for the multistep Boris method."
struct MultistepBoris{K} <: AbstractBorisAlgorithm
   "Number of steps"
   k::Int
   "Coefficients for multistep E averaging"
   γ::NTuple{K, Float64}
   "History of electric field"
   E_history::Vector{MVector{3, Float64}}
   "current tracking index"
   head::Ref{Int}
   "number of saved E"
   len::Ref{Int}
   "Boris method for startup and for storing intermediate variables"
   boris::BorisMethod

   "Constructor for the multistep Boris method."
   function MultistepBoris(k::Int)
      if k == 2
         γ = (3/2, -1/2)
      elseif k == 3
         γ = (23/12, -16/12, 5/12)
      elseif k == 4
         γ = (55/24, -59/24, 37/24, -9/24)
      else
         throw(ArgumentError("Invalid k value for MultistepBoris! k should be 2, 3, or 4."))
      end
      E_history = [zeros(MVector{3, Float64}) for _ in 1:k]
      new{k}(k, γ, E_history, Ref(0), Ref(0), BorisMethod())
   end
end

"Update the history of the electric field."
function update_E_history!(paramBoris::MultistepBoris, E)
    k = length(paramBoris.E_history)
    paramBoris.head[] = mod1(paramBoris.head[] + 1, k)
    paramBoris.E_history[paramBoris.head[]] .= E
    if paramBoris.len[] < k
        paramBoris.len[] += 1
    end
end

"Get the averaged electric field."
function get_Ẽ(paramBoris::MultistepBoris)
    k = length(paramBoris.E_history)
    Ẽ = MVector{3, Float64}(0.0, 0.0, 0.0)
    for i in 1:k
        hist_idx = mod1(paramBoris.head[] - i + 1, k)
        Ẽ .+= paramBoris.γ[i] .* paramBoris.E_history[hist_idx]
    end
    return Ẽ
end

"Core of the multistep Boris method."
function _multistep_boris!(xv, boris, param, dt, Ẽ, B)
    (; v⁻, v′, v⁺, t_rotate, s_rotate, v⁻_cross_t, v′_cross_s) = boris
    q2m = param[1]

    # t vector
    @inbounds for dim in 1:3
       t_rotate[dim] = q2m*B[dim]*0.5*dt
    end
    t_mag2 = sum(abs2, t_rotate)
    # s vector
    @inbounds for dim in 1:3
       s_rotate[dim] = 2*t_rotate[dim]/(1 + t_mag2)
    end
    # v-
    @inbounds for dim in 1:3
       v⁻[dim] = xv[dim + 3] + q2m*Ẽ[dim]*0.5*dt
    end
    # v′
    cross!(v⁻, t_rotate, v⁻_cross_t)
    @inbounds for dim in 1:3
       v′[dim] = v⁻[dim] + v⁻_cross_t[dim]
    end
    # v+
    cross!(v′, s_rotate, v′_cross_s)
    @inbounds for dim in 1:3
       v⁺[dim] = v⁻[dim] + v′_cross_s[dim]
    end
    # v[n+1/2]
    @inbounds for dim in 1:3
       xv[dim + 3] = v⁺[dim] + q2m*Ẽ[dim]*0.5*dt
    end

    return
end

"Update velocity using the multistep Boris method."
function update_velocity_multistep!(xv, paramBoris::MultistepBoris, param, dt, t)
    q2m, _, Efunc, Bfunc = param
    E = Efunc(xv, t)
    B = Bfunc(xv, t)

    update_E_history!(paramBoris, E)

    if paramBoris.len[] < paramBoris.k
        # Use standard Boris for startup
        update_velocity!(xv, paramBoris.boris, param, dt, t)
    else
        # Use multistep Boris
        Ẽ = get_Ẽ(paramBoris)
        _multistep_boris!(xv, paramBoris.boris, param, dt, Ẽ, B)
    end
end
