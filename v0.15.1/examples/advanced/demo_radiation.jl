import DisplayAs #hide
using TestParticle: qᵢ as e, mₑ, c, ϵ₀
using CairoMakie
CairoMakie.activate!(type = "png") #hide

"""
    calc_gamma(β)

Calculate the Lorentz factor γ given β (v/c).
"""
function calc_gamma(β)
   if β >= 1.0
      return Inf
   end
   return 1.0 / sqrt(1.0 - β^2)
end

"""
    radiated_power(B, β, α)

Calculate the instantaneous power radiated by an electron in a magnetic field.
Based on the relativistic generalization of the Larmor formula.

Parameters:
B: Magnetic field strength [Tesla]
β: Normalized velocity (v/c)
α: Pitch angle [radians] (angle between v and B)

Returns:
P     : Radiated Power [Watts]
"""
function radiated_power(B, β, α)
   gamma = calc_gamma(β)

   # Formula: P = (e^4 * B^2 * γ^2 * v_perp^2) / (6 * π * ϵ₀ * m^2 * c^3)
   # Substitute v_perp = v * sin(α) = β * c * sin(α)
   # Result: P = (e^4 * B^2 * γ^2 * β^2 * sin(α)^2) / (6 * π * ϵ₀ * m^2 * c)

   numerator = e^4 * B^2 * gamma^2 * β^2 * sin(α)^2
   denominator = 6 * π * ϵ₀ * (mₑ^2) * c

   return numerator / denominator
end

function main()
   fig = Figure(size = (1200, 900), fontsize = 18)

   # --- Plot: Power vs Magnetic Field B ---
   # Fixed parameters
   beta_fixed_1 = 0.9
   alpha_fixed_1 = π / 2 # Perpendicular

   B_range = range(0.1, 10.0, length = 200) # 0.1 T to 10 T
   P_B = radiated_power.(B_range, beta_fixed_1, alpha_fixed_1)

   ax1 = Axis(fig[1, 1],
      title = L"Power vs. Magnetic Field ($B$)",
      xlabel = "Magnetic Field B [T]",
      ylabel = "Power [W]",
      xscale = log10,
      yscale = log10
   )
   lines!(ax1, B_range, P_B, color = :blue, linewidth = 3,
      label = L"\beta = 0.9, \alpha = 90^\circ")
   axislegend(ax1, position = :lt)
   text!(ax1, B_range[10], P_B[10], text = L"P \propto B^2", align = (:left, :bottom))

   # --- Plot: Power vs Velocity (β) ---
   # Fixed parameters
   B_fixed_2 = 1.0 # 1 Tesla
   alpha_fixed_2 = π / 2

   # Range close to 1 to show relativistic effects
   beta_range = range(0.1, 0.999, length = 500)
   P_beta = radiated_power.(B_fixed_2, beta_range, alpha_fixed_2)

   ax2 = Axis(fig[1, 2],
      title = L"Power vs. Velocity ($\beta$)",
      xlabel = L"Velocity $\beta = v/c$",
      ylabel = "Power [W]",
      yscale = log10
   )
   lines!(ax2, beta_range, P_beta, color = :red, linewidth = 3,
      label = L"B = 1\text{ T}, \alpha = 90^\circ")
   axislegend(ax2, position = :lt)

   # --- Plot: Power vs Pitch Angle (α) ---
   # Fixed parameters
   B_fixed_3 = 1.0
   beta_fixed_3 = 0.9

   alpha_range = range(0, π, length = 360)
   alpha_deg = rad2deg.(alpha_range)
   P_alpha = radiated_power.(B_fixed_3, beta_fixed_3, alpha_range)

   ax3 = Axis(fig[2, 1],
      title = L"Power vs. Pitch Angle ($\alpha$)",
      xlabel = "Pitch Angle [Degrees]",
      ylabel = "Power [W]",
      xticks = 0:30:180
   )
   lines!(ax3, alpha_deg, P_alpha, color = :green,
      linewidth = 3, label = L"B=1\text{ T}, \beta=0.9")
   # Add a visual marker for max power
   vlines!(ax3, [90], color = :black, linestyle = :dash)
   axislegend(ax3)

   # --- Plot: Heatmap (Log Power) vs B and β ---
   # Here we vary both B and β
   B_hm = range(0.1, 5.0, length = 100)
   beta_hm = range(0.5, 0.999, length = 100)

   # Calculate matrix
   P_matrix = log10.(radiated_power.(B_hm, beta_hm', π / 2))

   ax4 = Axis(fig[2, 2],
      title = "Log10(Power) Heatmap",
      xlabel = "Magnetic Field B [T]",
      ylabel = L"Velocity $\beta$"
   )
   hm = heatmap!(ax4, B_hm, beta_hm, P_matrix, colormap = :inferno)
   Colorbar(fig[2, 3], hm, label = "Log10(Power [W])")

   Label(fig[0, :], "Electron Cyclotron/Synchrotron Radiation Loss",
      fontsize = 24, font = :bold)

   fig
end

f = main()

f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
