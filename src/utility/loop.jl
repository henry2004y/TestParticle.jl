# Magnetic field of a current loop.

"""
    getB_loop(r, R, a, I, n)

Calculate the magnetic field `B` [T] at point `r` from a current loop with current `I` [A],
radius `a` [m], centered at `R`, and normal vector `n`.

# Arguments

  - `r::AbstractVector`: Position vector [x, y, z].
  - `R::AbstractVector`: Position of the loop center [x, y, z].
  - `a::Real`: Radius of the loop.
  - `I::Real`: Current in the loop.
  - `n::AbstractVector`: Normal vector of the loop (direction of the B-field at the center).
"""
function getB_loop(
      r::AbstractVector, R::AbstractVector, a::Real, I::Real, n::AbstractVector)
   # Normalize the normal vector
   n_hat = normalize(n)

   # Relative position from center
   r_rel = r - R

   # Project r_rel onto the normal vector to get z component in local coordinates
   z_local = dot(r_rel, n_hat)

   # Vector component perpendicular to n (rho vector)
   rho_vec = r_rel - z_local * n_hat
   rho = norm(rho_vec)

   # Handle the singularity on the wire itself (rho = a, z = 0)
   # and the center (rho = 0) separately if needed.
   # But standard formulas usually handle rho=0 if careful.

   if rho < 1e-10 * a # On the axis
      # B is purely in n direction
      # B = \mu_0 I a^2 / (2 (a^2 + z^2)^(3/2))
      B_mag = μ₀ * I * a^2 / (2 * (a^2 + z_local^2)^1.5)
      return B_mag * n_hat
   end

   # Cylindrical components calculation
   # B_rho and B_z in local coordinates
   # Formulas from Jackson or similar standard texts

   # k^2 = 4 a rho / ((a + rho)^2 + z^2)
   denom_sq = (a + rho)^2 + z_local^2
   k_sq = 4 * a * rho / denom_sq

   K_val = Elliptic.K(k_sq)
   E_val = Elliptic.E(k_sq)

   # Common factor
   factor = μ₀ * I / (2 * π * sqrt(denom_sq))

   # B_z (local)
   # B_z = factor * (K + (a^2 - rho^2 - z^2)/((a-rho)^2 + z^2) * E)
   denom_diff_sq = (a - rho)^2 + z_local^2
   B_z_local = factor * (K_val + (a^2 - rho^2 - z_local^2) / denom_diff_sq * E_val)

   # B_rho (local)
   # B_rho = factor * (z / rho) * (-K + (a^2 + rho^2 + z^2)/((a-rho)^2 + z^2) * E)
   if abs(z_local) < 1e-15
      B_rho_local = 0.0
   else
      B_rho_local = factor * (z_local / rho) *
                    (-K_val + (a^2 + rho^2 + z_local^2) / denom_diff_sq * E_val)
   end

   # Transform back to global Cartesian coordinates
   # B = B_rho * rho_hat + B_z * n_hat
   # rho_hat = rho_vec / rho

   rho_hat = rho_vec / rho

   B = B_rho_local * rho_hat + B_z_local * n_hat

   return B
end
