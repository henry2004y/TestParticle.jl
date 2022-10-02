# Precompiling workloads

@precompile_all_calls begin
   # numerical field parameters
   x = range(-10, 10, length=15)
   y = range(-10, 10, length=20)
   z = range(-10, 10, length=25)
   B = fill(0.0, 3, length(x), length(y), length(z)) # [T]
   E = fill(0.0, 3, length(x), length(y), length(z)) # [V/m]

   B[3,:,:,:] .= 10e-9
   E[3,:,:,:] .= 5e-10

   Δx = x[2] - x[1]
   Δy = y[2] - y[1]
   Δz = z[2] - z[1]

   mesh = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
      (x[1], y[1], z[1]),
      (Δx, Δy, Δz))

   # numerical field
   param = prepare(mesh, E, B)
   # analytical field
   param = prepare(getE_dipole, getB_dipole)
end