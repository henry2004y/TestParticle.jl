# Precompiling workloads

@setup_workload begin
   @compile_workload begin
      # numerical field parameters
      x = range(-10, 10, length = 4)
      y = range(-10, 10, length = 6)
      z = range(-10, 10, length = 8)
      B = fill(0.0, 3, length(x), length(y), length(z)) # [T]
      E = fill(0.0, 3, length(x), length(y), length(z)) # [V/m]

      B[3, :, :, :] .= 10e-9
      E[3, :, :, :] .= 5e-10

      Δx = x[2] - x[1]
      Δy = y[2] - y[1]
      Δz = z[2] - z[1]

      mesh = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
         (x[1], y[1], z[1]),
         (Δx, Δy, Δz))

      vdf = Maxwellian([0.0, 0.0, 0.0], 1e-9, 1e6)
      v = sample(vdf)
      vdf = BiMaxwellian([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1e-9, 1e-9, 1e6)
      v = sample(vdf)
      # numerical field
      param = prepare(x, y, z, E, B)
      param = prepare(mesh, E, B)
      # analytical field
      param = prepare(getE_dipole, getB_dipole)
      # Boris pusher
      stateinit = let x0 = [0.0, 0.0, 0.0], v0 = [0.0, 1e5, 0.0]
         [x0..., v0...]
      end
      tspan = (0.0, 1.0)
      dt = 0.5
      prob = TraceProblem(stateinit, tspan, param)
      sol = solve(prob; dt, savestepinterval = 100)
      sol = solve(prob, EnsembleThreads(); dt, savestepinterval = 100)
      # guiding center
      X = get_gc(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0)
   end
end
