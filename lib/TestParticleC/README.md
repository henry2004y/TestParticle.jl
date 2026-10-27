# TestParticleC

`TestParticleC` is a lightweight C-ABI wrapper around `BorisPushers.jl` providing
statically callable C symbols for high-performance particle pushing in any language
with a C-compatible ABI or FFI (such as C, C++, Fortran, Rust, Python, or Zig)
in MHD and PIC simulation codes.

## Features & C API (`include/testparticle.h`)

### 1. Standard Boris Pushers
- **Synchronized Push**:
  - `testparticle_boris_push_f64(x, v, q2m, dt, E, B)`
  - `testparticle_boris_push_f32(x, v, q2m, dt, E, B)`
- **Staggered Leapfrog Push**:
  - `testparticle_boris_step_staggered_f64(x, v_half, q2m, dt, E, B)`
  - `testparticle_boris_step_staggered_f32(x, v_half, q2m, dt, E, B)`
- **Pure Velocity Rotation**:
  - `testparticle_update_velocity_boris_f64(v, E, B, qdt2m)`
  - `testparticle_update_velocity_boris_f32(v, E, B, qdt2m)`

### 2. High-Order & Multistep Hyper-Boris Pushers (Zenitani & Kato, 2025)
- **General Multistep Push** (`order` $\in \{2, 4, 6\}$, `n` subcycles):
  - `testparticle_multistep_boris_push_f64(x, v, q2m, dt, E, B, n, order)`
  - `testparticle_multistep_boris_step_staggered_f64(x, v_half, q2m, dt, E, B, n, order)`
  - Corresponding `_f32` variants.
- **Specialized 4th & 6th Order Hyper-Boris** (eliminates branch checks):
  - `testparticle_hyper_boris4_push_f64(x, v, q2m, dt, E, B, n)`
  - `testparticle_hyper_boris6_push_f64(x, v, q2m, dt, E, B, n)`
  - Corresponding `_f32` variants.

### 3. Adaptive Time Stepping & Velocity Resynchronization
- **Gyroperiod-tracking step size proposal**:
  - `testparticle_boris_propose_dt_f64(q2m, B, safety)`
- **Half-step velocity resynchronization** (preserves time-reversibility across $\Delta t$ changes):
  - `testparticle_boris_resync_velocity_half_f64(v_half, q2m, dt_prev, dt_next, E, B)`
- **Integrated adaptive step**:
  - `testparticle_adaptive_boris_step_f64(x, v_half, dt_current, &dt_next, q2m, E, B, safety)`

### 4. Vectorized Batch Pushers (AoS, SoA, and Multi-Species)
- **Array of Structures (AoS)**:
  - `testparticle_boris_push_batch_f64(n, x, v, q2m, dt, E, B)`
  - `testparticle_multistep_boris_push_batch_f64(n, x, v, q2m, dt, E, B, n_sub, order)`
- **Structure of Arrays (SoA)** (for SIMD / vector PIC codes):
  - `testparticle_boris_push_batch_soa_f64(n, x, y, z, vx, vy, vz, q2m_arr, q2m_scalar, dt, Ex, Ey, Ez, Bx, By, Bz)`
  - `testparticle_boris_step_staggered_batch_soa_f64(...)`
- **Multi-species plasma arrays**:
  - `testparticle_boris_push_batch_multispecies_f64(n, x, v, q2m_arr, dt, E, B)`

### 5. Trajectory Tracing with Dynamic Field Callbacks
- Full-trajectory orbit integration with C function pointer field evaluations:
  - `testparticle_trace_boris_f64(x, v, t_start, t_end, dt, q2m, field_fn, user_data, max_steps, &n_taken, traj_out)`

---

## Compiling to Shared Library with JuliaC

To compile into a standalone shared library (`libtestparticle.so` on Linux,
`libtestparticle.dylib` on macOS, or `libtestparticle.dll` on Windows):

```bash
juliac --output-lib libtestparticle.so \
       --project=lib/TestParticleC \
       lib/TestParticleC/src/TestParticleC.jl
```

---

## Example C Usage

```c
#include "testparticle.h"
#include <stdio.h>

int main() {
    double x[3] = {0.0, 0.0, 0.0};
    double v[3] = {1.0e5, 0.0, 0.0};
    double E[3] = {0.0, 0.0, 0.0};
    double B[3] = {0.0, 0.0, 1.0};
    double q2m = 1.0e8;
    double dt = 1.0e-9;

    // Standard Boris push
    testparticle_boris_push_f64(x, v, q2m, dt, E, B);

    // Hyper-Boris 4th-order push with 2 subcycles
    testparticle_hyper_boris4_push_f64(x, v, q2m, dt, E, B, 2);

    printf("New position: (%e, %e, %e)\n", x[0], x[1], x[2]);
    printf("New velocity: (%e, %e, %e)\n", v[0], v[1], v[2]);
    return 0;
}
```

Link against `-ltestparticle` when compiling:
```bash
gcc -O3 main.c -Iinclude -L. -ltestparticle -o main
```

---

## Example Fortran Usage (`ISO_C_BINDING`)

```fortran
program test_boris
    use, intrinsic :: iso_c_binding
    implicit none

    interface
        subroutine testparticle_boris_push_f64(x, v, q2m, dt, E, B) &
                bind(C, name="testparticle_boris_push_f64")
            import :: c_double
            real(c_double), intent(inout) :: x(3), v(3)
            real(c_double), value :: q2m, dt
            real(c_double), intent(in) :: E(3), B(3)
        end subroutine

        subroutine testparticle_hyper_boris4_push_f64(x, v, q2m, dt, E, B, n) &
                bind(C, name="testparticle_hyper_boris4_push_f64")
            import :: c_double, c_int32_t
            real(c_double), intent(inout) :: x(3), v(3)
            real(c_double), value :: q2m, dt
            real(c_double), intent(in) :: E(3), B(3)
            integer(c_int32_t), value :: n
        end subroutine
    end interface

    real(c_double) :: x(3) = [0.0d0, 0.0d0, 0.0d0]
    real(c_double) :: v(3) = [1.0d5, 0.0d0, 0.0d0]
    real(c_double) :: E(3) = [0.0d0, 0.0d0, 0.0d0]
    real(c_double) :: B(3) = [0.0d0, 0.0d0, 1.0d0]

    call testparticle_boris_push_f64(x, v, 1.0d8, 1.0d-9, E, B)
    call testparticle_hyper_boris4_push_f64(x, v, 1.0d8, 1.0d-9, E, B, 2)

    print *, "New pos:", x
    print *, "New vel:", v
end program
```
