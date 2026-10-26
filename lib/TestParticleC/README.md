# TestParticleC

`TestParticleC` is a lightweight C-ABI wrapper around `BorisPushers.jl` providing
statically callable C symbols for high-performance particle pushing in C, C++, and Fortran
MHD and PIC codes.

## C API (`include/testparticle.h`)

- **Single Particle Push**:
  - `testparticle_boris_push_f64(x, v, q2m, dt, E, B)`
  - `testparticle_boris_push_f32(x, v, q2m, dt, E, B)`
- **Staggered Velocity Push**:
  - `testparticle_boris_step_staggered_f64(x, v_half, q2m, dt, E, B)`
  - `testparticle_boris_step_staggered_f32(x, v_half, q2m, dt, E, B)`
- **Batch Particle Push (Vectorized for PIC)**:
  - `testparticle_boris_push_batch_f64(n, x, v, q2m, dt, E, B)`
  - `testparticle_boris_push_batch_f32(n, x, v, q2m, dt, E, B)`
- **Pure Velocity Rotation**:
  - `testparticle_update_velocity_boris_f64(v, E, B, qdt2m)`
  - `testparticle_update_velocity_boris_f32(v, E, B, qdt2m)`

## Compiling to Shared Library with JuliaC

To build `libtestparticle.so` (Linux), `libtestparticle.dylib` (macOS), or
`libtestparticle.dll` (Windows):

```bash
juliac --output-lib libtestparticle.so \
       --project=lib/TestParticleC \
       lib/TestParticleC/src/TestParticleC.jl
```

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

    testparticle_boris_push_f64(x, v, q2m, dt, E, B);

    printf("New position: (%e, %e, %e)\n", x[0], x[1], x[2]);
    printf("New velocity: (%e, %e, %e)\n", v[0], v[1], v[2]);
    return 0;
}
```

Link against `-ltestparticle` when compiling:
```bash
gcc -O3 main.c -Iinclude -L. -ltestparticle -o main
```

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
    end interface

    real(c_double) :: x(3) = [0.0d0, 0.0d0, 0.0d0]
    real(c_double) :: v(3) = [1.0d5, 0.0d0, 0.0d0]
    real(c_double) :: E(3) = [0.0d0, 0.0d0, 0.0d0]
    real(c_double) :: B(3) = [0.0d0, 0.0d0, 1.0d0]

    call testparticle_boris_push_f64(x, v, 1.0d8, 1.0d-9, E, B)
    print *, "New pos:", x
    print *, "New vel:", v
end program
```
