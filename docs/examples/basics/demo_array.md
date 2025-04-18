# ---
# title: Array Shape Conversion
# id: demo_array
# date: 2025-04-17
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.11.4
# description: F-style and C-style array conversion
# ---

This example shows how to convert between row-major arrays in Python/C and column-major arrays in Julia/Fortran.
Let's say you have created an 4D array in Numpy of size `(nx, ny, nz, 3)`:

```python
import numpy as np

nx, ny, nz = 2, 4, 6
# Calculate the total number of elements
total_elements = nx * ny * nz * 3

# Create a 1D array with the desired range
one_d_array = np.arange(1, total_elements + 1)

# Reshape the 1D array into the desired 4D shape
B = one_d_array.reshape((nx, ny, nz, 3))
# Validation
print(B[0,1,2,1]) # should be 26

file_name = "my_field.npz"
np.savez(file_name, data=B)
```

In Julia, you can load the Numpy array and convert to a Julia-style field array that TestParticle.jl takes:

```julia
using NPZ

vars = npzread("my_field.npz")
B = vars["data"] # size (nx, ny, nz, 3)
B = permutedims(B, (4,1,2,3)) # size (3, nx, ny, nz)
# Validation
println(B[2,1,2,3]) # should be 26
```

Note that the shape of the Julia-style array is now `(3, nx, ny, nz)`, instead of `(3, nz, ny, nx)`. The key thing here is that the array storage ordering is an implementation detail that is not related to how we do math. When we use the `NPZ` package to load the Numpy array, the storage ordering has already been modified. We perform one permutation of the vector dimension to move it to the first index such that the vector components becomes continuous and adapt to the expection of TestParticle.jl.