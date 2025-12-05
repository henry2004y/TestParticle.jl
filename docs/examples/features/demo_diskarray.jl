# # Interpolate Over Large HDF5 Arrays
#
# This example demonstrates how to interpolate over large HDF5 arrays using `DiskArrays.jl`.
# This is useful when the field data is too large to fit in memory.
# The implementation is based on [this Discourse post](https://discourse.julialang.org/t/interpolate-over-large-hdf5-arrays/127079/12).

using HDF5, DiskArrays
import DiskArrays: eachchunk, haschunks, readblock!, writeblock!, GridChunks, Chunked,
                   Unchunked
using Interpolations

# ## Implement HDF5DiskArray
#
# First, we define a wrapper around HDF5 Dataset that implements the `DiskArrays` interface.

struct HDF5DiskArray{T, N, CS} <: AbstractDiskArray{T, N}
   ds::HDF5.Dataset
   cs::CS
end

Base.size(x::HDF5DiskArray) = size(x.ds)
haschunks(x::HDF5DiskArray{<:Any, <:Any, Nothing}) = Unchunked()
haschunks(x::HDF5DiskArray) = Chunked()
eachchunk(x::HDF5DiskArray{<:Any, <:Any, <:GridChunks}) = x.cs
readblock!(x::HDF5DiskArray, aout, r::AbstractUnitRange...) = aout .= x.ds[r...]
writeblock!(x::HDF5DiskArray, v, r::AbstractUnitRange...) = x.ds[r...] = v

function HDF5DiskArray(ds::HDF5.Dataset)
   chunks = try
      HDF5.get_chunk(ds)
   catch
      nothing
   end
   cs = isnothing(chunks) ? nothing : GridChunks(size(ds), chunks)
   HDF5DiskArray{eltype(ds), ndims(ds), typeof(cs)}(ds, cs)
end

# ## Create Example Data
#
# We create a dummy HDF5 file with some field data.

filename = "example_field.h5"

h5open(filename, "w") do fid
   g = create_group(fid, "mygroup")
   ## Create a dataset with chunking enabled
   dset_a = create_dataset(g, "A", Float32, (10, 10, 10), chunk = (5, 5, 5))
   values_a = [Float32(i + j + k) for i in 1:10, j in 1:10, k in 1:10]
   write(dset_a, values_a)
   println("Wrote dataset A")

   ## Adding Vector Field B (Bx, By, Bz)

   ## Generate data for the vector components
   values_bx = [Float32(sin(i) * j) for i in 1:10, j in 1:10, k in 1:10]
   values_by = [Float32(cos(j) * k) for i in 1:10, j in 1:10, k in 1:10]
   values_bz = [Float32(i * k / 10) for i in 1:10, j in 1:10, k in 1:10]

   ## Create datasets with the same chunking properties
   dset_bx = create_dataset(g, "Bx", Float32, (10, 10, 10), chunk = (5, 5, 5))
   dset_by = create_dataset(g, "By", Float32, (10, 10, 10), chunk = (5, 5, 5))
   dset_bz = create_dataset(g, "Bz", Float32, (10, 10, 10), chunk = (5, 5, 5))

   ## Write the data
   write(dset_bx, values_bx)
   write(dset_by, values_by)
   write(dset_bz, values_bz)

   println("Wrote datasets Bx, By, Bz")
end

println("Successfully created $filename")

# ## Interpolation Function
#
# We define a function `itp` to query the field at a single location `x`.

## Open the file and wrap the dataset
h5open(filename, "r") do fid
   da = HDF5DiskArray(fid["mygroup/A"])

   cached = DiskArrays.cache(da)
   itp = extrapolate(
      interpolate(cached, BSpline(Linear(Periodic(OnCell())))), Periodic(OnCell()))

   ## Evaluate at a point
   loc_int = (5.0, 5.0, 5.0)
   println("Value at $loc_int: ", itp(loc_int...))

   loc_float = (5.5, 5.5, 5.5)
   println("Value at $loc_float: ", itp(loc_float...))

   loc_out = (-0.5, 1.0, 1.0)
   val_periodic = itp(loc_out...)
   println("Value at $loc_out (Periodic): ", val_periodic)
   ## Check correctness of Periodic
   ## Note that we assume cell center values. -0.5 wraps to 9.5 (since period is 10).
   ## Value at 9.5, 1.0, 1.0. Data is i+j+k.
   println("Expected Periodic: ", 9.5 + 1 + 1)
end

##TODO: How to interpolate field vectors efficiently?

# ## Cleanup
# Remove the temporary file.

rm(filename)
