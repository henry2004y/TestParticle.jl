# # Interpolate over large HDF5 arrays
#
# This example demonstrates how to interpolate over large HDF5 arrays using `DiskArrays.jl`.
# This is useful when the field data is too large to fit in memory.
# The implementation is based on [this Discourse post](https://discourse.julialang.org/t/interpolate-over-large-hdf5-arrays/127079/12).

using HDF5, DiskArrays
import DiskArrays: eachchunk, haschunks, readblock!, writeblock!, GridChunks, Chunked, Unchunked
using Interpolations

# ## Implement HDF5DiskArray
#
# First, we define a wrapper around HDF5 Dataset that implements the `DiskArrays` interface.

struct HDF5DiskArray{T,N,CS} <: AbstractDiskArray{T,N}
   ds::HDF5.Dataset
   cs::CS
end

Base.size(x::HDF5DiskArray) = size(x.ds)
haschunks(x::HDF5DiskArray{<:Any,<:Any,Nothing}) = Unchunked()
haschunks(x::HDF5DiskArray) = Chunked()
eachchunk(x::HDF5DiskArray{<:Any,<:Any,<:GridChunks}) = x.cs
readblock!(x::HDF5DiskArray, aout, r::AbstractUnitRange...) = aout .= x.ds[r...]
writeblock!(x::HDF5DiskArray, v, r::AbstractUnitRange...) = x.ds[r...] = v

function HDF5DiskArray(ds::HDF5.Dataset)
   chunks = try
      HDF5.get_chunk(ds)
   catch
      nothing
   end
   cs = isnothing(chunks) ? nothing : GridChunks(size(ds), chunks)
   HDF5DiskArray{eltype(ds),ndims(ds),typeof(cs)}(ds,cs)
end

# ## Create Example Data
#
# We create a dummy HDF5 file with some field data.

filename = "example_field.h5"
h5open(filename, "w") do fid
   g = create_group(fid, "mygroup")
   ## Create a dataset with chunking enabled
   dset = create_dataset(g, "myvector", Float32, (10,10,10), chunk=(5,5,5))
   values = [Float32(i + j + k) for i in 1:10, j in 1:10, k in 1:10]
   write(dset, values)
end

# ## Interpolation Function
#
# We define a function to query the field at a single location `x`.
# We implement the interpolation logic manually to read only the necessary chunks.
# This avoids loading the entire array into memory.

function get_value(da, x; order=Linear(), bc=Flat())
   ## Handle Periodic BC by wrapping x
   if bc == Periodic()
      sz = size(da)
      x = map((val, s) -> mod(val - 1, s) + 1, x, sz)
      ## After wrapping, we use Flat extrapolation on the wrapped coordinate
      bc = Flat()
   end

   ## Determine indices needed for interpolation (neighbors)
   inds_raw = map(val -> floor(Int, val):ceil(Int, val), x)

   ## Clamp indices to array bounds to find readable block
   sz = size(da)
   inds_read = map((r, s) -> max(1, min(s, first(r))):max(1, min(s, last(r))), inds_raw, sz)

   ## Read data block
   ## da[inds_read...] reads the block from HDF5 via DiskArrays
   data = da[inds_read...]

   ## Construct interpolator on the loaded block
   ## Handle singleton dimensions (e.g. at boundaries or integer coordinates) by using NoInterp
   orders = map(s -> s==1 ? NoInterp() : BSpline(order), size(data))
   itp = interpolate(data, orders)
   itp = extrapolate(itp, bc)

   ## Calculate local coordinates relative to the loaded block
   local_x = map((val, r) -> val - first(r) + 1, x, inds_read)

   return itp(local_x...)
end

# Open the file and wrap the dataset
fid = h5open(filename, "r")
da = HDF5DiskArray(fid["mygroup/myvector"])

# Evaluate at a point
loc_int = (5.0, 5.0, 5.0)
println("Value at $loc_int: ", get_value(da, loc_int))

loc_float = (5.5, 5.5, 5.5)
println("Value at $loc_float: ", get_value(da, loc_float))

# ## Periodic Boundary Conditions
#
# We can specify boundary conditions, e.g., `Periodic()`.
# Our function handles coordinate wrapping manually.

loc_out = (-0.5, 1.0, 1.0)
val_periodic = get_value(da, loc_out, order=Linear(), bc=Periodic())
println("Value at $loc_out (Periodic): ", val_periodic)

# Check correctness of Periodic:
# -0.5 wraps to 9.5 (since period is 10).
# Value at 9.5, 1.0, 1.0.
# Data is i+j+k.
# 9.5 + 1 + 1 = 11.5.
println("Expected Periodic: ", 9.5 + 1 + 1)

# ## Cleanup
#
# Close the file and remove it.

close(fid)
rm(filename)
