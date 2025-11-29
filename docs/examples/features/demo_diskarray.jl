# # Interpolate over large HDF5 arrays
#
# This example demonstrates how to interpolate over large HDF5 arrays using `DiskArrays.jl` and `DiskArrayTools.jl`.
# This is useful when the field data is too large to fit in memory.
# Note that this approach is optimized for batch processing and might be slow for random access patterns typical in particle tracing.
# However, it provides a viable solution for handling larger-than-memory datasets.
#
# The implementation is based on [this Discourse post](https://discourse.julialang.org/t/interpolate-over-large-hdf5-arrays/127079/12).

using HDF5, DiskArrays
import DiskArrays: eachchunk, haschunks, readblock!, writeblock!, GridChunks, Chunked, Unchunked
using DiskArrayTools: InterpolatedDiskArray
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
   # Create a dataset with chunking enabled
   dset = create_dataset(g, "myvector", Float32, (10,10,10), chunk=(5,5,5))
   values = [Float32(i + j + k) for i in 1:10, j in 1:10, k in 1:10]
   write(dset, values)
end

# ## Interpolation
#
# Now we wrap the dataset and create an interpolated view.

fid = h5open(filename, "r")
da = HDF5DiskArray(fid["mygroup/myvector"])

# Define the target grid for interpolation
# Here we upsample the grid.
target_indices = (range(1,10,91), range(1,10,91), 1:10)
newsize = length.(target_indices)
newchunks = GridChunks(newsize, (50,50,5))

# Create the interpolated disk array
# order=Quadratic() specifies the interpolation order.
di = InterpolatedDiskArray(da, newchunks, target_indices..., order=Quadratic())

# We can perform operations on this interpolated array.
# For example, computing the sum over the last dimension.
# This operation accesses the array in chunks, which is efficient.
result_sum = sum(di, dims=3)

println("Sum of interpolated array (subset): ", result_sum[1:5])

# We can also access specific elements, though random access might be slower.
println("Value at [1,1,1]: ", di[1,1,1])

# ## Cleanup
#
# Close the file and remove it.

close(fid)
rm(filename)
