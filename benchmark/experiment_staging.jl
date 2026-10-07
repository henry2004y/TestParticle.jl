# Experiment: is pinning the host staging buffer worth wiring into EnsembleKernel?
#
# With `raw_output = true` the remaining host cost is 4 to 12 ns per particle, which at
# N = 10^6 is 4 to 12 ms. That is two transfers of the N x 6 state array, one in and one
# out. `zeros` gives pageable memory, which CUDA copies through an internal staging
# buffer at roughly half the bandwidth of pinned memory. This measures the difference at
# the sizes the benchmark actually uses, so we know whether the extra complexity pays.
#
# Run on a machine with an NVIDIA GPU:
#   julia --project=test -e "using CUDA; include(\"benchmark/experiment_staging.jl\")"

using Printf
using Statistics: median

const SIZES = [(1_000_000, 6), (5_000_000, 6)]

"""
    pinned_array(T, dims)

Try to obtain a pinned host array of `T` and `dims`. The CUDA.jl spelling for this has
changed across versions, so every known route is attempted and the one that works is
reported. Returns `nothing` when none of them does.
"""
function pinned_array(T, dims)
    nbytes = prod(dims) * sizeof(T)
    CUDA = Main.CUDA

    if isdefined(CUDA, :alloc) && isdefined(CUDA, :HostMemory)
        try
            buf = CUDA.alloc(CUDA.HostMemory, nbytes)
            ptr = Ptr{T}(UInt(pointer(buf)))
            wrapped = unsafe_wrap(Array{T}, ptr, dims; own = false)
            return (; array = wrapped, route = "alloc(HostMemory)")
        catch err
            println("  route alloc(HostMemory) failed: ", err)
        end
    end

    if isdefined(CUDA, :pin)
        try
            a = zeros(T, dims)
            CUDA.pin(a)
            return (; array = a, route = "pin")
        catch err
            println("  route pin failed: ", err)
        end
    end

    if isdefined(CUDA, :Mem)
        if isdefined(CUDA.Mem, :alloc) && isdefined(CUDA.Mem, :Host)
            try
                buf = CUDA.Mem.alloc(CUDA.Mem.Host, nbytes)
                ptr = Ptr{T}(UInt(pointer(buf)))
                wrapped = unsafe_wrap(Array{T}, ptr, dims; own = false)
                return (; array = wrapped, route = "Mem.alloc(Mem.Host)")
            catch err
                println("  route Mem.alloc(Mem.Host) failed: ", err)
            end
        end

        if isdefined(CUDA.Mem, :pin)
            try
                a = zeros(T, dims)
                CUDA.Mem.pin(a)
                return (; array = a, route = "Mem.pin")
            catch err
                println("  route Mem.pin failed: ", err)
            end
        end
    end

    return nothing
end

bandwidth_gbs(nbytes, seconds) = nbytes / (seconds * 1.0e9)

function median_transfer(f, repeats = 5)
    f()
    Main.CUDA.synchronize()
    return median([@elapsed begin
        f()
        Main.CUDA.synchronize()
    end for _ in 1:repeats])
end

function run()
    isdefined(Main, :CUDA) || error("load CUDA first: using CUDA")
    CUDA = Main.CUDA
    CUDA.functional() || error("no working CUDA device")

    println("Device: ", CUDA.name(CUDA.device()))
    println("Pin routes are tried in order; the first that works is used.\n")

    for T in (Float32, Float64)
        for dims in SIZES
            N = dims[1]
            dev = CUDA.zeros(T, dims)
            pageable = zeros(T, dims)

            h2d_p = median_transfer() do
                copyto!(dev, pageable)
            end
            d2h_p = median_transfer() do
                copyto!(pageable, dev)
            end

            pinned = pinned_array(T, dims)
            pinned === nothing && continue

            h2d_pin = median_transfer() do
                copyto!(dev, pinned.array)
            end
            d2h_pin = median_transfer() do
                copyto!(pinned.array, dev)
            end

            nbytes = prod(dims) * sizeof(T)
            total_p = h2d_p + d2h_p
            total_pin = h2d_pin + d2h_pin

            @printf(
                "%-8s N = %-9d (%5.1f MB each way), pin route: %s\n",
                T, N, nbytes / 1.0e6, pinned.route
            )
            @printf(
                "    pageable  H2D %7.2f ms (%5.1f GB/s) | D2H %7.2f ms (%5.1f GB/s)\n",
                h2d_p * 1.0e3, bandwidth_gbs(nbytes, h2d_p),
                d2h_p * 1.0e3, bandwidth_gbs(nbytes, d2h_p)
            )
            @printf(
                "    pinned    H2D %7.2f ms (%5.1f GB/s) | D2H %7.2f ms (%5.1f GB/s)\n",
                h2d_pin * 1.0e3, bandwidth_gbs(nbytes, h2d_pin),
                d2h_pin * 1.0e3, bandwidth_gbs(nbytes, d2h_pin)
            )
            @printf(
                "    both ways: %.2fx faster | host cost %.1f -> %.1f ns/particle\n\n",
                total_p / total_pin, total_p * 1.0e9 / N, total_pin * 1.0e9 / N
            )
        end
    end

    println(
        "The last line is what matters: `a` in section [4] of the benchmark is the\n" *
            "sum of both transfers plus initialization, so compare the per-particle\n" *
            "figures against the 4.3 ns (Float32) and 12.3 ns (Float64) measured there."
    )
    return nothing
end

run()
