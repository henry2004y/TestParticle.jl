module TestParticleRecursiveArrayToolsExt
import TestParticle: trace, trace!
using TestParticle: get_dv
using RecursiveArrayTools: ArrayPartition

# Note: This version works with the `ArrayPartition` structure where `u.x[1]` is position and `u.x[2]` is velocity.
function trace!(du, u::ArrayPartition, p, t)
    v = u.x[2]
    du.x[1] .= v
    du.x[2] .= get_dv(u.x[1], v, p, t)
    return nothing
end

function trace(u::ArrayPartition, p, t)
    v = u.x[2]
    dv = get_dv(u.x[1], v, p, t)
    return vcat(v, dv)
end

end
