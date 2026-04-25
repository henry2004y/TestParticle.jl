module TestParticleRecursiveArrayToolsExt

import TestParticle: trace, trace!
using TestParticle: get_dv
using RecursiveArrayTools: ArrayPartition

# Note: This version works with the `ArrayPartition` structure where `u.x[1]` is position and `u.x[2]` is velocity.
@inbounds function trace!(du, u::ArrayPartition, p, t)
    v = u.x[2]
    dv = get_dv(v, u.x[1], p, t)

    for i in eachindex(du.x[1])
        du.x[1][i] = v[i]
    end
    for i in eachindex(du.x[2])
        du.x[2][i] = dv[i]
    end
    return nothing
end

@inbounds function trace(u::ArrayPartition, p, t)
    v = u.x[2]
    dv = get_dv(v, u.x[1], p, t)
    return vcat(v, dv)
end

end
