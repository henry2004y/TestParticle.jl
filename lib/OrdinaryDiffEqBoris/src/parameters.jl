"""
    get_q2m(p)

Return the charge-to-mass ratio `q/m` from the parameter container `p`.

The Boris solvers place no constraint on the type of `p`; they only require
`get_q2m`, `get_EField` and `get_BField` to be defined for it. The default
methods below assume the layout `(q2m, m, E, B, ...)`, so any other package
carrying the fields in a different type can support the solvers by adding
methods to these three functions.
"""
get_q2m(p) = p[1]

"""
    get_EField(p)

Return the electric field function `E(x, t)` from the parameter container `p`.
"""
get_EField(p) = p[3]

"""
    get_BField(p)

Return the magnetic field function `B(x, t)` from the parameter container `p`.
"""
get_BField(p) = p[4]
