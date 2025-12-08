# Random number generation utilities

"""
    rand_gamma(shape; scale=1.0)

Generate a random number from a Gamma distribution with given `shape` (k) and `scale` (θ).
The probability density function is:
\$ f(x) = \\frac{1}{\\Gamma(k) \\theta^k} x^{k-1} e^{-x/\\theta} \$

This implementation uses the method by Marsaglia and Tsang (2000) for \$k \\ge 1\$,
and the property \$Gamma(k, 1) = Gamma(k+1, 1) \\cdot U^{1/k}\$ for \$k < 1\$.
"""
function rand_gamma(shape::Real; scale::Real = 1.0)
   if shape < 1.0
      # Use property: Gamma(shape, 1) = Gamma(shape+1, 1) * U^(1/shape)
      return rand_gamma(shape + 1.0; scale = 1.0) * rand()^(1.0 / shape) * scale
   end

   d = shape - 1.0 / 3.0
   c = 1.0 / sqrt(9.0 * d)
   while true
      x = randn()
      v = 1.0 + c * x
      if v <= 0.0
         continue
      end
      v = v * v * v
      u = rand()
      xsq = x * x
      if u < 1.0 - 0.0331 * xsq * xsq
         return d * v * scale
      end
      if log(u) < 0.5 * xsq + d * (1.0 - v + log(v))
         return d * v * scale
      end
   end
end

"""
    rand_chisq(ν)

Generate a random number from a Chi-squared distribution with `ν` degrees of freedom.
"""
function rand_chisq(ν::Real)
   return rand_gamma(ν / 2.0; scale = 2.0)
end

"""
    rand_gen_normal(β; scale=1.0)

Generate a random number from a Generalized Normal (Exponential Power) distribution with shape `β`.
PDF: \$ f(x) \\propto \\exp(-(|x|/\\alpha)^\\beta) \$ where `α` is the scale.
"""
function rand_gen_normal(β::Real; scale::Real = 1.0)
   # Sample G ~ Gamma(1/β, 1)
   g = rand_gamma(1.0 / β; scale = 1.0)
   # x = scale * g^(1/β) * sign
   s = rand() < 0.5 ? -1.0 : 1.0
   return s * scale * g^(1.0 / β)
end
