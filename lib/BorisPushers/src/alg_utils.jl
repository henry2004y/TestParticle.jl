alg_order(::Boris) = 2
isfsal(::Boris) = false

# `N` is the order of the gyrophase correction and therefore of the method:
# N = 2 is the multicycle Boris method, N = 4 and N = 6 the Hyper Boris ones.
alg_order(::MultistepBoris{N}) where {N} = N
isfsal(::MultistepBoris{N}) where {N} = false

alg_order(::AdaptiveBoris) = 2
isfsal(::AdaptiveBoris) = false

alg_order(::AdaptiveMultistepBoris{N}) where {N} = N
isfsal(::AdaptiveMultistepBoris{N}) where {N} = false

alg_order(::Union{Vay, HigueraCary, AdaptiveVay, AdaptiveHigueraCary}) = 2
isfsal(::Union{Vay, HigueraCary, AdaptiveVay, AdaptiveHigueraCary}) = false

SciMLBase.isadaptive(
    ::Union{AdaptiveBoris, AdaptiveMultistepBoris, AdaptiveVay, AdaptiveHigueraCary}
) = true
SciMLBase.isadaptive(::Union{Boris, MultistepBoris, Vay, HigueraCary}) = false
