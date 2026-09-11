alg_order(::Boris) = 2
isfsal(::Boris) = false

alg_order(::MultistepBoris{N}) where {N} = 2
isfsal(::MultistepBoris{N}) where {N} = false

alg_order(::AdaptiveBoris) = 2
isfsal(::AdaptiveBoris) = false

alg_order(::AdaptiveMultistepBoris{N}) where {N} = 2
isfsal(::AdaptiveMultistepBoris{N}) where {N} = false
