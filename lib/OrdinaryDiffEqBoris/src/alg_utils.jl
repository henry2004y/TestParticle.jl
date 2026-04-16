alg_order(alg::Boris) = 2
isfsal(alg::Boris) = false

alg_order(alg::MultistepBoris{N}) where {N} = 2
isfsal(alg::MultistepBoris{N}) where {N} = false

alg_order(alg::AdaptiveBoris) = 2
isfsal(alg::AdaptiveBoris) = false
