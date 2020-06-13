# To a compute a p.value for a two-sided Zi = |3.49|:

# Do either:
(p.value = pnorm(-3.49) + pnorm(3.49, lower.tail = FALSE))


# Or:
(p.value = 2*pnorm(-3.49))