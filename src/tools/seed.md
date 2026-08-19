# R-compatible RNG

`r_rng.jl` ports R 4.6's `Mersenne-Twister` and `L'Ecuyer-CMRG` uniform
streams used by `stats::runif()` and `jitter()`.

The implementation follows `R/src/main/RNG.c`:

- R's 50-step `69069 * seed + 1` initialization;
- the Mersenne-Twister 624-word state, twist, and tempering operations;
- the six-word `L'Ecuyer-CMRG` recurrence used by `future`/`furrr` workers;
- R's open-interval Mersenne-Twister `unif_rand()` fixup;
- `runif()` and the uniform-offset form of `jitter()`.

It is kept separate from Julia's `Random.MersenneTwister`. Callers that need
an R-compatible stream should construct `RSeed.RMersenneTwister` explicitly.

`r_rng()` accepts R's complete `.Random.seed` vector. This preserves the
exact stream position after R has already consumed random values:

```julia
rng = RangeBuilder.RSeed.r_rng(r_random_seed)
```

Mersenne-Twister and L'Ecuyer-CMRG states are accepted. A generator is
mutable, so each concurrent worker must use its own instance.
