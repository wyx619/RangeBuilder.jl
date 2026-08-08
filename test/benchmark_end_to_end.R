benchmark_points <- function(n) {
  halton <- function(i, base) {
    value <- 0
    factor <- 1 / base
    while (i > 0) {
      value <- value + (i %% base) * factor
      i <- i %/% base
      factor <- factor / base
    }
    value
  }
  i <- seq_len(n)
  cbind(vapply(i, halton, numeric(1), base=2),
        vapply(i, halton, numeric(1), base=3))
}

median_time <- function(f, repeats=7) {
  f()
  times <- replicate(repeats, unname(system.time(f())["elapsed"]))
  median(times)
}

library(alphahull)
for (n in c(500, 2000)) {
  points <- benchmark_points(n)
  alpha <- 0.07
  elapsed <- median_time(function() {
    hull <- ahull(points, alpha=alpha)
    areaahull(hull)
  })
  cat(sprintf("R ahull + areaahull n=%d: %.4f s\n", n, elapsed))
}
