args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: generate_golden.R <r-source-dir> <output-dir>")
rdir <- normalizePath(args[[1]])
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
library(interp)

for (name in c("rotation", "dummycoor", "delvor", "ashape", "complement", "inter", "lengthahull", "ahull")) {
  source(file.path(rdir, paste0(name, ".R")))
}

datasets <- list(
  square_center = matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0.4, 0.5), ncol = 2, byrow = TRUE),
  concave = matrix(c(0, 0, 2, 0, 2, 2, 1, 0.6, 0, 2, 0.5, 1.1, 1.5, 1.1), ncol = 2, byrow = TRUE),
  ring = { set.seed(20260808); a <- runif(24, 0, 2*pi); r <- sqrt(runif(24, 0.35^2, 1)); cbind(r*cos(a), r*sin(a)) }
)
alphas <- c(0.35, 0.7)

for (name in names(datasets)) {
  for (alpha in alphas) {
    x <- datasets[[name]]
    tag <- sprintf("%s_alpha_%0.2f", name, alpha)
    dv <- delvor(x)
    sh <- ashape(dv, alpha = alpha)
    cp <- complement(dv, alpha = alpha)
    ah <- ahull(dv, alpha = alpha)
    write.csv(x, file.path(outdir, paste0(tag, "_points.csv")), row.names = FALSE)
    write.csv(dv$mesh, file.path(outdir, paste0(tag, "_mesh.csv")), row.names = FALSE)
    write.csv(sh$edges, file.path(outdir, paste0(tag, "_ashape_edges.csv")), row.names = FALSE)
    write.csv(cp, file.path(outdir, paste0(tag, "_complement.csv")), row.names = FALSE)
    write.csv(ah$arcs, file.path(outdir, paste0(tag, "_ahull_arcs.csv")), row.names = FALSE)
    write.csv(ah$xahull, file.path(outdir, paste0(tag, "_ahull_points.csv")), row.names = FALSE)
    write.csv(data.frame(length = ah$length), file.path(outdir, paste0(tag, "_summary.csv")), row.names = FALSE)
  }
}
