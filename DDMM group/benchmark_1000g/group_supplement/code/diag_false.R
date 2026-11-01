## Diagnostic of the false discoveries of DDMM (effective rank, BH): for each false group,
## the largest r^2 between its variants and a causal variant, the distance to the nearest
## causal variant, and the largest |r| with a causal variant on the other chromosome
## (a signature of population structure).
##   BENCH=<folder> Rscript diag_false.R EUR MIX ALL -> results/false_discoveries.csv
suppressMessages(library(data.table))
B <- Sys.getenv("BENCH"); COH <- commandArgs(TRUE)
out <- list()
for (coh in COH) {
  dir <- file.path(B, "data", coh); d <- readRDS(file.path(dir, "design.rds")); map <- d$map; p <- nrow(map)
  part <- rbindlist(lapply(c(21, 22), function(ch) fread(file.path(dir, "partitions", sprintf("part_chr%d.txt", ch)))))[match(map$SNP, SNP)]
  ci <- match(d$causal$SNP, map$SNP); X <- scale(d$X * 1.0); X[is.na(X)] <- 0; Rc <- cor(X, X[, ci])
  same <- outer(map$CHR, map$CHR[ci], "=="); dist <- abs(outer(map$BP, map$BP[ci], "-")); dist[!same] <- Inf
  sel <- readRDS(file.path(dir, "ddmm_selections.rds"))
  for (k in 1:7) { r <- sprintf("res%d", 7 - k); col <- sprintf("res_%d", k)
    g <- integer(p); off <- 0L
    for (ch in c(21, 22)) { i <- which(map$CHR == ch); loc <- as.integer(factor(part[[col]][i])); g[i] <- loc + off; off <- off + max(loc) }
    for (yc in grep("_(a8|a4)_", names(sel), value = TRUE)) for (G in sel[[yc]]$resolutions[[r]]$BH$rank_effective) {
      m <- which(g == G); if (any(m %in% ci)) next
      out[[length(out) + 1]] <- data.table(cohort = coh, phenotype = yc, case = sub("^y_(a[48])_.*", "\\1", yc), resolution = r, group = G, size = length(m),
        max_r2_same_chr = max(Rc[m, ][same[m, , drop = FALSE]]^2), dist_kb = min(dist[m, ]) / 1e3,
        max_r2_other_chr = max(c(0, Rc[m, ][!same[m, , drop = FALSE]]^2))) } }
}
D <- rbindlist(out); fwrite(D, file.path(B, "results", "false_discoveries.csv"))
print(D[, .(false = .N, r2_med = round(median(max_r2_same_chr), 2), r2_ge_0.5 = round(mean(max_r2_same_chr >= 0.5), 2), r2_ge_0.2 = round(mean(max_r2_same_chr >= 0.2), 2),
            within_100kb = round(mean(dist_kb <= 100), 2), dist_med_kb = round(median(dist_kb)), other_chr_r2_max = round(max(max_r2_other_chr), 3)), by = .(cohort, case, resolution)])
