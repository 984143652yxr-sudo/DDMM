# Corrected DDMM and GhostKnockoff on public AD candidate variants.
# Usage: Rscript run_ad.R [21,22]  (default: all 22 chromosomes)
script <- gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
setwd(dirname(normalizePath(script)))
args <- commandArgs(TRUE)
if ('--help' %in% args) {
  cat('Rscript run_ad.R [chromosomes, e.g. 21,22]\nN_BIC=50000 controls nodewise BIC; outputs go to output/.\n')
  quit(status = 0)
}
suppressPackageStartupMessages({library(data.table); library(Matrix)})
source('ddmm_summary.R')
CHRS <- if (length(args)) as.integer(strsplit(args[1], ',')[[1]]) else 1:22
stopifnot(length(args) <= 1, length(CHRS) > 0, !anyNA(CHRS), all(CHRS %in% 1:22), !anyDuplicated(CHRS))
N_BIC <- as.numeric(Sys.getenv('N_BIC', '50000'))
stopifnot(length(N_BIC) == 1, is.finite(N_BIC), N_BIC > 1)
DDIR <- 'data/realdata_section5'; OUT <- 'output'
dir.create(DDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, showWarnings = FALSE)

# Source: https://github.com/GuJQ5/FWER_Knockoff_filter/tree/main/Section%205
commit <- '7125b66454137aea19c0a1b9ce2697aaa3f15689'
base <- paste0('https://raw.githubusercontent.com/GuJQ5/FWER_Knockoff_filter/', commit, '/Section%205/')
for (name in c('AD_Zscores_candidate.csv', 'Candidate_variants_info_other.csv', sprintf('Cor_%d.csv', CHRS))) {
  path <- file.path(DDIR, name)
  if (file.exists(path) && file.info(path)$size > 0) next
  tmp <- tempfile(tmpdir = DDIR)
  status <- tryCatch(download.file(paste0(base, name), tmp, mode = 'wb', quiet = TRUE),
                     error = function(e) { unlink(tmp); stop(conditionMessage(e)) })
  if (status != 0 || file.info(tmp)$size == 0) { unlink(tmp); stop('Download failed: ', name) }
  if (!file.rename(tmp, path)) stop('Could not save ', path)
}

# GhostKnockoff comparator, following the upstream Functions/GhostKnockoff_hello.R.
is_posdef <- function(A, tol = 1e-9) {
  min(eigen(A, symmetric = TRUE, only.values = TRUE)$values) > tol * 10
}

## Shrink s until (M+1)/M Sigma - diag(s) stays positive definite.
shrink_to_psd <- function(s, G, M) {
  s_eps <- 1e-8
  while (s_eps <= 0.1 && !is_posdef((M + 1) / M * G - diag(s * (1 - s_eps), length(s)))) {
    s_eps <- s_eps * 10
  }
  s * (1 - s_eps)
}

## create.solve_equi_M: s = min((M+1)/M lambda_min, 1) for every variant.
solve_equi_M <- function(Sigma, M) {
  G <- cov2cor(Sigma)
  lambda_min <- min(eigen(G, symmetric = TRUE, only.values = TRUE)$values)
  s <- rep(min((M + 1) / M * lambda_min, 1), nrow(G))
  shrink_to_psd(s, G, M) * diag(Sigma)
}

## create.solve_sdp_M: maximise sum(s) subject to 0 <= s <= 1 and
## (M+1)/M G - diag(s) positive semidefinite.
solve_sdp_M <- function(Sigma, M, gaptol = 1e-6, maxit = 1000) {
  if (!requireNamespace("Rdsdp", quietly = TRUE)) {
    stop("method = 'sdp' needs the Rdsdp package; install it or use method = 'equi'.")
  }
  G <- cov2cor(Sigma)
  p <- ncol(G)
  if (!is_posdef(G)) {
    warning("LD matrix is not positive definite; knockoffs may lose power.", immediate. = TRUE)
  }

  ## Linear constraints s >= 0 and s <= 1, then the positive-definite cone.
  Al1 <- -Matrix::Diagonal(p)
  Al2 <- Matrix::Diagonal(p)
  d_As <- c(diag(p))
  As <- Matrix::Diagonal(length(d_As), x = d_As)
  As <- As[which(Matrix::rowSums(As) > 0), ]

  A <- cbind(Al1, Al2, As)
  C <- matrix(c(rep(0, p), rep(1, p), c((M + 1) / M * G)), 1)
  K <- list(s = p, l = 2 * p)

  sol <- Rdsdp::dsdp(A, rep(1, p), C, K,
                     list(gaptol = gaptol, maxit = maxit, logsummary = 0,
                          outputstats = 0, print = 0))
  if (!identical(sol$STATS$stype, "PDFeasible")) {
    warning("SDP solver returned a non-feasible solution; knockoffs may lose power.")
  }

  s <- pmin(pmax(sol$y, 0), 1)
  s <- shrink_to_psd(s, G, M)
  s[s < 0] <- 0
  if (all(s == 0)) {
    warning("SDP construction failed; knockoffs will have no power.", immediate. = TRUE)
  }
  s * diag(Sigma)
}

## GhostKnockoff.prelim
ghostknockoff_prelim <- function(Sigma, M, method = c("equi", "sdp")) {
  method <- match.arg(method)
  p <- nrow(Sigma)
  SigmaInv <- solve(Sigma)
  s <- if (method == "equi") solve_equi_M(Sigma, M) else solve_sdp_M(Sigma, M)
  diag_s <- diag(s, p)

  Sigma_k <- 2 * diag_s - diag_s %*% SigmaInv %*% diag_s
  V.each <- as.matrix(forceSymmetric(Sigma_k - diag_s))
  V <- matrix(1, M, M) %x% V.each + diag(1, M) %x% diag_s
  V.left <- t(chol(as.matrix(forceSymmetric(V)) + diag(1e-8, M * p)))

  list(P.each = diag(1, p) - diag_s %*% SigmaInv, V.left = V.left, M = M)
}

## GhostKnockoff.fit, method = 'marginal': one draw of M knockoff z-score
## vectors, returned as a p x M matrix.
ghostknockoff_zscores <- function(z, prelim) {
  p <- length(z)
  noise <- matrix(prelim$V.left %*% rnorm(ncol(prelim$V.left)), nrow = p)
  as.vector(prelim$P.each %*% z) + noise
}

## GhostKnockoff.filter: multiple-knockoff statistics and q-values.
ghostknockoff_filter <- function(T_0, T_k) {
  T.temp <- cbind(T_0, T_k)
  T.temp[is.na(T.temp)] <- 0
  M <- ncol(T_k)

  ## Ties go to the knockoff, as in the reference which.max.alt.
  which.max.alt <- function(x) {
    i <- which(x == max(x))
    if (length(i) != 1) i[2] else i[1]
  }
  kappa <- apply(T.temp, 1, which.max.alt) - 1
  tau <- apply(T.temp, 1, max) -
    apply(T.temp, 1, function(x) median(x[-which.max(x)]))

  b <- order(tau, decreasing = TRUE)
  c_0 <- kappa[b] == 0
  n_orig <- cumsum(c_0)                                   # temp_0
  ratio <- (1 + (seq_along(b) - n_orig)) / (M * pmax(1, n_orig))

  q <- rep(1, length(tau))
  if (any(tau[b] > 0)) {
    bound <- max(which(tau[b] > 0))
    for (i in seq_along(b)) {
      idx <- i:min(length(b), bound)
      if (length(idx) > 0) q[b[i]] <- min(ratio[idx]) * c_0[i] + 1 - c_0[i]
    }
    q[q > 1] <- 1
  }
  list(kappa = kappa, tau = tau, q = q)
}

set.seed(1)
s0 <- 0.10; M <- 5
Z <- fread(file.path(DDIR, 'AD_Zscores_candidate.csv'))
GA <- fread(file.path(DDIR, 'Candidate_variants_info_other.csv'))
setnames(Z, c('hg38.chr', 'hg38.pos'), c('chr', 'pos'))
gene_of <- function(ch, ps) { g <- GA$HH_ref[GA$chr == ch & GA$pos == ps]; if (length(g)) g[1] else NA_character_ }
meta <- list()
for (ch in CHRS) {
  S <- fread(file.path(DDIR, sprintf('Cor_%d.csv', ch)))
  index <- S[[1]]; Sigma <- as.matrix(S[, -1, with = FALSE])
  Zc <- Z[chr == ch][index]
  z <- Zc$Z; z[!is.finite(z)] <- 0
  Sigma <- as.matrix(forceSymmetric((1 - s0) * Sigma + s0 * diag(nrow(Sigma))))
  pv <- ddmm_pvalues(Sigma, z, n = 1, q = 0, n_bic = N_BIC)
  prelim <- ghostknockoff_prelim(Sigma, M, method = 'equi')
  Zk <- ghostknockoff_zscores(z, prelim)
  gk <- ghostknockoff_filter(z^2, Zk^2)
  meta[[length(meta) + 1]] <- data.table(chr = ch, pos = Zc$pos, Z = z,
      p_ddmm = as.numeric(pv), T_ddmm = attr(pv, 'T'), kappa = gk$kappa, tau = gk$tau)
  cat('chr', ch, ':', length(z), 'variants;', sum(pv < 1), 'DDMM tested\n')
}
MM <- rbindlist(meta); setorder(MM, chr, pos)
MM[, gene := mapply(gene_of, chr, pos)]
ko_sel <- function(kappa, tau, level) {
  thresholds <- sort(unique(tau[kappa == 0 & tau > 0]))
  good <- vapply(thresholds, function(t)
    (1/M) * (1 + sum(kappa >= 1 & tau >= t)) / max(1, sum(kappa == 0 & tau >= t)) <= level, logical(1))
  if (!any(good)) return(logical(length(tau)))
  kappa == 0 & tau >= min(thresholds[good])
}
n_loci <- function(sel) {
  if (!any(sel)) return(0L)
  d <- MM[sel, .(chr, pos)][order(chr, pos)]
  sum(c(TRUE, diff(d$pos) > 1e6 | diff(d$chr) != 0))
}
summ <- list()
for (q in c(0.05, 0.10)) {
  dd <- p.adjust(MM$p_ddmm, 'BH') <= q
  gk <- ko_sel(MM$kappa, MM$tau, q)
  set(MM, j = sprintf('ddmm_q%02d', 100 * q), value = dd)
  set(MM, j = sprintf('gk_q%02d', 100 * q), value = gk)
  for (method in c('DDMM', 'GhostKnockoff')) {
    selected <- if (method == 'DDMM') dd else gk
    summ[[length(summ) + 1]] <- data.table(q, method, variants = sum(selected),
      loci = n_loci(selected), shared_with_GK = sum(selected & gk))
  }
}
SUM <- rbindlist(summ); print(SUM)
fwrite(MM, file.path(OUT, 'ad_pvalues.csv'))
fwrite(SUM, file.path(OUT, 'ad_summary.csv'))
fwrite(MM[ddmm_q05 | ddmm_q10 | gk_q05 | gk_q10], file.path(OUT, 'ad_discoveries.csv'))

# Two-panel comparison; highlight discoveries at q = 0.05.
MM[, x := seq_len(.N)]
bnd <- MM[, .(mid = mean(x)), by = chr]
plot_panel <- function(y, selected, label, title, cap) {
  y <- pmin(y, cap)
  plot(MM$x, y, col = ifelse(MM$chr %% 2, '#9aa7b5', '#d1d8df'), pch = 19,
       cex = 0.4, xaxt = 'n', ylim = c(0, cap * 1.15), xlab = '', ylab = label, main = title)
  points(MM$x[selected], y[selected], col = '#c9252d', pch = 19, cex = 0.8)
  axis(1, at = bnd$mid, labels = bnd$chr, tick = FALSE)
}
png(file.path(OUT, 'ad_comparison.png'), width = 1700, height = 1100, res = 150)
par(mfrow = c(2, 1), mar = c(3, 5, 3, 1))
plot_panel(-log10(pmax(MM$p_ddmm, 1e-300)), MM$ddmm_q05, expression(-log[10](p)),
           sprintf('Corrected DDMM: %d variants', sum(MM$ddmm_q05)), 30)
plot_panel(MM$tau, MM$gk_q05, expression(tau),
           sprintf('GhostKnockoff: %d variants', sum(MM$gk_q05)), 40)
dev.off()

# Locus support at DDMM q = 0.10. Same 1.2 Mb grouping / 1 Mb matching rules as before.
# Published GhostKnockoff loci: He et al. (2022), Supplementary Table 2.
he <- data.table(
  chr = c(1, 1, 1, 2, 2, 4, 5, 6, 6, 6, 8, 10, 10, 11, 11, 11, 14, 14, 15, 15,
          16, 16, 16, 17, 17, 17, 17, 17, 17, 19, 19, 19, 20, 21),
  pos = c(155419060, 161185602, 178992361, 44026309, 127135234, 11026080,
          87002714, 32637301, 40783137, 47479305, 27598736, 11678621, 29966853,
          47440232, 60212842, 86089237, 52710264, 92469490, 58889786, 63277703,
          17478817, 31121341, 70258841, 388402, 49391824, 58320645, 60239372,
          63482562, 73384739, 1046077, 44908684, 51224706, 56443204, 26161943),
  gene = c("POU5F1P4", "ADAMTS4", "FAM20B", "LRPPRC", "BIN1", "CLNK", "MIR4280",
           "HLA-DQA1", "LOC101929555", "CD2AP", "CLU", "ECHDC3", "JCAD",
           "RAPSN", "MIR6503", "PICALM", "PSMC6", "SLC24A4", "SLTM", "APH1B",
           "XYLT1", "KAT8", "AARS1", "LOC105371430", "LOC102724596",
           "TSPOAP1-AS1", "SCARNA20", "ACE", "SDK2", "ABCA7", "APOE", "CD33",
           "CASS4", "APP"))

## Schwartzentruber aggregates Kunkle with UK Biobank, and ADGC and ADSP are the
## in-house sets, so none of the three is independent of Kunkle.
studies <- c(Jansen = "Jansenetal_2019sept", Kunkle = "Kunkle_etal_Stage1",
             Schwartz = "Schwartzentruber_etal", Bis = "Bis_etal_WES",
             ADGC = "Stanford2021_GWAS_AD_ADGC",
             ADSP_WES = "Stanford2021_GWAS_AD_ADSP_WES",
             ADSP_WGS = "Stanford2021_GWAS_AD_ADSP_WGS",
             LeGuan = "LeGuan_etal_WES", IGAP = "IGAP_survival")
independent <- c("Jansen", "Kunkle", "Bis", "LeGuan", "IGAP")

nearest <- function(chr_i, pos_i, ref) {
  same <- ref[chr == chr_i]
  if (!nrow(same)) Inf else min(abs(same$pos - pos_i))
}
cand <- Z[chr %in% CHRS]
he <- he[chr %in% CHRS]
dd <- MM[ddmm_q10 == TRUE]
if (nrow(dd)) {
  dd[, locus := cumsum(c(TRUE, diff(pos) > 1.2e6)), by = chr]
  dd <- dd[, .SD[which.max(abs(Z))], by = .(chr, locus)]
}
he[, covered := mapply(nearest, chr, pos, MoreArgs = list(ref = cand)) <= 1e6]
he[, replicated := mapply(nearest, chr, pos, MoreArgs = list(ref = dd)) <= 1e6]
fwrite(he, file.path(OUT, 'ad_he_loci_replication.csv'))
if (nrow(dd)) {
  support <- as.matrix(cand[dd[, .(chr, pos)], on = .(chr, pos)][, ..studies])
  colnames(support) <- names(studies)
  support[is.na(support) | support == 0] <- NA
  loci <- data.table(chr = dd$chr, pos = dd$pos, gene = dd$gene, Z = round(dd$Z, 2),
     n_studies = rowSums(abs(support) >= 1.96, na.rm = TRUE),
     n_independent = rowSums(abs(support[, independent, drop = FALSE]) >= 1.96, na.rm = TRUE))
} else loci <- data.table(chr = integer(), pos = numeric(), gene = character(), Z = numeric(),
                         n_studies = integer(), n_independent = integer())
fwrite(loci, file.path(OUT, 'ad_ddmm_locus_support.csv'))
writeLines(c(paste('Chromosomes:', paste(CHRS, collapse = ',')), paste('N_BIC:', N_BIC),
             paste('Data commit:', commit), capture.output(sessionInfo())), file.path(OUT, 'run_info.txt'))
cat('Done. Results: ', normalizePath(OUT), '\n', sep = '')
