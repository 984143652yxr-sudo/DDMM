## ============================================================
##  score.R — scores of the group benchmark (cohorts EUR, MIX, ALL).
##  Rule of the KnockoffGWAS paper (Sesia et al. 2021). A discovery is the genomic interval
##  spanned by its variants, [smallest, largest position] within a chromosome. It is true
##  when the interval contains a causal variant. Locus power is the share of the 20 loci
##  with a causal variant inside some discovery, and FDP is the share of discoveries with no
##  causal variant inside (the FDP of a null phenotype is 1{any discovery}).
##  Methods: DDMM (Simes and the chi-squared group test at the effective rank, K = 3 and
##  full rank; BH and the truncated rule), KnockoffGWAS (knockoff+), and SuSiE as a
##  reference without a fixed resolution (credible sets at coverage 0.95 and 0.90 and
##  single variants with PIP > 0.9; for credible sets, which need not be contiguous, the
##  interval rule above is used, and FDP_member counts a set as true only when it holds a
##  causal variant).
##  Diagnostic. A false discovery is an LD tag when one of its variants has r^2 >= 0.8 with
##  a causal variant (de Bakker et al. 2005).
##    BENCH=<folder> Rscript score.R EUR MIX ALL  -> results/records.csv, summary.csv, contrasts.csv
## ============================================================
suppressMessages(library(data.table))
script_dir <- function() { a <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(gsub('~+~', ' ', sub('^--file=', '', a[1]), fixed = TRUE))) else 'code' }
B <- Sys.getenv('BENCH', normalizePath(file.path(script_dir(), '..')))
COHORTS <- commandArgs(TRUE); if (!length(COHORTS)) COHORTS <- c('EUR', 'MIX', 'ALL')
TAG_R2 <- 0.8; N_LOCI <- 20L

read_assignment <- function(dir, map) {
  part <- rbindlist(lapply(c(21, 22), function(ch) fread(file.path(dir, 'partitions', sprintf('part_chr%d.txt', ch)))))
  part <- part[match(map$SNP, part$SNP)]; stopifnot(!anyNA(part$SNP))
  setNames(lapply(1:7, function(k) {
    col <- sprintf('res_%d', k); g <- integer(nrow(map)); off <- 0L
    for (ch in c(21, 22)) { i <- which(map$CHR == ch); loc <- as.integer(factor(part[[col]][i])); g[i] <- loc + off; off <- off + max(loc) }
    g }), sprintf('res%d', 6:0))
}

## score a list of discoveries, each a vector of variant indices
score_sets <- function(sets, ctx, null_case) {
  nd <- length(sets)
  if (!nd) return(data.table(discoveries = 0L, TP = 0L, FDP = 0, FDP_member = 0, loci_found = 0L,
                             power = if (null_case) NA_real_ else 0, tags = 0L, FDP_untagged = 0,
                             width_kb = NA_real_, size = NA_real_))
  width <- vapply(sets, function(s) max(vapply(split(ctx$pos[s], ctx$chr[s]), function(b) diff(range(b)), 0)) / 1e3, 0)
  size <- lengths(sets)
  if (null_case) return(data.table(discoveries = nd, TP = 0L, FDP = 1, FDP_member = 1, loci_found = 0L, power = NA_real_,
                                   tags = 0L, FDP_untagged = 1, width_kb = median(width), size = median(size)))
  ## loci hit by each discovery under the interval rule
  hit <- t(vapply(sets, function(s) {
    h <- logical(N_LOCI)
    for (ch in unique(ctx$chr[s])) { b <- range(ctx$pos[s][ctx$chr[s] == ch])
      inside <- ctx$causal$CHR == ch & ctx$causal$BP >= b[1] & ctx$causal$BP <= b[2]
      h[unique(ctx$causal$Locus[inside])] <- TRUE }
    h }, logical(N_LOCI)))
  if (N_LOCI == 1) hit <- t(hit)
  true_iv <- rowSums(hit) > 0
  true_mem <- vapply(sets, function(s) any(ctx$is_causal[s]), NA)
  tag <- vapply(seq_along(sets), function(i) !true_iv[i] && any(ctx$tag[sets[[i]]]), NA)
  data.table(discoveries = nd, TP = sum(true_iv), FDP = mean(!true_iv), FDP_member = mean(!true_mem),
             loci_found = sum(colSums(hit) > 0), power = sum(colSums(hit) > 0) / N_LOCI,
             tags = sum(tag), FDP_untagged = mean(!true_iv & !tag), width_kb = median(width), size = median(size))
}

recs <- list()
add <- function(...) recs[[length(recs) + 1]] <<- data.table(...)
for (coh in COHORTS) {
  dir <- file.path(B, 'data', coh)
  design <- readRDS(file.path(dir, 'design.rds')); map <- design$map; p <- nrow(map)
  assignment <- read_assignment(dir, map)
  ci <- match(design$causal$SNP, map$SNP)
  Xs <- scale(design$X * 1.0); Xs[is.na(Xs)] <- 0
  r2 <- cor(Xs, Xs[, ci])^2; r2[outer(map$CHR, map$CHR[ci], `!=`)] <- 0
  ctx <- list(chr = map$CHR, pos = map$BP, is_causal = seq_len(p) %in% ci, tag = apply(r2, 1, max) >= TAG_R2,
              causal = data.frame(CHR = design$causal$CHR, BP = design$causal$BP, Locus = design$causal$Locus))
  members <- lapply(assignment, function(a) split(seq_len(p), a))
  sel <- if (file.exists(f <- file.path(dir, 'ddmm_selections.rds'))) readRDS(f) else list()
  sus <- if (file.exists(f <- file.path(dir, 'susie_selections.rds'))) readRDS(f) else list()
  ycols <- union(names(sel), names(sus))
  for (yc in ycols) {
    pr <- regmatches(yc, regexec('^y_(a[0-9]+|null)_r([0-9]+)$', yc))[[1]]; case <- pr[2]; rep <- as.integer(pr[3])
    nul <- case == 'null'
    for (r in names(assignment)) {
      m <- max(assignment[[r]])
      if (!is.null(sel[[yc]])) for (rule in c('BH', 'truncated')) for (t in names(sel[[yc]]$resolutions[[r]][[rule]]))
        add(cohort = coh, case = case, replicate = rep, resolution = r, groups = m, method = paste0('DDMM_', t), rule = rule,
            score_sets(members[[r]][sel[[yc]]$resolutions[[r]][[rule]][[t]]], ctx, nul))
      f <- file.path(B, 'inputs', 'knockoffgwas', coh, paste0('s', 2019 + rep), sprintf('%s_%s_discoveries.txt', yc, r))
      if (file.exists(f)) {
        d <- fread(f)
        sets <- if (nrow(d)) lapply(seq_len(nrow(d)), function(i) which(map$CHR == d$CHR[i] & map$BP >= d$BP.min[i] & map$BP <= d$BP.max[i])) else list()
        add(cohort = coh, case = case, replicate = rep, resolution = r, groups = m, method = 'KnockoffGWAS', rule = 'knockoff+',
            score_sets(sets, ctx, nul))
      }
    }
    if (!is.null(sus[[yc]])) for (s in c('CS95', 'CS90', 'PIP90')) {
      sets <- if (s == 'PIP90') as.list(sus[[yc]]$PIP90) else lapply(sus[[yc]][[s]], as.integer)
      add(cohort = coh, case = case, replicate = rep, resolution = 'adaptive', groups = NA_integer_, method = paste0('SuSiE_', s),
          rule = s, score_sets(sets, ctx, nul))
    }
  }
}
records <- rbindlist(recs, fill = TRUE)
dir.create(file.path(B, 'results'), showWarnings = FALSE)
fwrite(records, file.path(B, 'results', 'records.csv'))

boot <- function(x, nb = 2000) { x <- x[is.finite(x)]; if (length(x) < 2) return(c(NA_real_, NA_real_)); set.seed(1)
  quantile(replicate(nb, mean(sample(x, replace = TRUE))), c(.025, .975), names = FALSE) }
summ <- records[, {
  sig <- case[1] != 'null'
  pw <- if (sig) boot(power) else c(NA_real_, NA_real_); fd <- boot(FDP)
  ci <- binom.test(sum(discoveries > 0), .N)$conf.int
  .(reps = .N, power = if (sig) mean(power) else NA_real_, power_lo = pw[1], power_hi = pw[2],
    FDR = mean(FDP), FDR_lo = fd[1], FDR_hi = fd[2], FDR_member = mean(FDP_member), FDR_untagged = mean(FDP_untagged),
    discoveries = mean(discoveries), tags = mean(tags), width_kb = median(width_kb, na.rm = TRUE), size = median(size, na.rm = TRUE),
    any_discovery = mean(discoveries > 0), any_lo = ci[1], any_hi = ci[2])
}, by = .(cohort, case, resolution, groups, method, rule)]
fwrite(summ, file.path(B, 'results', 'summary.csv'))

## paired power differences on the same phenotypes: DDMM effective rank (BH) minus
## KnockoffGWAS at the same resolution, and minus SuSiE CS95
w <- dcast(records[case != 'null' & ((method == 'DDMM_rank_effective' & rule == 'BH') | method == 'KnockoffGWAS')],
           cohort + case + replicate + resolution ~ method, value.var = 'power')
if (!'KnockoffGWAS' %in% names(w)) w[, KnockoffGWAS := NA_real_]
w <- w[!is.na(KnockoffGWAS) & !is.na(DDMM_rank_effective)]
pair <- function(a, b) list(reps = length(a), diff = mean(a - b), se = sd(a - b) / sqrt(length(a)),
                         ahead = sum(a > b), tied = sum(a == b), behind = sum(a < b))
contr <- w[, pair(DDMM_rank_effective, KnockoffGWAS), by = .(cohort, case, resolution)][, versus := 'KnockoffGWAS']
su <- records[case != 'null' & method == 'SuSiE_CS95', .(cohort, case, replicate, susie = power)]
ws <- merge(records[case != 'null' & method == 'DDMM_rank_effective' & rule == 'BH', .(cohort, case, replicate, resolution, ddmm = power)],
            su, by = c('cohort', 'case', 'replicate'))
contr <- rbind(contr, ws[, pair(ddmm, susie), by = .(cohort, case, resolution)][, versus := 'SuSiE_CS95'])
fwrite(contr, file.path(B, 'results', 'contrasts.csv'))

options(width = 220)
print(summ[rule %in% c('BH', 'knockoff+', 'CS95', 'CS90', 'PIP90') & (method %in% c('DDMM_rank_effective', 'KnockoffGWAS') | grepl('SuSiE', method)),
           .(cohort, case, resolution, groups, method, reps, power = round(power, 3), FDR = round(FDR, 3), any = round(any_discovery, 2), width = round(width_kb, 1))])
