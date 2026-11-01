# Simulate phenotypes on bundled 1000G genotypes, refit DDMM and SuSiE, and score.
# Rscript run_1000g.R --cohort EUR --replicates 1 --out results/1000g_EUR
script <- normalizePath(gsub('~+~', ' ', sub('^--file=', '',
  grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE))
root <- dirname(script)
source(file.path(root, 'R', 'utils.R'))
cmd <- parse_command_line(commandArgs(TRUE), flags = c('help', 'prepare-only'),
  options = c('cohort', 'replicates', 'amplitudes', 'noise-sd', 'seed', 'q', 'out'))
o <- cmd$options
if (isTRUE(o$help)) {
  cat('Rscript run_1000g.R --cohort EUR --replicates 1 --out results/1000g_EUR\n',
      '--cohort EUR|MIX|ALL     One bundled genotype panel (default EUR)\n',
      '--replicates 1:20       Noise replicate IDs (default 1)\n',
      '--amplitudes 8,4        Positive integer effect multipliers (default 8,4)\n',
      '--noise-sd 1            Gaussian noise SD\n',
      '--seed 1000             Noise seed = this value + replicate ID\n',
      '--q 0.1                 DDMM BH/truncated FDR level\n',
      '--out DIR               Separate run folder; reuse to resume identical settings\n',
      '--prepare-only          Generate phenotypes without fitting\n', sep = '')
  quit(status = 0)
}
if (length(cmd$positional)) stop('Use named options; see --help.')
cohort <- toupper(o$cohort %||% 'EUR')
if (!cohort %in% c('EUR', 'MIX', 'ALL')) stop('--cohort must be EUR, MIX or ALL.')
reps <- parse_replicates(o$replicates %||% '1')
amplitudes <- unique(parse_number_list(o$amplitudes %||% '8,4'))
noise_sd <- as.numeric(o$`noise-sd` %||% '1')
seed <- as.numeric(o$seed %||% '1000')
q <- as.numeric(o$q %||% '0.1')
if (!length(amplitudes) || any(!is.finite(amplitudes) | amplitudes <= 0 |
                             amplitudes != floor(amplitudes))) stop('Use positive integer amplitudes.')
if (!is.finite(noise_sd) || noise_sd <= 0) stop('--noise-sd must be positive.')
if (!is.finite(seed) || seed < 0 || seed != floor(seed) || seed + max(reps) > .Machine$integer.max)
  stop('--seed must be a nonnegative integer within the R seed range.')
if (!is.finite(q) || q <= 0 || q >= 1) stop('--q must be between 0 and 1.')
require_packages(c('data.table', if (!isTRUE(o$`prepare-only`)) c('glmnet', 'susieR')))
base <- file.path(root, 'benchmark_1000g', 'group_supplement')
input <- file.path(base, 'data', cohort)
out <- absolute_path(o$out %||% file.path('results', paste0('1000g_', cohort)), getwd())
files <- c(file.path(input, 'design.rds'), list.files(file.path(input, 'partitions'), full.names = TRUE),
           script, file.path(root, 'R', 'ddmm.R'),
           file.path(base, 'code', c('run_ddmm.R', 'run_susie.R', 'score.R')))
settings <- list(cohort = cohort, replicates = reps, amplitudes = amplitudes,
                 noise_sd = noise_sd, seed = seed, q = q,
                 input_code_md5 = unname(tools::md5sum(files)))
manifest <- file.path(out, 'settings.rds')
if (file.exists(manifest)) {
  if (!identical(readRDS(manifest), settings)) stop('Settings or inputs changed. Choose a new --out directory.')
} else {
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
    stop('Output directory is not empty. Choose a new --out directory.')
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  save_rds_atomic(settings, manifest)
}
work <- file.path(out, 'data', cohort)
dir.create(file.path(work, 'partitions'), recursive = TRUE, showWarnings = FALSE)
parts <- list.files(file.path(input, 'partitions'), full.names = TRUE)
stopifnot(all(file.copy(parts, file.path(work, 'partitions'), overwrite = TRUE)))
d <- readRDS(file.path(input, 'design.rds'))
n <- nrow(d$X)
G <- scale(d$X[, d$causal$SNP, drop = FALSE])
z <- as.numeric(scale(2 - d$samples$sex))
sex_term <- 10 * z / sqrt(n)
pheno <- data.table::data.table(FID = d$samples$sample, IID = d$samples$sample, sex = d$samples$sex)
genetic <- lapply(amplitudes, function(a) as.numeric(G %*% (a * d$causal$Sign * d$causal$Scale / sqrt(n))))
for (r in reps) {
  set.seed(seed + r)
  e <- rnorm(n, sd = noise_sd)
  for (k in seq_along(amplitudes))
    pheno[[sprintf('y_a%d_r%02d', amplitudes[k], r)]] <- round(sex_term + genetic[[k]] + e, 6)
  pheno[[sprintf('y_null_r%02d', r)]] <- round(sex_term + e, 6)
}
d$amplitudes <- amplitudes
d$replicates <- length(reps)
d$h2 <- setNames(vapply(genetic, function(g) var(g) / (var(g) + noise_sd^2 + var(sex_term)), 0), amplitudes)
save_rds_atomic(d, file.path(work, 'design.rds'))
data.table::fwrite(pheno, file.path(work, 'phenotypes.tab'), sep = '\t')
cat(sprintf('%s: %d individuals x %d SNPs; %d newly generated phenotypes\n',
            cohort, n, ncol(d$X), ncol(pheno) - 3L))
if (isTRUE(o$`prepare-only`)) quit(status = 0)
Sys.setenv(BENCH = out, DDMM_PKG = root)
run <- function(name, args) {
  status <- system2(file.path(R.home('bin'), 'Rscript'), shQuote(c(file.path(base, 'code', name), args)))
  if (status != 0) stop(name, ' failed; rerun the same command to resume.')
}
run('run_ddmm.R', c(work, as.character(q)))
run('run_susie.R', c(work, '40'))
run('score.R', cohort)
s <- data.table::fread(file.path(out, 'results', 'summary.csv'))
s <- s[rule %in% c('BH', 'CS95', 'CS90', 'PIP90')]
table <- s[, .(case, resolution, method, reps, power = round(power, 3),
               FDR = round(FDR, 3), any_discovery = round(any_discovery, 3))]
lines <- c('# 1000G simulation results', '',
           sprintf('%s; %d replicates per setting; DDMM q = %g.', cohort, length(reps), q), '',
           'Power = fraction of causal loci detected. FDR = mean false discovery proportion.',
           'Group/credible-set discoveries use genomic interval overlap; summary.csv also includes exact-member FDR.', '',
           paste0('| ', paste(names(table), collapse = ' | '), ' |'),
           paste0('| ', paste(rep('---', ncol(table)), collapse = ' | '), ' |'),
           apply(as.data.frame(table), 1, function(row) paste0('| ', paste(row, collapse = ' | '), ' |')))
writeLines(lines, file.path(out, 'results', 'tables.md'))
writeLines(capture.output(sessionInfo()), file.path(out, 'results', 'sessionInfo.txt'))
cat('\nNew fits: ', work, '\nScores and tables: ', file.path(out, 'results'), '\n', sep = '')
