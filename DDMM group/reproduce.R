# Recompute the 1000G group benchmark from bundled designs and selections.
# Usage: Rscript reproduce.R [--check-only]
args <- commandArgs(TRUE)
if (any(!args %in% '--check-only')) stop('Usage: Rscript reproduce.R [--check-only]')
script <- gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
root <- dirname(normalizePath(script, mustWork = TRUE))
B <- file.path(root, 'benchmark_1000g', 'group_supplement')
if (!requireNamespace('data.table', quietly = TRUE))
  stop("Install the scoring dependency: install.packages('data.table')")

# Check every bundled input, including comparator files which score.R otherwise skips.
manifest <- read.csv(file.path(root, 'input_checksums.csv'), stringsAsFactors = FALSE)
paths <- file.path(root, manifest$path)
missing <- !file.exists(paths)
if (any(missing)) stop('Missing input(s): ', paste(manifest$path[missing], collapse = ', '))
bad <- unname(tools::md5sum(paths)) != manifest$md5
if (any(bad)) stop('Input checksum mismatch: ', paste(manifest$path[bad], collapse = ', '))
for (cohort in c('EUR', 'MIX', 'ALL')) {
  folder <- file.path(B, 'data', cohort)
  d <- readRDS(file.path(folder, 'design.rds'))
  ph <- data.table::fread(file.path(folder, 'phenotypes.tab'))
  ys <- grep('^y_', names(ph), value = TRUE)
  dd <- readRDS(file.path(folder, 'ddmm_selections.rds'))
  su <- readRDS(file.path(folder, 'susie_selections.rds'))
  stopifnot(length(ys) == 60L, nrow(ph) == nrow(d$X), ncol(d$X) == nrow(d$map),
            setequal(ys, names(dd)), setequal(ys, names(su)))
  for (y in ys) for (resolution in sprintf('res%d', 6:0)) {
    replicate <- as.integer(sub('.*_r', '', y))
    discovery <- file.path(B, 'inputs', 'knockoffgwas', cohort, paste0('s', 2019 + replicate), sprintf('%s_%s_discoveries.txt', y, resolution))
    if (!file.exists(discovery)) stop('Missing comparator output: ', discovery)
    for (rule in c('BH', 'truncated'))
      stopifnot(setequal(names(dd[[y]]$resolutions[[resolution]][[rule]]),
                         c('simes', 'rank_effective', 'rank_3', 'rank_full')))
  }
  cat(sprintf('%s: %d individuals, %d variants, %d phenotypes; inputs complete.\n',
              cohort, nrow(d$X), ncol(d$X), length(ys)))
}
if ('--check-only' %in% args) quit(status = 0L)

# These variables apply only to the child analysis processes.
Sys.setenv(BENCH = B, DDMM_PKG = root, OPENBLAS_NUM_THREADS = '1',
           OMP_NUM_THREADS = '1', VECLIB_MAXIMUM_THREADS = '1')
dir.create(file.path(B, 'results'), showWarnings = FALSE)
run <- function(name, arguments = character()) {
  cat('\nRunning ', name, '\n', sep = '')
  status <- system2(file.path(R.home('bin'), 'Rscript'),
                    c(shQuote(file.path(B, 'code', name)), arguments))
  if (status != 0L) stop(name, ' failed with exit status ', status)
}
run('score.R', c('EUR', 'MIX', 'ALL'))
run('describe_design.R', c('EUR', 'MIX', 'ALL'))
run('diag_false.R', c('EUR', 'MIX', 'ALL'))
run('make_tables.R')
run('figure.R')

writeLines(capture.output(sessionInfo()), file.path(B, 'results', 'sessionInfo.txt'))
cat('\nDone. Tables, scores, and figures: ', file.path(B, 'results'), '\n', sep = '')
