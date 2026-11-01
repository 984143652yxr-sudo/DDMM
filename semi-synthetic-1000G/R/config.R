# Experiment settings, random-seed rules and the frozen design of a run.

# Panel codes enter every seed. Never reorder or reuse them.
PANEL_CODE <- c(AFR300 = 1L, EAS300 = 2L, EUR300 = 3L, MIX300 = 4L, MIX600 = 5L)

# Files whose content defines the simulated data and the fitted methods.
# Editing any of them changes the design fingerprint, so an existing run
# directory cannot silently mix old and new code.
SCIENCE_CODE <- c('R/data.R', 'R/simulate.R', 'R/methods.R', 'R/knockoffs.R')

EXPERIMENT_ALIASES <- c(ancestry = 'ancestry', ancestry300 = 'ancestry',
                        kinship = 'kinship600', kinship600 = 'kinship600')

# ---------------------------------------------------------------------------
# Experiment configuration (config/<name>.json)
# ---------------------------------------------------------------------------

load_experiment <- function(name) {
  path <- if (grepl('\\.json$', name)) name else {
    key <- EXPERIMENT_ALIASES[tolower(name)]
    if (is.na(key)) stop('Unknown experiment "', name, '". Use ancestry or kinship (alias kinship600).')
    file.path('config', paste0(key, '.json'))
  }
  if (!file.exists(path)) stop('Missing configuration file ', path)
  cfg <- jsonlite::read_json(path, simplifyVector = TRUE)
  cfg$config_file <- path
  required <- c('experiment', 'panels', 'amplitudes', 'null_case', 'methods',
                'projection', 'kinship_effect', 'sigma_g', 'sigma_epsilon',
                'bolt_pc_covariates', 'q', 'susie_L', 'susie_pip',
                'phenotype_seed_base', 'knockoff_seed_base', 'output', 'figure_name',
                'table_caption', 'table_labels', 'plot_labels')
  missing <- setdiff(required, names(cfg))
  if (length(missing)) stop(path, ' lacks: ', paste(missing, collapse = ', '))
  stopifnot(all(cfg$panels %in% names(PANEL_CODE)),
            cfg$projection %in% c('NONE', 'PCA'),
            all(cfg$methods %in% c('DDMM', 'SuSiE', 'Knockoff', 'BOLT')),
            all(is.finite(cfg$amplitudes)), all(cfg$amplitudes > 0))
  cfg
}

# Signal amplitudes plus the null case a = 0 (signal first, as before).
case_amplitudes <- function(cfg) c(cfg$amplitudes, if (isTRUE(cfg$null_case)) 0)

# ---------------------------------------------------------------------------
# Seeds. One stream per panel x replicate x amplitude; none is ever recycled.
#   phenotype : phenotype_seed_base + 1e8 * panel + 1e5 * replicate + 1000 * a
#   pilot Lasso CV folds : phenotype seed + 11
#   knockoff filter (column order, CV folds) : phenotype seed + 22
#   HMM knockoff copy (snpknock2) : knockoff_seed_base + 1e4 * panel + replicate
# ---------------------------------------------------------------------------

phenotype_seed <- function(cfg, panel_name, replicate, amplitude) {
  code <- round(amplitude * 1000)
  stopifnot(replicate >= 1, replicate < 1000, amplitude >= 0, code < 100000,
            abs(amplitude * 1000 - code) < 1e-7)
  seed <- cfg$phenotype_seed_base + PANEL_CODE[[panel_name]] * 100000000 +
    replicate * 100000 + code
  stopifnot(seed <= .Machine$integer.max)
  seed
}

knockoff_seed <- function(cfg, panel_name, replicate)
  as.integer(cfg$knockoff_seed_base + 10000L * PANEL_CODE[[panel_name]] + replicate)

# ---------------------------------------------------------------------------
# Input files read for each panel
# ---------------------------------------------------------------------------

panel_input_files <- function(cfg, panel_name) {
  d <- file.path('data', panel_name)
  files <- file.path(d, c('panel.rds', paste0('cache_', cfg$projection, '.rds')))
  if ('Knockoff' %in% cfg$methods)
    files <- c(files, file.path(d, c('hap.haps', 'hap.legend', 'hap.sample',
                                     'hmm_pcs.txt', 'partitions.txt')),
               'data/genetic_map_chr1.txt')
  if ('BOLT' %in% cfg$methods)
    files <- c(files, file.path(d, c('bolt.bed', 'bolt.bim', 'bolt.fam')),
               if (isTRUE(cfg$bolt_pc_covariates)) file.path(d, 'bolt_pcs.txt'),
               'data/genetic_map.txt')
  unique(files)
}

md5_of <- function(files) {
  missing <- files[!file.exists(files)]
  if (length(missing)) stop('Missing input file(s): ', paste(missing, collapse = ', '))
  as.list(tools::md5sum(files))
}

# ---------------------------------------------------------------------------
# Frozen design (design.json in every run directory)
# ---------------------------------------------------------------------------

DESIGN_FIELDS <- c('experiment', 'panels', 'amplitudes', 'null_case', 'methods',
                   'projection', 'kinship_effect', 'sigma_g', 'sigma_epsilon',
                   'bolt_pc_covariates', 'q', 'susie_L', 'susie_pip',
                   'phenotype_seed_base', 'knockoff_seed_base')

design_record <- function(cfg) {
  inputs <- unique(unlist(lapply(cfg$panels, panel_input_files, cfg = cfg)))
  c(list(layout = 'fits/<panel>/rep<NNNN>/a<amplitude>/{phenotype,fit}.rds'),
    cfg[DESIGN_FIELDS],
    list(snpknock2_options = paste(SNPKNOCK2_OPTIONS, collapse = ' '),
         code_md5 = md5_of(SCIENCE_CODE),
         data_md5 = md5_of(inputs)))
}

as_json <- function(x) jsonlite::toJSON(x, auto_unbox = TRUE, digits = NA, pretty = TRUE)

same_design <- function(a, b) {
  keep <- function(x) x[setdiff(names(x), c('created', 'adopted_from_legacy'))]
  identical(as.character(as_json(keep(a))), as.character(as_json(keep(b))))
}

design_differences <- function(saved, current) {
  keys <- union(names(saved), names(current))
  keys <- setdiff(keys, c('created', 'adopted_from_legacy'))
  keys[!vapply(keys, function(k) identical(as.character(as_json(saved[[k]])),
                                           as.character(as_json(current[[k]]))), logical(1))]
}

# Create design.json, or check that an existing run used the same design.
ensure_design <- function(output, cfg) {
  path <- file.path(output, 'design.json')
  current <- design_record(cfg)
  lock <- file.path(output, '.design.lock')     # several run.R may start together
  acquired <- FALSE
  for (i in 1:150) { if (acquire_lock(lock)) { acquired <- TRUE; break }; Sys.sleep(0.2) }
  if (!acquired) stop('Could not lock ', lock, ' (', describe_lock(lock), ').')
  on.exit(release_lock(lock), add = TRUE)
  if (file.exists(path)) {
    saved <- jsonlite::read_json(path, simplifyVector = TRUE)
    if (!same_design(saved, current))
      stop('The run directory ', output, ' was created with a different design (',
           paste(design_differences(saved, current), collapse = ', '),
           '). Resume with the original code/settings, or choose a new --out directory.')
    return(invisible(saved))
  }
  if (file.exists(file.path(output, 'design.rds'))) {
    current$adopted_from_legacy <- adopt_legacy_design(output, cfg)
  } else if (length(setdiff(list.files(output), c('run.log', '.design.lock')))) {
    stop('The directory ', output, ' is not empty and has no design.json. ',
         'Choose a new --out directory; existing results are never overwritten.')
  }
  current$created <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
  write_lines_atomic(as_json(current), path)
  invisible(current)
}

# ---------------------------------------------------------------------------
# Run directories created by the previous scripts (ancestry_readable.R and
# kinship600.R, now in legacy/) carry design.rds instead of design.json.
# They can be resumed or extended when settings and data are unchanged and
# the old code is the version checked to give identical phenotypes and fits.
# ---------------------------------------------------------------------------

LEGACY_VERIFIED_MD5 <- c(
  'methods.R'       = 'af43d1018ac70bd1d464520f92cf5182',
  'simulation.R'    = '91b171c15d17ff9cfc215f329be74b13',
  'hmm_knockoffs.R' = '43ee550a7a870e0fc507dad66417052b')

adopt_legacy_design <- function(output, cfg) {
  old <- readRDS(file.path(output, 'design.rds'))
  expected <- c(ancestry = 'independent_hmm_v3', kinship = 'kinship600_independent_hmm_v1')
  problems <- character()
  check <- function(ok, what) if (!isTRUE(ok)) problems <<- c(problems, what)
  check(identical(old$version, unname(expected[cfg$experiment])), 'design version')
  check(identical(as.character(old$panels), as.character(cfg$panels)), 'panels')
  check(identical(as.character(old$methods), as.character(cfg$methods)), 'methods')
  check(isTRUE(all.equal(as.numeric(old$amplitude), as.numeric(cfg$amplitudes))), 'amplitudes')
  check(old$seed_base == cfg$phenotype_seed_base, 'seed base')
  for (k in c('q', 'sigma_epsilon', 'susie_L', 'susie_pip'))
    check(isTRUE(all.equal(as.numeric(old[[k]]), as.numeric(cfg[[k]]))), k)
  if (cfg$kinship_effect) check(isTRUE(all.equal(old$sigma_g, cfg$sigma_g)), 'sigma_g')
  code <- old$input_md5[names(LEGACY_VERIFIED_MD5)]
  check(identical(unname(code), unname(LEGACY_VERIFIED_MD5)), 'legacy code version')
  data_files <- grep('^data/', names(old$input_md5), value = TRUE)
  check(all(file.exists(data_files)) &&
          identical(unname(old$input_md5[data_files]), unname(tools::md5sum(data_files))),
        'input data checksums')
  if (length(problems))
    stop('Cannot continue the legacy run in ', output, ': it differs in ',
         paste(problems, collapse = ', '), '. Use a new --out directory.')
  sprintf('design.rds (%s) checked on %s: same settings, data and verified code',
          old$version, format(Sys.Date()))
}
