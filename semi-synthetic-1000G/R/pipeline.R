# Running replicates: simulate -> fit four methods -> checkpoint.
#
# Run directory layout (one directory per run, e.g. results/ancestry):
#   design.json                     frozen settings, code and data checksums
#   run.log                         every message printed by run.R
#   replicates.csv                  one scored row per panel x replicate x a x method
#   fits/<panel>/rep<NNNN>/
#       knockoff.rds, hmm/          HMM knockoff copy (genotype-only input)
#       a<amplitude>/phenotype.rds  simulated y, beta, u, epsilon, seed (input)
#       a<amplitude>/fit.rds        fitted selections / p-values / PIPs (output)
#       a<amplitude>/bolt/          BOLT-LMM input, logs and association stats
#
# fit.rds is saved after every method, so an interrupted run resumes where it
# stopped. Each replicate is locked while it is processed (fits/<panel>/
# rep<NNNN>.lock); several run.R processes on the same directory therefore
# share the work instead of colliding.

PILOT_SEED_OFFSET <- 11
KNOCKOFF_SEED_OFFSET <- 22

load_panel <- function(panel_name, cfg) {
  list(panel = readRDS(file.path('data', panel_name, 'panel.rds')),
       cache = readRDS(file.path('data', panel_name, paste0('cache_', cfg$projection, '.rds'))))
}

replicate_directory <- function(output, panel_name, replicate)
  file.path(output, 'fits', panel_name, sprintf('rep%04d', replicate))

case_directory <- function(replicate_dir, amplitude)
  file.path(replicate_dir, paste0('a', amplitude))

fit_method <- function(method, cfg, panel, cache, generated, y, pilot, augmented, case_dir) {
  switch(method,
    DDMM = fit_ddmm(cache, pilot, cfg$q),
    SuSiE = fit_susie(cache, y, pilot, L = cfg$susie_L, cutoff = cfg$susie_pip),
    Knockoff = fit_knockoff(augmented, y, generated$seed + KNOCKOFF_SEED_OFFSET, cfg$q),
    BOLT = fit_bolt(panel, generated$y, cfg$bolt_pc_covariates, file.path(case_dir, 'bolt'), cfg$q))
}

# All cases of one panel x replicate. Returns a one-line progress note.
run_replicate <- function(cfg, data, replicate, output, prepare_only = FALSE) {
  panel <- data$panel
  cache <- data$cache
  replicate_dir <- replicate_directory(output, panel$panel, replicate)
  dir.create(replicate_dir, recursive = TRUE, showWarnings = FALSE)
  augmented <- NULL
  knockoff_matrix <- function() {
    if (is.null(augmented)) {
      copy <- get_knockoff(panel, replicate, replicate_dir, cfg)
      augmented <<- cbind(cache$Q %*% panel$D, cache$Q %*% copy)
    }
    augmented
  }
  if ('Knockoff' %in% cfg$methods && prepare_only) invisible(knockoff_matrix())
  fitted <- list()     # amplitude -> methods fitted now
  invalid <- list()    # method -> amplitudes with an invalid fit
  for (amplitude in case_amplitudes(cfg)) {
    case_dir <- case_directory(replicate_dir, amplitude)
    dir.create(case_dir, recursive = TRUE, showWarnings = FALSE)
    phenotype_path <- file.path(case_dir, 'phenotype.rds')
    if (!file.exists(phenotype_path))
      save_rds_atomic(simulate_phenotype(panel, amplitude, replicate, cfg), phenotype_path)
    if (prepare_only) next
    generated <- readRDS(phenotype_path)
    stopifnot(generated$seed == phenotype_seed(cfg, panel$panel, replicate, amplitude))
    fit_path <- file.path(case_dir, 'fit.rds')
    saved <- if (file.exists(fit_path)) readRDS(fit_path) else
      list(panel = panel$panel, replicate = replicate, a = amplitude,
           generated = generated, fits = list())
    pending <- setdiff(cfg$methods, names(saved$fits))
    if (!length(pending)) next
    if ('Knockoff' %in% pending) knockoff_matrix()
    y <- drop(cache$Q %*% generated$y)
    # Fix the generator kind so a resumed case equals an uninterrupted one.
    RNGkind("L'Ecuyer-CMRG")
    pilot <- fit_pilot(cache$X, y, generated$seed + PILOT_SEED_OFFSET)
    for (method in pending) {
      started <- proc.time()[['elapsed']]
      saved$fits[[method]] <- fit_method(method, cfg, panel, cache, generated, y,
                                         pilot, augmented, case_dir)
      saved$seconds[[method]] <- proc.time()[['elapsed']] - started
      save_rds_atomic(saved, fit_path)
    }
    fitted[[format(amplitude)]] <- pending
    for (m in pending) if (!isTRUE(saved$fits[[m]]$valid)) invalid[[m]] <- c(invalid[[m]], amplitude)
  }
  if (prepare_only) return('inputs prepared')
  progress_note(fitted, invalid)
}

# "fitted 10 cases (DDMM/SuSiE/Knockoff/BOLT); invalid: BOLT at a=3"
progress_note <- function(fitted, invalid) {
  if (!length(fitted)) return('already complete')
  sets <- vapply(fitted, paste, '', collapse = '/')
  note <- if (length(unique(sets)) == 1)
    sprintf('fitted %d case%s (%s)', length(sets), if (length(sets) > 1) 's' else '', sets[1])
  else paste0('fitted ', paste0('a=', names(sets), ' (', sets, ')', collapse = ', '))
  if (length(invalid)) note <- paste0(note, '; invalid: ', paste0(names(invalid), ' at a=',
    vapply(invalid, paste, '', collapse = ','), collapse = ', '))
  note
}

run_experiment <- function(cfg, replicates, output, panels = cfg$panels, prepare_only = FALSE) {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(paste0(output, '.lock')))
    stop('A legacy script holds ', paste0(output, '.lock'), '. Wait for it to finish ',
         '(or delete that folder if no R process is running).')
  log_to(file.path(output, 'run.log'))
  ensure_design(output, cfg)
  require_programs(if (prepare_only) intersect(cfg$methods, 'Knockoff') else cfg$methods)
  say(sprintf('%s %s: panels %s, replicates %s, output %s', if (prepare_only) 'Prepare' else 'Run',
              cfg$experiment, paste(panels, collapse = ','), format_ids(replicates), pretty_path(output)))
  skipped <- 0L
  for (panel_name in panels) {
    data <- load_panel(panel_name, cfg)
    for (replicate in replicates) {
      lock <- paste0(replicate_directory(output, panel_name, replicate), '.lock')
      dir.create(dirname(lock), recursive = TRUE, showWarnings = FALSE)
      if (!acquire_lock(lock)) {
        say(sprintf('%s rep %04d  skipped: in use by %s', panel_name, replicate, describe_lock(lock)))
        skipped <- skipped + 1L
        next
      }
      started <- proc.time()[['elapsed']]
      note <- tryCatch(run_replicate(cfg, data, replicate, output, prepare_only),
                       finally = release_lock(lock))
      say(sprintf('%s rep %04d  %s  [%.0f s]', panel_name, replicate, note,
                  proc.time()[['elapsed']] - started))
    }
  }
  if (skipped) say(skipped, ' replicate(s) were being processed by another run.R and were skipped; ',
                   'the process that finishes last writes the complete report.')
  invisible(skipped)
}

format_ids <- function(ids) {
  if (length(ids) > 1 && all(diff(ids) == 1)) sprintf('%d:%d', min(ids), max(ids))
  else paste(ids, collapse = ',')
}

# ---------------------------------------------------------------------------
# Progress of a run directory
# ---------------------------------------------------------------------------

show_status <- function(cfg, output) {
  cat('Run directory:', pretty_path(output), '\n')
  if (!dir.exists(output)) { cat('  (does not exist yet)\n'); return(invisible()) }
  design <- file.path(output, if (file.exists(file.path(output, 'design.json'))) 'design.json' else 'design.rds')
  cat('Design file:  ', if (file.exists(design)) basename(design) else 'none', '\n')
  if (dir.exists(paste0(output, '.lock')))
    cat('NOTE: legacy lock', paste0(basename(output), '.lock'), 'present (a legacy script is running or was killed).\n')
  fit_files <- list.files(file.path(output, 'fits'), pattern = '^fit[.]rds$', recursive = TRUE, full.names = TRUE)
  phenotypes <- list.files(file.path(output, 'fits'), pattern = '^phenotype[.]rds$', recursive = TRUE)
  locks <- list.files(file.path(output, 'fits'), pattern = '[.]lock$', recursive = TRUE, include.dirs = TRUE)
  n_cases <- length(case_amplitudes(cfg))
  cat(sprintf('Cases per replicate: %d (a = %s)\n', n_cases, paste(case_amplitudes(cfg), collapse = ', ')))
  if (!length(fit_files)) {
    cat(sprintf('Prepared phenotypes: %d; no fitted cases yet.\n', length(phenotypes)))
    return(invisible())
  }
  rows <- do.call(rbind, lapply(fit_files, function(path) {
    x <- readRDS(path)
    data.frame(panel = x$panel, replicate = x$replicate, a = x$a,
               complete = all(cfg$methods %in% names(x$fits)),
               invalid = sum(!vapply(x$fits, function(f) isTRUE(f$valid), logical(1))))
  }))
  by_rep <- aggregate(cbind(complete, cases = 1) ~ panel + replicate, rows, sum)
  by_rep$done <- by_rep$complete == n_cases
  overview <- do.call(rbind, lapply(split(by_rep, by_rep$panel), function(z) data.frame(
    panel = z$panel[1], complete_replicates = sum(z$done),
    partial_replicates = sum(!z$done),
    replicate_ids = format_ids(sort(z$replicate[z$done])))))
  overview$invalid_fits <- tapply(rows$invalid, rows$panel, sum)[overview$panel]
  print(overview, row.names = FALSE)
  if (length(locks)) cat('Replicates in progress:', paste(basename(locks), collapse = ' '), '\n')
  invisible(rows)
}
