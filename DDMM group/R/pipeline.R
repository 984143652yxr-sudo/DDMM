# Running cases: phenotype -> per-variant test -> group tests at every
# resolution -> BH at level q within each resolution -> scoring.
#
# Run directory (e.g. results/knockoffzoom):
#   design.json        frozen settings, ranks, code and data checksums
#   run.log            messages of every run.R call
#   records.csv        one scored row per case x resolution x track
#   fits/shipped.rds   the published phenotype
#   fits/rep0001.rds   replicate 1: its signal and null cases
# A fits file holds, per case, sigma, the tested variants, every p-value and
# the selected groups, so a report never refits anything.

CODE_FILES <- c('R/data.R', 'R/partitions.R', 'R/ddmm.R', 'R/scoring.R')

# ---------------------------------------------------------------------------
# Experiment configuration
# ---------------------------------------------------------------------------

load_experiment <- function(name) {
  path <- if (grepl('\\.json$', name)) name else file.path('config', paste0(tolower(name), '.json'))
  if (!file.exists(path)) stop('Unknown experiment "', name, '". Use knockoffgwas, or the path of a config ',
                               'file such as config/knockoffgwas.json.')
  cfg <- jsonlite::read_json(path, simplifyVector = FALSE)
  cfg$chromosomes <- unlist(cfg$chromosomes)
  cfg$config_file <- path
  cfg
}

resolution_names <- function(cfg) vapply(cfg$resolutions, `[[`, '', 'name')

# Default report directory of a run directory: reports/<run name>, or the
# sibling reports/ folder when the run directory itself is called results
# (as in legacy/knockoffzoom/results).
default_report_dir <- function(input, root)
  if (basename(input) == 'results') file.path(dirname(input), 'reports') else file.path(root, 'reports', basename(input))

# ---------------------------------------------------------------------------
# Context: genotypes, partitions, group operators, nodewise cache
# ---------------------------------------------------------------------------

nodewise_cache_file <- function(cfg) {
  beds <- file.path(cfg$data_dir, 'data/genotypes', sprintf('example_chr%d.bed', cfg$chromosomes))
  key <- substr(paste(tools::md5sum(beds), collapse = ''), 1, 12)
  file.path('cache', sprintf('nodewise_%s_path%d.rds', key, NODEWISE_PATH$nlambda))
}

load_nodewise <- function(cache, file) {
  if (!file.exists(file)) return(0L)
  saved <- readRDS(file)
  have <- which(!vapply(saved, is.null, logical(1)))
  for (j in have) if (is.null(cache$fits[[j]])) cache$fits[[j]] <- saved[[j]]
  length(have)
}

# Merge this process's nodewise fits into the shared cache file.
save_nodewise <- function(cache, file) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  lock <- paste0(file, '.lock')
  for (i in 1:150) { if (acquire_lock(lock)) break; Sys.sleep(0.2) }
  on.exit(release_lock(lock), add = TRUE)
  load_nodewise(cache, file)
  save_rds_atomic(cache$fits, file)
}

build_context <- function(cfg, ranks, q) {
  say('Loading genotypes and building partitions (', cfg$experiment, ')')
  geno <- load_genotypes(cfg)
  truth <- load_truth(cfg, geno)
  parts <- build_partitions(cfg, geno)
  started <- proc.time()[['elapsed']]
  operators <- lapply(parts$resolutions, function(res)
    if (res$single_variant) NULL else group_operators(geno, res$assignment))
  say(sprintf('Group operators for %d resolutions built in %.0f s (genotype-only)',
              sum(!vapply(parts$resolutions, `[[`, NA, 'single_variant')),
              proc.time()[['elapsed']] - started))
  nodewise <- new_nodewise_cache(geno$xcen)
  file <- nodewise_cache_file(cfg)
  have <- load_nodewise(nodewise, file)
  if (have) say('Nodewise cache: ', have, ' variants already computed (', pretty_path(file), ')')
  list(cfg = cfg, ranks = ranks, q = q, geno = geno, truth = truth, parts = parts,
       operators = operators, nodewise = nodewise, nodewise_file = file)
}

# ---------------------------------------------------------------------------
# One case
# ---------------------------------------------------------------------------

analyse_case <- function(ctx, case, replicate) {
  geno <- ctx$geno
  y <- case_phenotype(case, replicate, ctx$cfg, geno, ctx$truth)
  fit <- variant_tests(y, geno$xcen, ctx$nodewise)
  residual <- y - as.numeric(geno$Xc %*% fit$beta)
  causal <- if (case == 'null') NULL else ctx$truth$causal
  tracks <- c('simes', track_id(ctx$ranks))

  rows <- list(); detail <- list()
  for (res in ctx$parts$resolutions) {
    if (res$single_variant) {
      p <- matrix(fit$p, length(fit$p), length(tracks), dimnames = list(NULL, tracks))
      used <- matrix(1L, length(fit$p), length(ctx$ranks), dimnames = list(NULL, ctx$ranks))
    } else {
      gt <- group_tests(ctx$operators[[res$name]], fit$beta, residual, fit$sigma^2, ctx$ranks)
      p <- cbind(simes = simes(fit$p, res$assignment, res$sizes), gt$p)
      colnames(p) <- tracks
      used <- gt$rank
    }
    selected <- lapply(tracks, function(t) bh(p[, t], ctx$q))
    names(selected) <- tracks
    for (t in tracks) {
      iv <- group_intervals(selected[[t]], res$assignment, geno)
      rank <- sub('^rank_', '', t)
      rows[[length(rows) + 1]] <- cbind(
        data.frame(case = case, replicate = replicate, resolution = res$name, track = t,
                   groups = length(res$sizes),
                   mean_rank = if (t == 'simes') NA_real_ else mean(used[, rank])),
        score_intervals(iv, causal, ctx$truth$n_loci))
    }
    detail[[res$name]] <- list(p = p, selected = selected)
  }
  list(rows = do.call(rbind, rows),
       detail = list(sigma = fit$sigma, tested = fit$S, p_variant = fit$p, resolutions = detail))
}

# ---------------------------------------------------------------------------
# Design freeze, replicate loop, status
# ---------------------------------------------------------------------------

design_record <- function(cfg, ranks, q) {
  data_files <- c(file.path(cfg$data_dir, 'data/genotypes', sprintf('example_chr%d.bed', cfg$chromosomes)),
                  file.path(cfg$data_dir, 'data/phenotypes', c('causal.txt', 'phenotypes.tab')))
  if (cfg$partition == 'published')
    data_files <- c(data_files, file.path(cfg$data_dir, sprintf(cfg$partition_files, cfg$chromosomes)))
  list(experiment = cfg$experiment, partition = cfg$partition,
       resolutions = resolution_names(cfg), ranks = ranks, q = q,
       phenotype = cfg$phenotype,
       code_md5 = as.list(tools::md5sum(CODE_FILES)),
       data_md5 = as.list(tools::md5sum(data_files)))
}

as_json <- function(x) jsonlite::toJSON(x, auto_unbox = TRUE, digits = NA, pretty = TRUE)

ensure_design <- function(output, cfg, ranks, q) {
  path <- file.path(output, 'design.json')
  current <- design_record(cfg, ranks, q)
  if (file.exists(path)) {
    saved <- jsonlite::read_json(path, simplifyVector = TRUE)
    saved$created <- NULL
    keys <- union(names(saved), names(current))
    differ <- keys[!vapply(keys, function(k) identical(as.character(as_json(saved[[k]])),
                                                       as.character(as_json(current[[k]]))), logical(1))]
    if (length(differ))
      stop('The run directory ', output, ' was created with a different design (',
           paste(differ, collapse = ', '), '). Use the original settings or a new --out directory.')
    return(invisible(saved))
  }
  if (length(setdiff(list.files(output), 'run.log')))
    stop('The directory ', output, ' is not empty and has no design.json. Choose a new --out directory.')
  current$created <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
  write_lines_atomic(as_json(current), path)
  invisible(current)
}

fit_file <- function(output, replicate)
  file.path(output, 'fits', if (replicate == 0) 'shipped.rds' else sprintf('rep%04d.rds', replicate))

run_cases <- function(ctx, output, replicates, include_null = TRUE) {
  dir.create(file.path(output, 'fits'), recursive = TRUE, showWarnings = FALSE)
  skipped <- 0L
  for (replicate in replicates) {
    path <- fit_file(output, replicate)
    cases <- if (replicate == 0) 'shipped' else c('signal', if (include_null) 'null')
    if (file.exists(path) && all(cases %in% names(readRDS(path)))) {
      say(sprintf('%-8s already complete', label_replicate(replicate))); next
    }
    lock <- paste0(path, '.lock')
    if (!acquire_lock(lock)) {
      say(sprintf('%-8s skipped: in use by %s', label_replicate(replicate), describe_lock(lock)))
      skipped <- skipped + 1L; next
    }
    started <- proc.time()[['elapsed']]
    tryCatch({
      saved <- if (file.exists(path)) readRDS(path) else list()
      for (case in setdiff(cases, names(saved))) {
        saved[[case]] <- analyse_case(ctx, case, replicate)
        save_rds_atomic(saved, path)
      }
      save_nodewise(ctx$nodewise, ctx$nodewise_file)
    }, finally = release_lock(lock))
    say(sprintf('%-8s %s  [%.0f s]', label_replicate(replicate),
                progress_note(saved, ctx), proc.time()[['elapsed']] - started))
  }
  if (skipped) say(skipped, ' replicate(s) were being processed by another run.R and were skipped.')
  invisible(skipped)
}

label_replicate <- function(replicate) if (replicate == 0) 'shipped' else sprintf('rep %04d', replicate)

# e.g. "signal: 105 tested, rank_effective finds 15/20 loci at res2; null: 0 discoveries"
progress_note <- function(saved, ctx) {
  paste(vapply(names(saved), function(case) {
    r <- saved[[case]]$rows
    tested <- length(saved[[case]]$detail$tested)
    if (case == 'null') sprintf('null: %d tested, max discoveries %d', tested, max(r$discoveries))
    else sprintf('%s: %d tested, loci found %d-%d of %d', case, tested,
                 min(r$loci_found), max(r$loci_found), ctx$truth$n_loci)
  }, ''), collapse = '; ')
}

collect_records <- function(output) {
  files <- list.files(file.path(output, 'fits'), pattern = '^(shipped|rep[0-9]+)[.]rds$', full.names = TRUE)
  if (!length(files)) return(NULL)
  rows <- do.call(rbind, lapply(files, function(f) do.call(rbind, lapply(readRDS(f), `[[`, 'rows'))))
  rows[order(rows$replicate, match(rows$case, c('shipped', 'signal', 'null'))), ]
}

show_status <- function(output) {
  cat('Run directory:', pretty_path(output), '\n')
  records <- if (dir.exists(output)) collect_records(output)
  if (is.null(records)) { cat('  no fitted cases yet\n'); return(invisible()) }
  design <- jsonlite::read_json(file.path(output, 'design.json'), simplifyVector = TRUE)
  cat('Experiment:', design$experiment, '  ranks:', paste(design$ranks, collapse = ', '), '  q =', design$q, '\n')
  done <- unique(records[, c('case', 'replicate')])
  for (case in unique(done$case)) {
    ids <- sort(done$replicate[done$case == case])
    cat(sprintf('  %-8s %3d case(s)%s\n', case, length(ids),
                if (case == 'shipped') '' else paste0('  replicates ', format_ids(ids))))
  }
  locks <- list.files(file.path(output, 'fits'), pattern = '[.]lock$', include.dirs = TRUE)
  if (length(locks)) cat('In progress:', paste(sub('[.]rds[.]lock$', '', locks), collapse = ' '), '\n')
}

format_ids <- function(ids) {
  if (length(ids) > 1 && all(diff(ids) == 1)) sprintf('%d:%d', min(ids), max(ids))
  else paste(ids, collapse = ',')
}
