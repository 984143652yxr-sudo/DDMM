# From fitted checkpoints (or a results CSV) to scored replicates and to the
# summary table used by the figures and LaTeX tables.

RESULT_COLUMNS <- c('experiment', 'panel', 'a', 'replicate', 'method',
                    'valid', 'R', 'TP', 'FP', 'power', 'FDP', 'error')

# Score every fit.rds of a run directory: one row per case x method.
collect_run <- function(output, cfg) {
  files <- list.files(file.path(output, 'fits'), pattern = '^fit[.]rds$',
                      recursive = TRUE, full.names = TRUE)
  if (!length(files)) return(NULL)
  rows <- lapply(files, function(path) {
    saved <- readRDS(path)
    if (!length(saved$fits)) return(NULL)
    truth <- which(saved$generated$beta != 0)
    do.call(rbind, lapply(names(saved$fits), function(method)
      cbind(data.frame(experiment = cfg$experiment, panel = saved$panel, a = saved$a,
                       replicate = saved$replicate, method = method),
            score_fit(saved$fits[[method]], truth))))
  })
  records <- do.call(rbind, rows)
  if (is.null(records)) return(NULL)
  tidy_records(records, cfg)
}

# Read a run directory or a CSV (including the historical
# results/independent_trial_all/replicates.csv, which holds both experiments).
read_records <- function(input, cfg) {
  if (dir.exists(input)) {
    records <- collect_run(input, cfg)
    if (is.null(records) && file.exists(file.path(input, 'replicates.csv')))
      records <- utils::read.csv(file.path(input, 'replicates.csv'), stringsAsFactors = FALSE)
    if (is.null(records)) {
      prepared <- length(list.files(input, pattern = '^phenotype[.]rds$', recursive = TRUE))
      stop('No fitted results in ', input, ' (', prepared, ' prepared phenotypes). ',
           'Run the methods first: Rscript run.R ', cfg$experiment, ' --replicates 1:50 --out ', input)
    }
  } else if (file.exists(input)) {
    records <- utils::read.csv(input, stringsAsFactors = FALSE)
  } else stop('Input not found: ', input)
  if ('level' %in% names(records)) records <- records[records$level == 'SNP', ]
  if ('experiment' %in% names(records)) records <- records[records$experiment == cfg$experiment, ]
  records <- records[records$panel %in% cfg$panels & records$method %in% cfg$methods, ]
  if (!nrow(records)) stop('The input has no ', cfg$experiment, ' results for panels ',
                           paste(cfg$panels, collapse = ','), ': ', input)
  records$experiment <- cfg$experiment
  tidy_records(records, cfg)
}

tidy_records <- function(records, cfg) {
  records$error[is.na(records$error)] <- ''
  records$valid <- as.logical(records$valid)
  records <- records[, RESULT_COLUMNS]
  key <- records[, c('panel', 'replicate', 'a', 'method')]
  if (anyDuplicated(key))
    stop('Duplicate panel/replicate/a/method rows: do not concatenate cumulative CSVs.')
  records[order(match(records$panel, cfg$panels), records$replicate, -records$a,
                match(records$method, cfg$methods)), ]
}

# ---------------------------------------------------------------------------
# Summary over replicates: valid fits only.
#   power   mean of TP/20 (signal cases)
#   FDR     mean FDP; in the null case this is P(any selection) (= FWER)
#   *_lo/_hi  95% percentile bootstrap interval of the mean (3000 resamples)
#             for signal cases; exact Clopper-Pearson interval for P(R > 0)
#             in the null case.
# Each cell has its own bootstrap seed, so intervals do not change when
# other cells (methods, panels, replicates elsewhere) are added.
# ---------------------------------------------------------------------------

cell_seed <- function(...) {
  codes <- utf8ToInt(paste(..., sep = '|'))
  as.integer((sum(codes * seq_along(codes)) * 7919 + 106060) %% 2147483647)
}

bootstrap_mean_ci <- function(x, seed, B = 3000) {
  if (length(x) < 2) return(c(NA_real_, NA_real_))
  set.seed(seed)
  unname(stats::quantile(replicate(B, mean(sample(x, replace = TRUE))), c(.025, .975)))
}

summarize_records <- function(records, cfg) {
  methods <- intersect(cfg$methods, unique(records$method))
  panels <- intersect(cfg$panels, unique(records$panel))
  out <- list()
  for (panel in panels) for (a in sort(unique(records$a[records$panel == panel]), decreasing = TRUE))
    for (method in methods) {
      all <- records[records$panel == panel & records$a == a & records$method == method, ]
      if (!nrow(all)) next
      z <- all[all$valid, ]
      if (a > 0) {
        power_ci <- bootstrap_mean_ci(z$power, cell_seed(panel, a, method, 'power'))
        fdr_ci <- bootstrap_mean_ci(z$FDP, cell_seed(panel, a, method, 'FDP'))
      } else {
        power_ci <- c(NA_real_, NA_real_)
        fdr_ci <- if (nrow(z)) stats::binom.test(sum(z$R > 0), nrow(z))$conf.int[1:2] else c(NA, NA)
      }
      out[[length(out) + 1]] <- data.frame(
        panel = panel, a = a, method = method,
        attempted = nrow(all), valid = nrow(z),
        mean_selected = if (nrow(z)) mean(z$R) else NA_real_,
        power = if (a > 0 && nrow(z)) mean(z$power) else NA_real_,
        power_lo = power_ci[1], power_hi = power_ci[2],
        FDR = if (nrow(z)) mean(z$FDP) else NA_real_,
        FDR_lo = fdr_ci[1], FDR_hi = fdr_ci[2])
    }
  summary <- do.call(rbind, out)
  summary[order(match(summary$panel, panels), summary$a == 0, summary$a,
                match(summary$method, methods)), ]
}

# Compact "power / FDP" view for the console and the LaTeX table.
cell_text <- function(power, fdr, valid) ifelse(valid == 0 | is.na(fdr), '--',
  ifelse(is.na(power), sprintf('-- / %.3f', fdr), sprintf('%.3f / %.3f', power, fdr)))

wide_view <- function(summary, cfg) {
  methods <- intersect(cfg$methods, summary$method)
  rows <- unique(summary[, c('panel', 'a')])
  view <- rows
  for (m in methods) {
    s <- summary[summary$method == m, ]
    s <- s[match(paste(rows$panel, rows$a), paste(s$panel, s$a)), ]
    view[[m]] <- cell_text(s$power, s$FDR, s$valid)
  }
  view
}

replicate_counts <- function(records) {
  tab <- aggregate(replicate ~ panel + a, records, function(r) length(unique(r)))
  names(tab)[3] <- 'replicates'
  tab
}

# ---------------------------------------------------------------------------
# LaTeX table: methods x panels (several panels) or amplitudes x methods
# (one panel). Cells are power / mean FDP among valid fits.
# ---------------------------------------------------------------------------

latex_table <- function(summary, cfg, records) {
  labels <- unlist(cfg$table_labels)
  methods <- intersect(cfg$methods, summary$method)
  panels <- intersect(cfg$panels, summary$panel)
  signal <- summary[summary$a > 0, ]
  n_rep <- max(replicate_counts(records)$replicates)
  cell <- function(s) if (!nrow(s) || s$valid == 0) '--' else sprintf('%.3f/%.3f', s$power, s$FDR)
  lines <- character()
  if (length(panels) > 1) {
    for (a in sort(unique(signal$a))) {
      lines <- c(lines, '\\begin{table}[htbp]\\centering\\small',
        sprintf('\\caption{%s, $a=%g$: exact-SNP power / mean FDP among valid fits (%d replicates per panel). Valid/attempted counts are in fit\\_counts.csv.}', cfg$table_caption, a, n_rep),
        sprintf('\\label{tab:%s-a%g}', cfg$figure_name, a),
        paste0('\\begin{tabular}{l', strrep('r', length(panels)), '}\\toprule'),
        paste0('Method & ', paste(sub('[0-9]+$', '', panels), collapse = ' & '), '\\\\\\midrule'))
      for (m in methods) lines <- c(lines, paste0(labels[[m]], ' & ', paste(vapply(panels, function(p)
        cell(signal[signal$panel == p & signal$a == a & signal$method == m, ]), ''), collapse = ' & '), '\\\\'))
      lines <- c(lines, '\\bottomrule\\end{tabular}', '\\end{table}', '')
    }
  } else {
    lines <- c('\\begin{table}[htbp]\\centering\\small',
      sprintf('\\caption{%s: exact-SNP power / mean FDP among valid fits (%d replicates).}', cfg$table_caption, n_rep),
      sprintf('\\label{tab:%s}', cfg$figure_name),
      paste0('\\begin{tabular}{r', strrep('l', length(methods)), '}\\toprule'),
      paste0('$a$ & ', paste(labels[methods], collapse = ' & '), '\\\\\\midrule'))
    for (a in sort(unique(signal$a))) lines <- c(lines, paste0(a, ' & ', paste(vapply(methods, function(m)
      cell(signal[signal$a == a & signal$method == m, ]), ''), collapse = ' & '), '\\\\'))
    lines <- c(lines, '\\bottomrule\\end{tabular}\\end{table}')
  }
  lines
}

# ---------------------------------------------------------------------------
# Write every report file into report_dir.
# A report directory is only (re)written when it is new, empty, or an earlier
# report.R output for the same input; older reports are never overwritten.
# ---------------------------------------------------------------------------

REPORT_MARKER <- 'Written by report.R (semi-synthetic-1000G package).'

check_report_directory <- function(report_dir, input) {
  if (!dir.exists(report_dir) || !length(list.files(report_dir))) return(invisible(TRUE))
  source_file <- file.path(report_dir, 'report_source.txt')
  previous <- if (file.exists(source_file)) readLines(source_file, warn = FALSE) else character()
  same_input <- paste('Input:', pretty_path(input)) %in% previous
  if (!(REPORT_MARKER %in% previous && same_input))
    stop('The report directory ', report_dir, ' already holds other files (an older report ',
         'or a report of a different input). Choose a new --output directory.')
  invisible(TRUE)
}

write_report <- function(records, cfg, report_dir, input) {
  check_report_directory(report_dir, input)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  summary <- summarize_records(records, cfg)
  files <- c(
    replicates = file.path(report_dir, 'replicates.csv'),
    summary = file.path(report_dir, 'summary.csv'),
    fit_counts = file.path(report_dir, 'fit_counts.csv'),
    failures = file.path(report_dir, 'failures.csv'),
    table = file.path(report_dir, paste0(cfg$figure_name, '_table.tex')),
    figure_png = file.path(report_dir, paste0(cfg$figure_name, '.png')),
    figure_pdf = file.path(report_dir, paste0(cfg$figure_name, '.pdf')),
    null_png = file.path(report_dir, paste0(cfg$figure_name, '_null.png')),
    null_pdf = file.path(report_dir, paste0(cfg$figure_name, '_null.pdf')),
    source = file.path(report_dir, 'report_source.txt'))
  write_csv_atomic(records, files[['replicates']])
  write_csv_atomic(summary, files[['summary']])
  write_csv_atomic(summary[, c('panel', 'a', 'method', 'attempted', 'valid')], files[['fit_counts']])
  write_csv_atomic(records[!records$valid, c('panel', 'a', 'replicate', 'method', 'error')], files[['failures']])
  write_lines_atomic(latex_table(summary, cfg, records), files[['table']])
  draw_main_figure(summary, cfg, files[['figure_pdf']], files[['figure_png']])
  has_null <- any(summary$a == 0)
  if (has_null) draw_null_figure(summary, cfg, files[['null_pdf']], files[['null_png']])
  else files <- files[!names(files) %in% c('null_png', 'null_pdf')]
  counts <- replicate_counts(records)
  write_lines_atomic(c(
    REPORT_MARKER,
    paste('Report written:', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
    paste('Input:', pretty_path(input)),
    paste('Experiment:', cfg$experiment, '-', cfg$title),
    paste('Model:', cfg$model),
    paste('Method/case rows:', nrow(records)),
    'Replicates per panel and amplitude:',
    paste0('  ', counts$panel, '  a=', counts$a, '  ', counts$replicates)), files[['source']])
  list(summary = summary, files = files)
}
