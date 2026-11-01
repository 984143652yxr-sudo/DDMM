# Summaries over replicates and the shipped-phenotype table.

TRACK_LABELS <- function(tracks, comparator = 'knockoff') {
  vapply(tracks, function(t) {
    if (t == 'simes') 'Simes' else if (t == 'knockoff') comparator
    else rank_label(sub('^rank_', '', t))
  }, '')
}

ordered_tracks <- function(tracks) {
  fixed <- c('rank_effective', 'rank_full', 'simes', 'knockoff')
  ints <- sort(as.integer(sub('^rank_', '', grep('^rank_[0-9]+$', tracks, value = TRUE))))
  c(intersect(c('rank_effective', paste0('rank_', ints), 'rank_full', 'simes', 'knockoff'), tracks),
    setdiff(tracks, c(fixed, paste0('rank_', ints))))
}

cell_seed <- function(...) {
  codes <- utf8ToInt(paste(..., sep = '|'))
  as.integer((sum(codes * seq_along(codes)) * 7919 + 106070) %% 2147483647)
}

bootstrap_ci <- function(x, seed, B = 3000) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(c(NA_real_, NA_real_))
  set.seed(seed)
  unname(stats::quantile(replicate(B, mean(sample(x, replace = TRUE))), c(.025, .975)))
}

# One row per case (signal, null) x resolution x track.
#   power, FDR = mean locus power and mean FDP over replicates, with 95%
#   percentile bootstrap intervals; any_discovery = P(at least one discovery),
#   with a Clopper-Pearson interval (in the null case this equals the FDR).
summarize_replicates <- function(records, resolutions) {
  records <- records[records$case %in% c('signal', 'null'), ]
  if (!nrow(records)) return(NULL)
  out <- list()
  for (case in intersect(c('signal', 'null'), records$case))
    for (res in intersect(resolutions, records$resolution))
      for (track in ordered_tracks(unique(records$track))) {
        z <- records[records$case == case & records$resolution == res & records$track == track, ]
        if (!nrow(z)) next
        any_hit <- sum(z$discoveries > 0)
        ci_any <- stats::binom.test(any_hit, nrow(z))$conf.int
        pw <- if (case == 'signal') bootstrap_ci(z$power, cell_seed(case, res, track, 'power')) else c(NA, NA)
        fd <- bootstrap_ci(z$FDP, cell_seed(case, res, track, 'FDP'))
        out[[length(out) + 1]] <- data.frame(
          case = case, resolution = res, track = track, groups = z$groups[1],
          mean_rank = mean(z$mean_rank), replicates = nrow(z),
          discoveries = mean(z$discoveries),
          power = if (case == 'signal') mean(z$power) else NA_real_,
          power_lo = pw[1], power_hi = pw[2],
          FDR = mean(z$FDP), FDR_lo = fd[1], FDR_hi = fd[2],
          any_discovery = any_hit / nrow(z), any_lo = ci_any[1], any_hi = ci_any[2])
      }
  do.call(rbind, out)
}

# Paired comparison of the tracks on the same phenotypes. For every resolution
# and every other track: mean difference in locus power and in FDP against the
# reference track (the effective rank when it was run), the standard error of
# that paired difference over replicates, and in how many replicates the
# reference found more, as many, or fewer loci.
paired_contrasts <- function(records, resolutions, reference = 'rank_effective') {
  sig <- records[records$case == 'signal', ]
  if (!nrow(sig) || !(reference %in% sig$track)) return(NULL)
  out <- list()
  for (res in intersect(resolutions, sig$resolution))
    for (t in setdiff(ordered_tracks(unique(sig$track)), reference)) {
      a <- sig[sig$resolution == res & sig$track == reference, c('replicate', 'power', 'FDP')]
      b <- sig[sig$resolution == res & sig$track == t, c('replicate', 'power', 'FDP')]
      m <- merge(a, b, by = 'replicate', suffixes = c('_ref', '_other'))
      if (nrow(m) < 2) next
      dp <- m$power_ref - m$power_other; dq <- m$FDP_ref - m$FDP_other
      out[[length(out) + 1]] <- data.frame(
        resolution = res, reference = reference, versus = t, replicates = nrow(m),
        power_diff = mean(dp), power_diff_se = stats::sd(dp) / sqrt(nrow(m)),
        FDR_diff = mean(dq), FDR_diff_se = stats::sd(dq) / sqrt(nrow(m)),
        reference_ahead = sum(dp > 0), tied = sum(dp == 0), reference_behind = sum(dp < 0))
    }
  do.call(rbind, out)
}

shipped_table <- function(records, comparator_rows, resolutions) {
  s <- rbind(records[records$case == 'shipped', ], comparator_rows)
  if (!nrow(s)) return(NULL)
  s <- s[order(match(s$resolution, resolutions), match(s$track, ordered_tracks(unique(s$track)))), ]
  s[, c('resolution', 'track', 'groups', 'mean_rank', 'discoveries', 'TP', 'FDP',
        'loci_found', 'power', 'width_kb')]
}

# "power / FDR" grid: resolutions x tracks
wide_cells <- function(tab, value = function(r) sprintf('%.3f / %.3f', r$power, r$FDR), resolutions) {
  tracks <- ordered_tracks(unique(tab$track))
  view <- data.frame(resolution = intersect(resolutions, tab$resolution))
  for (t in tracks) view[[t]] <- vapply(view$resolution, function(res) {
    r <- tab[tab$resolution == res & tab$track == t, ]
    if (!nrow(r)) '' else value(r)
  }, '')
  view
}

latex_grid <- function(view, caption, label, labels) {
  tracks <- setdiff(names(view), 'resolution')
  c('\\begin{table}[htbp]\\centering\\small',
    sprintf('\\caption{%s}', caption), sprintf('\\label{%s}', label),
    paste0('\\begin{tabular}{l', strrep('c', length(tracks)), '}\\toprule'),
    paste0('Resolution & ', paste(labels[tracks], collapse = ' & '), '\\\\\\midrule'),
    vapply(seq_len(nrow(view)), function(i)
      paste0(view$resolution[i], ' & ', paste(gsub(' / ', '/', unlist(view[i, tracks])), collapse = ' & '), '\\\\'), ''),
    '\\bottomrule\\end{tabular}', '\\end{table}')
}

# ---------------------------------------------------------------------------
# Every report file for one run directory.
# A report directory is (re)written only when new, empty, or an earlier
# report.R output for the same input, so older reports are never overwritten.
# ---------------------------------------------------------------------------

REPORT_MARKER <- 'Written by report.R (ddmm_group package).'

check_report_directory <- function(report_dir, input) {
  if (!dir.exists(report_dir) || !length(list.files(report_dir))) return(invisible(TRUE))
  source_file <- file.path(report_dir, 'report_source.txt')
  previous <- if (file.exists(source_file)) readLines(source_file, warn = FALSE) else character()
  if (!(REPORT_MARKER %in% previous && paste('Input:', pretty_path(input)) %in% previous))
    stop('The report directory ', report_dir, ' already holds other files. Choose a new --output directory.')
  invisible(TRUE)
}

write_report <- function(output, report_dir, cfg) {
  check_report_directory(report_dir, output)
  records <- collect_records(output)
  if (is.null(records)) stop('No fitted cases in ', output, '. Run run.R first.')
  layout <- readRDS(file.path(output, 'layout.rds'))
  design <- jsonlite::read_json(file.path(output, 'design.json'), simplifyVector = TRUE)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  res <- layout$resolutions
  all_tracks <- c(unique(records$track), 'knockoff')
  labels <- TRACK_LABELS(all_tracks, paste(cfg$comparator$name, '(published)'))
  written <- character()
  file_out <- function(name) { f <- file.path(report_dir, name); written <<- c(written, name); f }

  write_csv_atomic(records, file_out('records.csv'))

  summary <- summarize_replicates(records, res)
  if (!is.null(summary)) {
    write_csv_atomic(summary, file_out('replicates_summary.csv'))
    sig <- summary[summary$case == 'signal', ]
    if (nrow(sig)) write_lines_atomic(latex_grid(wide_cells(sig, resolutions = res),
      sprintf('%s: locus power / mean FDP over %d replicate phenotypes, BH at $q=%g$ within each resolution.',
              cfg$title, max(sig$replicates), design$q),
      paste0('tab:', cfg$experiment, '-replicates'), labels), file_out('replicates_table.tex'))
    draw_replicates(summary, res, design$q, labels, file_out('replicates.pdf'), file_out('replicates.png'))
    contrasts <- paired_contrasts(records, res)
    if (!is.null(contrasts)) write_csv_atomic(contrasts, file_out('replicates_contrasts.csv'))
  }

  if (any(records$case == 'shipped')) {
    shipped <- shipped_table(records, comparator_scores(cfg, layout), res)
    write_csv_atomic(shipped, file_out('shipped.csv'))
    write_lines_atomic(latex_grid(
      wide_cells(transform(shipped, FDR = FDP), resolutions = res),
      sprintf('%s, shipped phenotype: locus power / FDP at target FDR %g. %s is the published result.',
              cfg$title, design$q, cfg$comparator$name),
      paste0('tab:', cfg$experiment, '-shipped'), labels), file_out('shipped_table.tex'))
    draw_shipped(transform(shipped, FDR = FDP), res, design$q, labels,
                 file_out('shipped.pdf'), file_out('shipped.png'))
    # Discovery map: the effective-rank test (or the first rank run) against the
    # published knockoff result.
    fit <- readRDS(fit_file(output, 0))$shipped$detail$resolutions
    tracks <- track_id(if ('effective' %in% design$ranks) 'effective' else design$ranks[1])
    intervals <- lapply(tracks, function(t) {
      out <- lapply(res, function(r) group_intervals(fit[[r]]$selected[[t]], layout$assignments[[r]], layout))
      names(out) <- res; out
    })
    names(intervals) <- tracks
    intervals$knockoff <- comparator_intervals(cfg, layout)
    draw_discovery_map(layout, intervals, c(tracks, 'knockoff'), labels,
                       file_out('shipped_discoveries.pdf'), file_out('shipped_discoveries.png'))
  }

  counts <- unique(records[, c('case', 'replicate')])
  write_lines_atomic(c(REPORT_MARKER,
    paste('Report written:', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
    paste('Input:', pretty_path(output)),
    paste('Experiment:', cfg$experiment, '-', cfg$title),
    paste('Ranks:', paste(design$ranks, collapse = ', '), '  target FDR:', design$q),
    sprintf('Cases: shipped %d, signal %d, null %d', sum(counts$case == 'shipped'),
            sum(counts$case == 'signal'), sum(counts$case == 'null'))),
    file_out('report_source.txt'))
  list(summary = summary, shipped = if (any(records$case == 'shipped')) shipped, files = written,
       contrasts = if (!is.null(summary)) paired_contrasts(records, res),
       labels = labels, resolutions = res)
}
