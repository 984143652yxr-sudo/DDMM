# DDMM group-level tests on the KnockoffGWAS toy data and its published
# partitions: per-variant debiased Lasso, Simes, and the debiased group
# chi-squared at a chosen rank, with BH within each resolution, scored against
# the 20 causal loci.
#
#   Rscript run.R knockoffgwas --shipped                      published phenotype vs published KnockoffGWAS
#   Rscript run.R knockoffgwas --replicates 1:20              20 redrawn phenotypes + 20 global nulls
#   Rscript run.R knockoffgwas --shipped --replicates 1:40 --rank effective,3,full
#   Rscript run.R knockoffgwas --replicates 1:20 --rank 3 --out results/knockoffgwas_K3
#   Rscript run.R knockoffgwas --status
#
# Options
#   --rank LIST       rank of the group chi-squared: effective (participation
#                     ratio of the group's LD eigenvalues), a fixed K such as 3,
#                     or full (all |g| directions); a list runs several tracks.
#                     Default: the ranks of an existing run directory, else effective
#   --shipped         analyse the phenotype shipped with the data
#   --replicates IDS  redrawn phenotypes: 1:20, 21:40, 7, 1,3,5
#   --no-null         skip the global-null case of each replicate
#   --out DIR         run directory (default results/<experiment>)
#   --report DIR      report directory (default reports/<run directory name>)
#   --no-report       do not write the report at the end
#   --status          show what the run directory holds and exit
#
# Repeating a command resumes it; several run.R processes may share one run
# directory (each replicate is processed by one of them).

script <- normalizePath(gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE))
usage <- function() {
  lines <- readLines(script)
  cat(sub('^# ?', '', lines[seq_len(which(!startsWith(lines, '#'))[1] - 1)]), sep = '\n')
}

options(warn = 1)
caller_directory <- getwd()
root <- dirname(script)
setwd(root)
for (file in c('utils.R', 'data.R', 'partitions.R', 'ddmm.R', 'scoring.R', 'pipeline.R',
               'summarize.R', 'plots.R')) source(file.path('R', file))

main <- function() {
  command <- parse_command_line(commandArgs(TRUE),
    flags = c('shipped', 'no-null', 'no-report', 'status', 'help'),
    options = c('rank', 'replicates', 'out', 'report'))
  opt <- command$options
  if (isTRUE(opt$help) || !length(command$positional)) { usage(); return(invisible()) }
  require_packages(c('jsonlite', 'glmnet', 'snpStats', 'data.table'))
  cfg <- load_experiment(experiment_argument(command$positional[1]))
  if (cfg$partition == 'adjclust') require_packages('adjclust')
  output <- if (is.null(opt$out)) absolute_path(cfg$output, root) else absolute_path(opt$out, caller_directory)

  replicates <- c(if (isTRUE(opt$shipped)) 0L, if (!is.null(opt$replicates)) parse_replicates(opt$replicates))
  if (isTRUE(opt$status) || !length(replicates)) {
    show_status(output)
    if (!isTRUE(opt$status)) cat('\nTo run: Rscript run.R', command$positional[1], '--shipped --replicates 1:40\n')
    return(invisible())
  }
  # Without --rank: the ranks of an existing run directory, otherwise effective.
  design_file <- file.path(output, 'design.json')
  ranks <- parse_ranks(opt$rank %||% if (file.exists(design_file))
    paste(jsonlite::read_json(design_file, simplifyVector = TRUE)$ranks, collapse = ',') else 'effective')
  q <- cfg$q

  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  log_to(file.path(output, 'run.log'))
  ensure_design(output, cfg, ranks, q)
  say(sprintf('Run %s: %s, ranks %s, q = %g, output %s', cfg$experiment,
              paste(c(if (0 %in% replicates) 'shipped phenotype',
                      if (any(replicates > 0)) paste('replicates', format_ids(replicates[replicates > 0]))),
                    collapse = ' + '),
              paste(ranks, collapse = ','), q, pretty_path(output)))
  include_null <- !isTRUE(opt$`no-null`)
  layout_file <- file.path(output, 'layout.rds')
  pending <- Filter(function(r) !replicate_complete(output, r, include_null), replicates)
  if (length(pending) || !file.exists(layout_file)) {
    ctx <- build_context(cfg, ranks, q)
    if (!file.exists(layout_file)) save_rds_atomic(make_layout(ctx), layout_file)
    run_cases(ctx, output, replicates, include_null = include_null)
  } else say('All requested cases are already in ', pretty_path(output), '; nothing to fit.')

  records <- collect_records(output)
  write_csv_atomic(records, file.path(output, 'records.csv'))
  if (!isTRUE(opt$`no-report`)) {
    report_dir <- if (is.null(opt$report)) default_report_dir(output, root) else absolute_path(opt$report, caller_directory)
    report <- write_report(output, report_dir, cfg)
    say('Report -> ', pretty_path(report_dir), '/  (', paste(report$files, collapse = ', '), ')')
    print_overview(report)
  }
}

# A config file given by path is read relative to where run.R was called from.
experiment_argument <- function(x) if (grepl('\\.json$', x)) absolute_path(x, caller_directory) else x

replicate_complete <- function(output, replicate, include_null) {
  path <- fit_file(output, replicate)
  cases <- if (replicate == 0) 'shipped' else c('signal', if (include_null) 'null')
  file.exists(path) && all(cases %in% names(readRDS(path)))
}

print_overview <- function(report) {
  if (!is.null(report$shipped)) {
    cat('\nShipped phenotype, locus power / FDP:\n')
    print(wide_cells(transform(report$shipped, FDR = FDP), resolutions = report$resolutions), row.names = FALSE)
  }
  if (!is.null(report$summary)) for (case in intersect(c('signal', 'null'), report$summary$case)) {
    s <- report$summary[report$summary$case == case, ]
    cat(sprintf('\n%s, %d replicates, %s:\n', if (case == 'signal') 'Signal' else 'Global null',
                max(s$replicates), if (case == 'signal') 'mean power / mean FDP' else 'P(any discovery)'))
    value <- if (case == 'signal') function(r) sprintf('%.3f / %.3f', r$power, r$FDR)
             else function(r) sprintf('%.2f', r$any_discovery)
    print(wide_cells(s, value, report$resolutions), row.names = FALSE)
  }
  k <- report$contrasts
  if (!is.null(k)) k <- k[k$versus != 'simes', ]
  if (!is.null(k) && nrow(k)) {
    cat('\nEffective rank minus the other ranks, paired over replicates: power difference (SE):\n')
    k$cell <- sprintf('%+.3f (%.3f)', k$power_diff, k$power_diff_se)
    view <- data.frame(resolution = unique(k$resolution))
    for (t in unique(k$versus)) view[[t]] <- k$cell[k$versus == t][match(view$resolution, k$resolution[k$versus == t])]
    print(view, row.names = FALSE)
  }
}

log_line <- function(text) if (!is.null(.log$file))
  cat(format(Sys.time(), '%H:%M:%S'), '  ', text, '\n', sep = '', file = .log$file, append = TRUE)
tryCatch(main(), error = function(e) {
  message('Error: ', conditionMessage(e)); log_line(paste('Error:', conditionMessage(e))); quit(status = 1)
}, interrupt = function(e) {
  message('Interrupted. Finished cases are saved; repeat the same command to resume.')
  log_line('Interrupted by the user.'); quit(status = 130)
})
