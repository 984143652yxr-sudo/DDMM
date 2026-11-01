# Run the semi-synthetic 1000 Genomes experiments end to end:
# simulate phenotypes -> fit DDMM, SuSiE, HMM Knockoff+, BOLT-LMM -> report.
#
#   Rscript run.R ancestry --replicates 1:50        replicates 1..50, then report
#   Rscript run.R ancestry --replicates 51:100      add replicates 51..100
#   Rscript run.R kinship  --replicates 1:50
#   Rscript run.R ancestry --status                 progress of the run directory
#
# Options
#   --replicates SPEC   global replicate IDs: 1:50, 51:100, 7, 1,3,5 (1..999)
#   --out DIR           run directory (default: results/ancestry, results/kinship)
#   --panels A,B        only these panels (e.g. --panels AFR300,MIX300)
#   --report DIR        report directory (default: reports/<run directory name>)
#   --no-report         do not write the report at the end
#   --prepare-only      simulate phenotypes and knockoff copies, fit nothing
#   --amplitudes LIST   override the signal amplitudes (needs a new --out)
#   --methods LIST      override the methods, e.g. DDMM,SuSiE (needs a new --out)
#   --status            show progress and exit
#
# Interrupted runs resume where they stopped when the same command is
# repeated. Several run.R processes may share one run directory; each
# replicate is processed by exactly one of them.

script <- normalizePath(gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE))
usage <- function() {
  lines <- readLines(script)
  cat(sub('^# ?', '', lines[seq_len(which(!startsWith(lines, '#'))[1] - 1)]), sep = '\n')
}

options(warn = 1)
caller_directory <- getwd()          # user paths are relative to where you call from
root <- dirname(script)              # all package paths are relative to the package
setwd(root)
for (file in c('utils.R', 'config.R', 'simulate.R', 'methods.R', 'knockoffs.R',
               'pipeline.R', 'summarize.R', 'plots.R', 'cache.R')) source(file.path('R', file))

main <- function() {
  command <- parse_command_line(commandArgs(TRUE),
    flags = c('status', 'no-report', 'prepare-only', 'help'),
    options = c('replicates', 'out', 'panels', 'report', 'amplitudes', 'methods'))
  opt <- command$options
  if (isTRUE(opt$help) || !length(command$positional)) { usage(); return(invisible()) }

  require_packages('jsonlite')
  cfg <- load_experiment(command$positional[1])
  if (!is.null(opt$amplitudes)) cfg$amplitudes <- parse_number_list(opt$amplitudes)
  if (!is.null(opt$methods)) cfg$methods <- parse_name_list(opt$methods)
  output <- absolute_path(opt$out %||% cfg$output, caller_directory)
  panels <- if (is.null(opt$panels)) cfg$panels else parse_name_list(opt$panels)
  if (!all(panels %in% cfg$panels))
    stop('Panels for ', cfg$experiment, ' are ', paste(cfg$panels, collapse = ', '))

  if (isTRUE(opt$status) || is.null(opt$replicates)) {
    show_status(cfg, output)
    if (is.null(opt$replicates) && !isTRUE(opt$status))
      cat('\nTo run replicates: Rscript run.R', command$positional[1], '--replicates 1:50\n')
    return(invisible())
  }

  require_packages(c('glmnet', if ('SuSiE' %in% cfg$methods) 'susieR',
                     if ('Knockoff' %in% cfg$methods) c('knockoff', 'snpStats')))
  require_programs(cfg$methods)
  for (panel in cfg$panels) build_cache(panel, cfg$projection)
  replicates <- parse_replicates(opt$replicates)
  run_experiment(cfg, replicates, output, panels, prepare_only = isTRUE(opt$`prepare-only`))
  if (isTRUE(opt$`prepare-only`)) return(invisible())

  records <- collect_run(output, cfg)
  if (is.null(records)) { say('No fitted cases yet in ', pretty_path(output)); return(invisible()) }
  write_csv_atomic(records, file.path(output, 'replicates.csv'))
  say('Scored ', nrow(records), ' method/case rows -> ', pretty_path(file.path(output, 'replicates.csv')))
  if (!isTRUE(opt$`no-report`)) {
    report_dir <- absolute_path(opt$report %||% file.path(root, 'reports', basename(output)), caller_directory)
    report <- write_report(records, cfg, report_dir, output)
    say('Report -> ', pretty_path(report_dir), '/  (', paste(basename(report$files), collapse = ', '), ')')
    cat('\nPower / mean FDP (valid fits):\n')
    print(wide_view(report$summary, cfg), row.names = FALSE)
  }
}

# Errors are reported as one plain message (no R call stack).
log_line <- function(text) if (!is.null(.log$file))
  cat(format(Sys.time(), '%H:%M:%S'), '  ', text, '\n', sep = '', file = .log$file, append = TRUE)
tryCatch(main(), error = function(e) {
  message('Error: ', conditionMessage(e))
  log_line(paste('Error:', conditionMessage(e)))
  quit(status = 1)
}, interrupt = function(e) {
  message('Interrupted. Finished methods are saved; repeat the same command to resume.')
  log_line('Interrupted by the user.')
  quit(status = 130)
})
