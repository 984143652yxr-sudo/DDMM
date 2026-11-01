# Summaries, LaTeX table and figures from fitted results.
# Reads only; never changes the run directory or CSV it reads.
#
#   Rscript report.R ancestry                       report results/ancestry -> reports/ancestry
#   Rscript report.R kinship                        report results/kinship  -> reports/kinship
#   Rscript report.R ancestry --input results/ancestry_independent_hmm
#   Rscript report.R kinship  --input results/independent_trial_all/replicates.csv \
#                             --output reports/kinship_historical
#
# Options
#   --input PATH    run directory (reads fits/**/fit.rds) or a replicates CSV
#                   (default: the experiment's run directory)
#   --output DIR    report directory (default: reports/<input folder name>;
#                   for a CSV: reports/<folder of the CSV>_<experiment>)
#
# Output files: replicates.csv, summary.csv, fit_counts.csv, failures.csv,
# <figure>_table.tex, <figure>.png/.pdf (power and mean FDP),
# <figure>_null.png/.pdf (null case), report_source.txt.

script <- normalizePath(gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE))
usage <- function() {
  lines <- readLines(script)
  cat(sub('^# ?', '', lines[seq_len(which(!startsWith(lines, '#'))[1] - 1)]), sep = '\n')
}

options(warn = 1)
caller_directory <- getwd()
root <- dirname(script)
setwd(root)
for (file in c('utils.R', 'config.R', 'methods.R', 'summarize.R', 'plots.R'))
  source(file.path('R', file))

main <- function() {
  command <- parse_command_line(commandArgs(TRUE), flags = 'help', options = c('input', 'output'))
  opt <- command$options
  if (isTRUE(opt$help) || !length(command$positional)) { usage(); return(invisible()) }
  require_packages('jsonlite')

  cfg <- load_experiment(command$positional[1])
  input <- absolute_path(opt$input %||% cfg$output, caller_directory)
  default_output <- if (dir.exists(input)) basename(input) else
    paste0(basename(dirname(input)), '_', cfg$experiment)
  report_dir <- absolute_path(opt$output %||% file.path(root, 'reports', default_output), caller_directory)
  if (identical(report_dir, if (dir.exists(input)) input else dirname(input)))
    stop('Choose a report directory different from the input.')

  records <- read_records(input, cfg)
  report <- write_report(records, cfg, report_dir, input)

  cat('Read:   ', pretty_path(input), '\n')
  cat('Wrote:  ', pretty_path(report_dir), '/\n', sep = '')
  cat(paste0('          ', basename(report$files), collapse = '\n'), '\n\n')
  counts <- replicate_counts(records)
  cat('Replicates per panel:', paste0(unique(counts$panel), ' ',
      tapply(counts$replicates, counts$panel, max)[unique(counts$panel)], collapse = ', '), '\n')
  cat('Power / mean FDP (valid fits):\n')
  print(wide_view(report$summary, cfg), row.names = FALSE)
}

tryCatch(main(), error = function(e) {
  message('Error: ', conditionMessage(e))
  quit(status = 1)
})
