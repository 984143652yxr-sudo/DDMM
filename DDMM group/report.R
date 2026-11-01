# Tables and figures from a run directory (reads only; refits nothing).
#
#   Rscript report.R knockoffgwas                         results/knockoffgwas -> reports/knockoffgwas
#   Rscript report.R knockoffgwas --input results/my_run --output reports/my_run
#
# Writes records.csv; for replicates replicates_summary.csv, replicates_table.tex,
# replicates.png/.pdf (power, FDR, global-null rejection by resolution) and
# replicates_contrasts.csv (effective rank against the other tracks, paired
# over the same phenotypes); for the
# shipped phenotype shipped.csv, shipped_table.tex, shipped.png/.pdf and
# shipped_discoveries.png/.pdf (discovery map: the effective-rank test against
# the published knockoff result); and report_source.txt.

script <- normalizePath(gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE))
usage <- function() {
  lines <- readLines(script)
  cat(sub('^# ?', '', lines[seq_len(which(!startsWith(lines, '#'))[1] - 1)]), sep = '\n')
}

options(warn = 1)
caller_directory <- getwd()
root <- dirname(script)
setwd(root)
for (file in c('utils.R', 'ddmm.R', 'scoring.R', 'pipeline.R', 'summarize.R', 'plots.R'))
  source(file.path('R', file))

main <- function() {
  command <- parse_command_line(commandArgs(TRUE), flags = 'help', options = c('input', 'output'))
  opt <- command$options
  if (isTRUE(opt$help) || !length(command$positional)) { usage(); return(invisible()) }
  require_packages(c('jsonlite', 'data.table'))
  x <- command$positional[1]
  cfg <- load_experiment(if (grepl('\\.json$', x)) absolute_path(x, caller_directory) else x)
  input <- if (is.null(opt$input)) absolute_path(cfg$output, root) else absolute_path(opt$input, caller_directory)
  if (!file.exists(file.path(input, 'design.json'))) stop('Not a run directory: ', input)
  report_dir <- if (is.null(opt$output)) default_report_dir(input, root) else absolute_path(opt$output, caller_directory)
  report <- write_report(input, report_dir, cfg)
  cat('Read: ', pretty_path(input), '\nWrote:', pretty_path(report_dir), '/\n')
  cat(paste0('       ', report$files, collapse = '\n'), '\n')
}

tryCatch(main(), error = function(e) { message('Error: ', conditionMessage(e)); quit(status = 1) })
