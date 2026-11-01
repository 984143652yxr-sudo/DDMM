# Download the public AD candidate data; optional argument: 21,22 (default all).
script <- gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
setwd(dirname(normalizePath(script)))
args <- commandArgs(TRUE)
chr <- if (length(args)) as.integer(strsplit(args[1], ',')[[1]]) else 1:22
stopifnot(length(chr) > 0, all(chr %in% 1:22))
commit <- '7125b66454137aea19c0a1b9ce2697aaa3f15689'
base <- paste0('https://raw.githubusercontent.com/GuJQ5/FWER_Knockoff_filter/', commit, '/Section%205/')
destination <- 'data/realdata_section5'
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
files <- c('AD_Zscores_candidate.csv', 'Candidate_variants_info_other.csv', sprintf('Cor_%d.csv', unique(chr)))
for (name in files) {
  path <- file.path(destination, name)
  if (file.exists(path) && file.info(path)$size > 0) next
  tmp <- tempfile(tmpdir = destination)
  status <- tryCatch(download.file(paste0(base, name), tmp, mode = 'wb', quiet = TRUE),
                     error = function(e) { unlink(tmp); stop(conditionMessage(e)) })
  if (status != 0 || file.info(tmp)$size == 0) { unlink(tmp); stop('Download failed: ', name) }
  if (!file.rename(tmp, path)) stop('Could not save ', path)
  cat('Downloaded', name, '\n')
}
