# Restore duplicate matrices from the compact panel without recomputing eigenvectors.
read_panel <- function(name) {
  d <- readRDS(file.path('data', name, 'panel.rds'))
  if (!is.null(d$U)) {
    d$Q_none <- t(d$U)
    d$Q_pca <- d$Q_none[-(1:3), , drop = FALSE]
    d$K <- tcrossprod(sweep(d$U, 2, sqrt(d$values), '*'))
  }
  d
}

# External tools require unpacked files. Expand them only for the selected methods.
prepare_external_inputs <- function(cfg) {
  names <- c(if ('Knockoff' %in% cfg$methods) 'hap.haps',
             if ('BOLT' %in% cfg$methods) c('bolt.bed', 'bolt.bim'))
  for (panel in cfg$panels) for (name in names) {
    target <- file.path('data', panel, name)
    if (file.exists(target)) next
    packed <- paste0(target, '.gz')
    if (!file.exists(packed)) stop('Missing input: ', packed)
    tmp <- tempfile(tmpdir = dirname(target))
    input <- gzfile(packed, 'rb'); output <- file(tmp, 'wb')
    tryCatch({
      repeat {
        bytes <- readBin(input, 'raw', n = 1024L * 1024L)
        if (!length(bytes)) break
        writeBin(bytes, output)
      }
    }, error = function(e) { unlink(tmp); stop(conditionMessage(e)) },
       finally = { close(input); close(output) })
    if (!file.rename(tmp, target)) stop('Could not expand ', packed)
  }
  invisible(TRUE)
}
