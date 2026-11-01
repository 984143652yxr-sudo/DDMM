script <- gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
setwd(dirname(normalizePath(script)))
args <- commandArgs(TRUE); kind <- if (length(args)) args[1] else 'admixture'
stopifnot(kind %in% c('admixture', 'family'))
out <- file.path('results', kind)
x <- read.csv(file.path(out, 'replicates.csv'))
summary <- aggregate(cbind(FDP, power, MSPE, L1) ~ method, x, mean)
write.csv(summary, file.path(out, 'summary.csv'), row.names = FALSE)
pdf(file.path(out, 'comparison.pdf'), width = 11, height = 7)
par(mfrow = c(2, 2), mar = c(8, 4, 2, 1))
for (metric in c('FDP', 'power', 'MSPE', 'L1')) {
  vals <- summary[[metric]]
  barplot(vals, names.arg = summary$method, las = 2, col = '#9bbcd1', cex.names = 0.7,
          ylab = if (metric == 'FDP') 'Mean FDP' else metric, main = kind)
  if (metric == 'FDP') abline(h = 0.1, lty = 2, col = '#aa5544')
}
dev.off()
