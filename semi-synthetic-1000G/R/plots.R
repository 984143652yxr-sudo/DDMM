# Figures (base R graphics, no extra packages). Colours: Okabe-Ito.

METHOD_COLOURS <- c(DDMM = '#0072B2', SuSiE = '#D55E00', Knockoff = '#009E73', BOLT = '#999999')
METHOD_SYMBOLS <- c(DDMM = 16, SuSiE = 17, Knockoff = 18, BOLT = 19)

save_both <- function(draw, pdf_path, png_path, width, height) {
  grDevices::pdf(pdf_path, width = width, height = height); draw(); grDevices::dev.off()
  grDevices::png(png_path, width = width * 300, height = height * 300, res = 300); draw(); grDevices::dev.off()
}

legend_strip <- function(methods, labels) {
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  width <- max(graphics::strwidth(labels[methods], cex = .85)) * 1.15
  graphics::legend('center', legend = labels[methods], col = METHOD_COLOURS[methods],
                   pch = METHOD_SYMBOLS[methods], lty = 1, horiz = TRUE, bty = 'n', cex = .85,
                   text.width = width)
}

panel_label <- function(panel) sub('[0-9]+$', '', panel)

# Dots with 95% intervals at each panel (several panels, one amplitude).
dot_panel <- function(s, metric, panels, methods, ylab, ylim, reference = NULL) {
  graphics::plot(NA, xlim = c(.6, length(panels) + .4), ylim = ylim, xaxt = 'n',
                 xlab = '', ylab = ylab, bty = 'l')
  graphics::axis(1, seq_along(panels), panel_label(panels))
  if (!is.null(reference)) graphics::abline(h = reference, lty = 2, col = 'grey50')
  for (i in seq_along(methods)) {
    z <- s[s$method == methods[i], ]
    z <- z[match(panels, z$panel), ]
    x <- seq_along(panels) + (i - (length(methods) + 1) / 2) * .13
    lo <- z[[paste0(metric, '_lo')]]; hi <- z[[paste0(metric, '_hi')]]
    ok <- is.finite(lo) & is.finite(hi)
    graphics::segments(x[ok], lo[ok], x[ok], hi[ok], col = METHOD_COLOURS[methods[i]])
    graphics::points(x, z[[metric]], pch = METHOD_SYMBOLS[methods[i]], col = METHOD_COLOURS[methods[i]])
  }
}

# Curves over amplitude with shaded 95% intervals (one panel).
curve_panel <- function(s, metric, methods, ylab, ylim, reference = NULL) {
  xr <- range(s$a); if (diff(xr) == 0) xr <- xr + c(-.2, .2)
  graphics::plot(NA, xlim = xr, ylim = ylim, xlab = 'Signal amplitude a', ylab = ylab, bty = 'l')
  if (!is.null(reference)) graphics::abline(h = reference, lty = 2, col = 'grey50')
  for (m in methods) {
    z <- s[s$method == m, ]; z <- z[order(z$a), ]
    lo <- z[[paste0(metric, '_lo')]]; hi <- z[[paste0(metric, '_hi')]]
    ok <- is.finite(lo) & is.finite(hi)
    if (sum(ok) > 1) graphics::polygon(c(z$a[ok], rev(z$a[ok])), c(lo[ok], rev(hi[ok])),
                                       col = grDevices::adjustcolor(METHOD_COLOURS[m], alpha.f = .12), border = NA)
    graphics::lines(z$a, z[[metric]], type = 'o', pch = METHOD_SYMBOLS[m], col = METHOD_COLOURS[m], lwd = 1.5)
  }
}

# Power (left) and mean FDP (right) for the signal cases.
draw_main_figure <- function(summary, cfg, pdf_path, png_path) {
  s <- summary[summary$a > 0, ]
  methods <- intersect(cfg$methods, s$method)
  panels <- intersect(cfg$panels, s$panel)
  labels <- unlist(cfg$plot_labels)
  fdr_top <- max(.15, s$FDR, s$FDR_hi, na.rm = TRUE) * 1.1
  curves <- length(panels) == 1 && length(unique(s$a)) > 1
  draw <- function() {
    graphics::layout(matrix(c(1, 1, 2, 3), 2, byrow = TRUE), heights = c(.15, 1))
    graphics::par(cex = .85, las = 1)
    legend_strip(methods, labels)
    graphics::par(mar = c(if (curves) 4.5 else 3.5, 4, .5, 1))
    if (curves) {
      curve_panel(s, 'power', methods, 'Power', c(0, 1))
      curve_panel(s, 'FDR', methods, 'Mean FDP', c(0, fdr_top), reference = cfg$q)
    } else {
      dot_panel(s, 'power', panels, methods, 'Power', c(0, 1))
      dot_panel(s, 'FDR', panels, methods, 'Mean FDP', c(0, fdr_top), reference = cfg$q)
    }
  }
  save_both(draw, pdf_path, png_path, width = 8, height = if (curves) 3.6 else 3.2)
}

# Null case a = 0: probability of selecting at least one SNP (= FDR when no
# SNP is causal), with exact binomial intervals. Dashed line: q.
draw_null_figure <- function(summary, cfg, pdf_path, png_path) {
  s <- summary[summary$a == 0, ]
  methods <- intersect(cfg$methods, s$method)
  panels <- intersect(cfg$panels, s$panel)
  labels <- unlist(cfg$plot_labels)
  top <- max(.2, s$FDR_hi, na.rm = TRUE) * 1.1
  draw <- function() {
    graphics::layout(matrix(1:2, 2), heights = c(.18, 1))
    graphics::par(cex = .85, las = 1)
    legend_strip(methods, labels)
    graphics::par(mar = c(3.5, 4, .5, 1))
    dot_panel(s, 'FDR', panels, methods, 'P(any selection), a = 0', c(0, top), reference = cfg$q)
  }
  save_both(draw, pdf_path, png_path, width = max(7, 1.6 * length(panels) + 2), height = 3.2)
}
