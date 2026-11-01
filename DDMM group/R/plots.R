# Figures (base R graphics). Colours: Paul Tol's bright palette, as in the
# earlier group_comparison figures.

track_colour <- function(tracks) {
  base <- c(simes = '#4477AA', rank_effective = '#228833', rank_3 = '#CCBB44',
            rank_full = '#AA3377', knockoff = '#EE6677')
  extra <- c('#66CCEE', '#BBBBBB', '#000000')
  out <- base[tracks]
  out[is.na(out)] <- extra[seq_len(sum(is.na(out)))]
  names(out) <- tracks
  out
}
track_symbol <- function(tracks) {
  s <- c(simes = 15, rank_effective = 16, rank_3 = 17, rank_full = 18, knockoff = 4)[tracks]
  s[is.na(s)] <- 1
  names(s) <- tracks
  s
}

save_both <- function(draw, pdf_path, png_path, width, height) {
  grDevices::pdf(pdf_path, width = width, height = height); draw(); grDevices::dev.off()
  grDevices::png(png_path, width = width * 300, height = height * 300, res = 300); draw(); grDevices::dev.off()
}

legend_strip <- function(tracks, labels) {
  graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()
  width <- max(graphics::strwidth(labels[tracks], cex = .8)) * 1.1
  graphics::legend('center', legend = labels[tracks], col = track_colour(tracks),
                   pch = track_symbol(tracks), lty = 1, horiz = TRUE, bty = 'n', cex = .8,
                   text.width = width)
}

# One metric over resolutions (coarse -> fine), a line per track.
metric_panel <- function(tab, metric, lo, hi, resolutions, tracks, ylab, ylim, reference = NULL, groups = NULL) {
  x <- seq_along(resolutions)
  graphics::plot(NA, xlim = c(.7, length(x) + .3), ylim = ylim, xaxt = 'n', xlab = '', ylab = ylab, bty = 'l')
  # Tick labels drawn with text() so that none is dropped on narrow panels.
  graphics::axis(1, x, labels = FALSE, tck = -.02)
  usr <- graphics::par('usr')
  graphics::text(x, usr[3] - .045 * diff(usr[3:4]), adj = c(.5, 1), xpd = NA, cex = .72,
                 if (is.null(groups)) resolutions else paste0(resolutions, '\n', groups))
  graphics::title(xlab = if (is.null(groups)) 'Resolution' else 'Resolution (number of groups)',
                  line = 3.1, cex.lab = .85)
  if (!is.null(reference)) graphics::abline(h = reference, lty = 2, col = 'grey50')
  offset <- (seq_along(tracks) - (length(tracks) + 1) / 2) * .06
  for (i in seq_along(tracks)) {
    z <- tab[tab$track == tracks[i], ]
    z <- z[match(resolutions, z$resolution), ]
    xx <- x + offset[i]
    if (!is.null(lo)) {
      ok <- is.finite(z[[lo]]) & is.finite(z[[hi]])
      graphics::segments(xx[ok], z[[lo]][ok], xx[ok], z[[hi]][ok], col = track_colour(tracks[i]))
    }
    graphics::lines(xx, z[[metric]], col = track_colour(tracks[i]), lwd = 1.3)
    graphics::points(xx, z[[metric]], col = track_colour(tracks[i]), pch = track_symbol(tracks[i]))
  }
}

group_counts <- function(tab, resolutions) {
  g <- tab$groups[match(resolutions, tab$resolution)]
  ifelse(is.na(g), '', as.character(g))
}

# Replicate study: (a) power and (b) mean FDP under signal, (c) P(any discovery) under the null.
draw_replicates <- function(summary, resolutions, q, labels, pdf_path, png_path) {
  tracks <- ordered_tracks(unique(summary$track))
  sig <- summary[summary$case == 'signal', ]; nul <- summary[summary$case == 'null', ]
  panels <- (nrow(sig) > 0) * 2 + (nrow(nul) > 0)
  top_fdr <- max(.2, sig$FDR_hi, sig$FDR, na.rm = TRUE) * 1.08
  top_null <- max(.3, nul$any_hi, na.rm = TRUE) * 1.05
  draw <- function() {
    graphics::layout(matrix(c(rep(1, panels), seq_len(panels) + 1), 2, byrow = TRUE), heights = c(.13, 1))
    graphics::par(cex = .85, las = 1)
    legend_strip(tracks, labels)
    graphics::par(mar = c(4.2, 4.2, 2, .8))
    if (nrow(sig)) {
      metric_panel(sig, 'power', 'power_lo', 'power_hi', resolutions, tracks, 'Locus power', c(0, 1),
                   groups = group_counts(sig, resolutions))
      graphics::title(sprintf('a  Signal, %d replicates', max(sig$replicates)), adj = 0, font.main = 2, cex.main = .95)
      metric_panel(sig, 'FDR', 'FDR_lo', 'FDR_hi', resolutions, tracks, 'Mean FDP', c(0, top_fdr), q,
                   groups = group_counts(sig, resolutions))
      graphics::title('b  Realised FDR', adj = 0, font.main = 2, cex.main = .95)
    }
    if (nrow(nul)) {
      metric_panel(nul, 'any_discovery', 'any_lo', 'any_hi', resolutions, tracks,
                   'P(any discovery)', c(0, top_null), q, groups = group_counts(nul, resolutions))
      graphics::title(sprintf('%s  Global null, %d replicates', if (nrow(sig)) 'c' else 'a', max(nul$replicates)),
                      adj = 0, font.main = 2, cex.main = .95)
    }
  }
  save_both(draw, pdf_path, png_path, width = 3.9 * panels, height = 3.9)
}

# Shipped phenotype: power and FDP per resolution, DDMM tracks and the published knockoffs.
draw_shipped <- function(tab, resolutions, q, labels, pdf_path, png_path) {
  tracks <- ordered_tracks(unique(tab$track))
  tab$FDR <- tab$FDP
  draw <- function() {
    graphics::layout(matrix(c(1, 1, 2, 3), 2, byrow = TRUE), heights = c(.13, 1))
    graphics::par(cex = .85, las = 1)
    legend_strip(tracks, labels)
    graphics::par(mar = c(4.2, 4.2, 2, .8))
    metric_panel(tab, 'power', NULL, NULL, resolutions, tracks, 'Locus power', c(0, 1),
                 groups = group_counts(tab, resolutions))
    graphics::title('a  Shipped phenotype: power', adj = 0, font.main = 2, cex.main = .95)
    metric_panel(tab, 'FDR', NULL, NULL, resolutions, tracks, 'FDP', c(0, max(.25, tab$FDP) * 1.08), q,
                 groups = group_counts(tab, resolutions))
    graphics::title('b  Shipped phenotype: FDP', adj = 0, font.main = 2, cex.main = .95)
  }
  save_both(draw, pdf_path, png_path, width = 9, height = 3.9)
}

# Discovery map of the shipped phenotype (the Manhattan-style figure of the group
# section): one panel row per method, one column per chromosome, a row of bars
# per resolution; filled bars = true discoveries, hatched = false, dotted lines =
# causal variants. Underneath, the structure of the groups: the LD dendrogram
# when the partitions come from adjclust, or the nesting of the published
# partitions otherwise.
draw_discovery_map <- function(layout, intervals, tracks, labels, pdf_path, png_path) {
  chrs <- sort(unique(layout$chr))
  resolutions <- layout$resolutions
  M <- length(resolutions)
  causal <- layout$causal
  structure_panel <- if (!is.null(layout$trees)) 'dendrogram' else if (partitions_nested(layout)) 'nesting' else 'none'
  rows <- length(tracks) + (structure_panel != 'none')
  is_tp <- function(ch, a, b) any(causal$CHR == ch & causal$BP >= a & causal$BP <= b)
  draw <- function() {
    graphics::layout(matrix(seq_len(rows * length(chrs)), rows, length(chrs), byrow = TRUE),
                     heights = c(rep(1.12, length(tracks)), if (structure_panel != 'none') .92))
    for (ti in seq_along(tracks)) for (ci in seq_along(chrs)) {
      ch <- chrs[ci]; col <- track_colour(tracks[ti])
      xr <- range(layout$pos[layout$chr == ch] / 1e6); cz <- causal$BP[causal$CHR == ch] / 1e6
      graphics::par(mar = c(.4, if (ci == 1) 9.5 else 1, 2.8, 1))
      graphics::plot(NA, xlim = xr, ylim = c(.4, M + .6), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')
      graphics::abline(v = cz, col = 'firebrick', lty = 3, lwd = .8)
      for (m in seq_len(M)) {
        row <- m                       # coarsest at the bottom, finest on top
        iv <- intervals[[tracks[ti]]][[resolutions[m]]]
        iv <- iv[iv$chr == ch, , drop = FALSE]
        for (i in seq_len(nrow(iv))) {
          tp <- is_tp(ch, iv$start[i], iv$end[i])
          hw <- max((iv$end[i] - iv$start[i]) / 1e6, diff(xr) * .003) / 2
          mid <- (iv$start[i] + iv$end[i]) / 2e6
          graphics::rect(mid - hw, row - .34, mid + hw, row + .34, col = if (tp) col else NA, border = col,
                         lwd = if (tp) .4 else 1, density = if (tp) NULL else 22)
        }
        graphics::abline(h = row - .5, col = 'grey90', lwd = .4)
      }
      graphics::mtext(sprintf('Chromosome %d', ch), side = 3, line = .9, cex = .95, font = 2)
      if (ci == 1) {
        graphics::axis(2, at = M:1, labels = layout$row_labels, las = 1, cex.axis = .8, tick = FALSE, line = -.3)
        graphics::mtext(letters[ti], side = 3, line = 1.3, adj = -.2, font = 2, cex = 1.3)
        graphics::mtext(labels[tracks[ti]], side = 3, line = .2, adj = 0, cex = .8, font = 3, col = col)
        graphics::mtext('Resolution', side = 2, line = 8, cex = .85)
        graphics::legend('topright', c('causal variant', 'true discovery', 'false discovery'),
                         col = c('firebrick', col, col), lty = c(3, NA, NA), pch = c(NA, 15, 0), bty = 'n', cex = .7)
      }
    }
    if (structure_panel != 'none') for (ci in seq_along(chrs)) {
      ch <- chrs[ci]; xr <- range(layout$pos[layout$chr == ch] / 1e6); cz <- causal$BP[causal$CHR == ch] / 1e6
      graphics::par(mar = c(3.4, if (ci == 1) 9.5 else 1, 2.2, 1))
      if (structure_panel == 'dendrogram') draw_ld_dendrogram(layout, ch, xr) else draw_partition_nesting(layout, ch, xr, ci == 1)
      graphics::abline(v = cz, col = 'firebrick', lty = 3, lwd = .6)
      graphics::axis(1, cex.axis = .9)
      graphics::mtext(sprintf('Chromosome %d position (Mb)', ch), side = 1, line = 2.3, cex = .85)
      if (ci == 1) {
        graphics::mtext(letters[length(tracks) + 1], side = 3, line = .7, adj = -.2, font = 2, cex = 1.3)
        graphics::mtext(if (structure_panel == 'dendrogram') 'LD dendrogram (adjclust on |r|) that defines the DDMM groups'
                        else 'Nesting of the published partitions, tested by both methods',
                        side = 3, line = -.1, adj = 0, cex = .72, font = 3)
        if (structure_panel == 'dendrogram') graphics::mtext('1 - |r|', side = 2, line = 8, cex = .8)
      }
    }
  }
  save_both(draw, pdf_path, png_path, width = 14,
            height = (560 * length(tracks) + if (structure_panel != 'none') 460 else 0) / 150)
}

draw_ld_dendrogram <- function(layout, ch, xr) {
  hc <- layout$trees[[as.character(ch)]]
  leaf <- layout$pos[layout$chr == ch] / 1e6
  node <- numeric(nrow(hc$merge))
  xof <- function(id) if (id < 0) leaf[-id] else node[id]
  for (k in seq_len(nrow(hc$merge))) node[k] <- (xof(hc$merge[k, 1]) + xof(hc$merge[k, 2])) / 2
  graphics::plot(NA, xlim = xr, ylim = c(0, max(hc$height) * 1.05), xaxt = 'n', yaxt = 'n',
                 xlab = '', ylab = '', bty = 'n')
  for (k in seq_len(nrow(hc$merge))) {
    l <- xof(hc$merge[k, 1]); r <- xof(hc$merge[k, 2]); h <- hc$height[k]
    hl <- if (hc$merge[k, 1] < 0) 0 else hc$height[hc$merge[k, 1]]
    hr <- if (hc$merge[k, 2] < 0) 0 else hc$height[hc$merge[k, 2]]
    graphics::segments(c(l, r, l), c(hl, hr, h), c(l, r, r), c(h, h, h), col = 'grey55', lwd = .25)
  }
}

# Every group of a finer resolution lies inside one group of the next coarser one.
partitions_nested <- function(layout) {
  a <- layout$assignments[layout$resolutions]
  length(a) > 1 && all(vapply(seq_len(length(a) - 1), function(k)
    all(tapply(a[[k]], a[[k + 1]], function(z) length(unique(z))) == 1), TRUE))
}

# Tree of the nested partitions: variants at height 0, each resolution one level
# up (finest group resolution at height 1, coarsest at the top); a group is
# joined to the group that contains it at the next coarser resolution.
draw_partition_nesting <- function(layout, ch, xr, label_levels) {
  on <- layout$chr == ch
  pos <- layout$pos[on] / 1e6
  levels <- rev(layout$resolutions)                       # fine -> coarse
  assign <- lapply(layout$assignments[levels], function(a) a[on])
  graphics::plot(NA, xlim = xr, ylim = c(0, length(levels) - .5), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')
  x_child <- pos                                          # current level: one x per unit
  unit <- seq_along(pos)                                  # unit index of each variant
  for (L in seq_along(levels)[-1]) {
    parent_of_variant <- assign[[L]]
    parents <- unique(parent_of_variant)
    unit_parent <- tapply(parent_of_variant, unit, `[`, 1)
    x_parent <- tapply(x_child, unit_parent[as.character(seq_along(x_child))], mean)
    h0 <- L - 2; h1 <- L - 1
    graphics::segments(x_child, h0, x_child, h1, col = 'grey55', lwd = if (L == 2) .15 else .4)
    span <- tapply(x_child, unit_parent[as.character(seq_along(x_child))], range)
    for (p in names(span)) graphics::segments(span[[p]][1], h1, span[[p]][2], h1, col = 'grey55', lwd = .4)
    unit <- match(parent_of_variant, as.integer(names(x_parent)))
    x_child <- as.numeric(x_parent)
  }
  if (label_levels) graphics::axis(2, at = seq_along(levels)[-1] - 1, labels = levels[-1], las = 1,
                                   cex.axis = .7, tick = FALSE, line = -.6)
}
