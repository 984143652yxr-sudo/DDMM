## Tables of the group benchmark (cohorts EUR, MIX, ALL) from results/summary.csv,
## contrasts.csv and design_summary.csv:
##   results/tables.md            everything in Markdown
##   results/table_main.tex       DDMM (effective rank, BH) against KnockoffGWAS per resolution,
##                                with SuSiE (CS95, CS90, PIP > 0.9) as reference rows
##   results/table_tracks.tex     the four DDMM statistics, FDR with LD tags counted as true,
##                                and the null any-discovery rate
##   results/table_data.tex       the three cohorts and their partitions
##   BENCH=<folder> Rscript make_tables.R
suppressMessages(library(data.table))
B <- Sys.getenv("BENCH"); R <- file.path(B, "results")
S <- fread(file.path(R, "summary.csv")); CT <- fread(file.path(R, "contrasts.csv")); DZ <- fread(file.path(R, "design_summary.csv"))
COH <- c("EUR", "MIX", "ALL"); RES <- sprintf("res%d", 6:0)
NN <- c(EUR = "$n=503$", MIX = "$n=504$", ALL = "$n=2{,}504$")
f2 <- function(x) if (!length(x)) "--" else ifelse(is.na(x), "--", sprintf("%.2f", x))
sgn <- function(x) ifelse(is.na(x), "--", sprintf("$%s%.2f$", ifelse(x >= 0, "+", "-"), abs(x)))
get <- function(coh, cs, res, meth, rl) S[cohort == coh & case == cs & resolution == res & method == meth & rule == rl]
pf <- function(d) if (!nrow(d)) "--" else sprintf("%s / %s", f2(d$power), f2(d$FDR))
dif <- function(coh, cs, res, vs = "KnockoffGWAS") { d <- CT[cohort == coh & case == cs & resolution == res & versus == vs]
  if (!nrow(d)) "--" else sprintf("%s (%s)", sgn(d$diff), f2(d$se)) }
DD <- c("DDMM_rank_effective", "BH"); KG <- c("KnockoffGWAS", "knockoff+")
SU <- list(CS95 = c("SuSiE_CS95", "CS95"), CS90 = c("SuSiE_CS90", "CS90"), PIP90 = c("SuSiE_PIP90", "PIP90"))
SU_LAB <- c(CS95 = "SuSiE, 95\\% CS", CS90 = "SuSiE, 90\\% CS", PIP90 = "SuSiE, PIP $>0.9$")
wkb <- function(x) if (!length(x) || is.na(x)) "--" else if (x < 1) "$<1$" else sprintf("%.0f", x)

## ---------------- main table ----------------
tex <- c("\\begin{tabular}{llccccccc}", "\\toprule",
  " & & & \\multicolumn{3}{c}{Amplitude 8} & \\multicolumn{3}{c}{Amplitude 4}\\\\", "\\cmidrule(lr){4-6}\\cmidrule(lr){7-9}",
  "Cohort & Resolution (groups, width) & All groups & DDMM & KnockoffGWAS & Difference & DDMM & KnockoffGWAS & Difference\\\\", "\\midrule")
md <- c("# Group benchmark on 1000 Genomes haplotypes: EUR, MIX, ALL", "",
  "Power / realised FDR over 20 phenotypes, q = 0.10. DDMM = effective-rank chi-squared + BH; KGWAS = KnockoffGWAS knockoff+ (fresh knockoffs per replicate); same adjclust partitions closed under r^2 >= 0.8. Width = median group width (kb). All FDP = FDP of declaring every group. Diff = paired DDMM - KGWAS power (SE). SuSiE rows: one fit, credible sets at coverage 0.95 / 0.90 and single variants with PIP > 0.9, width = median credible-set width.", "")
for (coh in COH) {
  md <- c(md, sprintf("## %s", coh), "", "| Resolution | Groups | Width kb | All FDP | a8 DDMM | a8 KGWAS | a8 diff | a4 DDMM | a4 KGWAS | a4 diff | null any DDMM / KGWAS |", "|---|---|---|---|---|---|---|---|---|---|---|")
  for (i in seq_along(RES)) { r <- RES[i]; z <- DZ[cohort == coh & resolution == r]
    cells <- c(pf(get(coh, "a8", r, DD[1], DD[2])), pf(get(coh, "a8", r, KG[1], KG[2])), dif(coh, "a8", r),
               pf(get(coh, "a4", r, DD[1], DD[2])), pf(get(coh, "a4", r, KG[1], KG[2])), dif(coh, "a4", r))
    nul <- sprintf("%s / %s", f2(get(coh, "null", r, DD[1], DD[2])$any_discovery), f2(get(coh, "null", r, KG[1], KG[2])$any_discovery))
    md <- c(md, sprintf("| %s | %d | %s | %s | %s | %s |", r, z$groups, wkb(z$median_width_kb), f2(z$trivial_FDP), paste(cells, collapse = " | "), nul))
    tex <- c(tex, sprintf("%s & %s (%d, %s kb) & %s & %s\\\\", if (i == 1) coh else if (i == 2) NN[coh] else "", r, z$groups,
                          wkb(z$median_width_kb), f2(z$trivial_FDP), paste(cells, collapse = " & ")))
  }
  for (s in names(SU)) {
    a8 <- get(coh, "a8", "adaptive", SU[[s]][1], SU[[s]][2]); a4 <- get(coh, "a4", "adaptive", SU[[s]][1], SU[[s]][2])
    w <- if (nrow(a8)) wkb(a8$width_kb) else "--"
    md <- c(md, sprintf("| %s | -- | %s | -- | %s | | | %s | | | %s |", s, w, pf(a8), pf(a4), f2(get(coh, "null", "adaptive", SU[[s]][1], SU[[s]][2])$any_discovery)))
    tex <- c(tex, sprintf(" & %s (%s) & -- & \\multicolumn{3}{c}{%s} & \\multicolumn{3}{c}{%s}\\\\", SU_LAB[s], if (s == "PIP90") "single SNP" else paste(w, "kb"), pf(a8), pf(a4)))
  }
  md <- c(md, ""); if (coh != COH[length(COH)]) tex <- c(tex, "\\midrule")
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(R, "table_main.tex"))

## ---------------- DDMM statistics ----------------
TR <- c(rank_effective = "$\\chi^2_{K_g}$", rank_3 = "$\\chi^2_{3}$", rank_full = "$\\chi^2_{d_g}$", simes = "Simes")
tt <- c("\\begin{tabular}{lllcccccc}", "\\toprule",
  sprintf("Cohort & Amp. & Resolution & %s & KnockoffGWAS & Null any, $\\chi^2_{K_g}$ / KnockoffGWAS\\\\", paste(TR, collapse = " & ")), "\\midrule")
md <- c(md, "## DDMM statistics (power / FDR), BH", "", "| Cohort | Amp | Res | effective | K=3 | full | Simes | KGWAS | eff FDR untagged | null any eff / KGWAS |", "|---|---|---|---|---|---|---|---|---|---|")
for (coh in COH) for (a in c("a8", "a4")) { for (i in seq_along(RES)) { r <- RES[i]
    cells <- c(vapply(names(TR), function(t) pf(get(coh, a, r, paste0("DDMM_", t), "BH")), ""), pf(get(coh, a, r, KG[1], KG[2])))
    unt <- f2(get(coh, a, r, DD[1], DD[2])$FDR_untagged)
    nul <- sprintf("%s / %s", f2(get(coh, "null", r, DD[1], DD[2])$any_discovery), f2(get(coh, "null", r, KG[1], KG[2])$any_discovery))
    md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s |", coh, sub("a", "", a), r, paste(cells, collapse = " | "), unt, nul))
    tt <- c(tt, sprintf("%s & %s & %s & %s & %s\\\\", if (i == 1 && a == "a8") coh else if (i == 2 && a == "a8") NN[coh] else "",
                        if (i == 1) sub("a", "", a) else "", r, paste(cells, collapse = " & "), nul)) }
  if (!(coh == COH[length(COH)] && a == "a4")) tt <- c(tt, "\\midrule") }
tt <- c(tt, "\\bottomrule", "\\end{tabular}")
writeLines(tt, file.path(R, "table_tracks.tex"))

## ---------------- data table ----------------
des <- lapply(COH, function(coh) readRDS(file.path(B, "data", coh, "design.rds")))
names(des) <- COH
cm <- function(d) paste(sprintf("%.1f", tapply(d$map$cM, d$map$CHR, function(x) diff(range(x)))), collapse = "+")
mb <- function(d) paste(sprintf("%.1f", tapply(d$map$BP, d$map$CHR, function(x) diff(range(x)) / 1e6)), collapse = "+")
copy <- function(d, causal = FALSE) { X <- scale(d$X * 1.0); X[is.na(X)] <- 0; Rm <- abs(cor(X)); diag(Rm) <- 0
  Rm[outer(d$map$CHR, d$map$CHR, "!=")] <- 0; has <- apply(Rm, 1, max) > 0.99
  if (causal) mean(has[match(d$causal$SNP, d$map$SNP)]) else mean(has) }
row <- function(lab, v) sprintf("%s & %s\\\\", lab, paste(v, collapse = " & "))
td <- c("\\begin{tabular}{lccc}", "\\toprule", " & EUR & MIX & ALL\\\\", "\\midrule",
  row("Individuals", c("503 of European ancestry", "168 each of European,", "all 2{,}504, five")),
  row("", c("", "African and East Asian ancestry", "continental groups")),
  row("$n$", c("$503$", "$504$", "$2{,}504$")),
  row("$p$ (chromosomes 21 and 22)", rep("$2{,}000$", 3)),
  row("Physical length (Mb)", vapply(des, mb, "")),
  row("Genetic length (cM)", vapply(des, cm, "")),
  row("Variants with a copy at $|r|>0.99$", sprintf("$%.0f\\%%$", 100 * vapply(des, copy, 0))),
  row("Causal variants with such a copy", sprintf("$%.0f\\%%$", 100 * vapply(des, copy, 0, causal = TRUE))),
  row("Heritability, amplitude 8 and 4", vapply(des, function(d) sprintf("$%.2f$, $%.2f$", d$h2[1], d$h2[2]), "")),
  row("Groups, res6 to res0", vapply(COH, function(coh) paste(DZ[cohort == coh][match(RES, resolution)]$groups, collapse = ", "), "")),
  row("Median width (kb), res6 to res1", vapply(COH, function(coh) paste(vapply(DZ[cohort == coh][match(RES[-7], resolution)]$median_width_kb, wkb, ""), collapse = ", "), "")),
  row("All-groups FDP, res6 to res0", vapply(COH, function(coh) paste(f2(DZ[cohort == coh][match(RES, resolution)]$trivial_FDP), collapse = ", "), "")),
  "\\bottomrule", "\\end{tabular}")
writeLines(td, file.path(R, "table_data.tex"))
md <- c(md, "", "## Cohorts", "", td)
writeLines(md, file.path(R, "tables.md")); cat(md, sep = "\n")
