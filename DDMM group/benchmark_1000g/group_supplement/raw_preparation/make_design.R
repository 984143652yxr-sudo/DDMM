# Design of the 1000 Genomes benchmark for one cohort: genotypes, causal
# variants and phenotypes, following the authors' toy-data generator
# (misc/make_phenotypes.R of KnockoffZoom/KnockoffGWAS) as closely as possible.
#
#   Rscript make_design.R data/EUR
#
# Causal variants: the authors' choose_causal(20, 40) with set.seed(2019), i.e.
# 20 loci evenly spaced along each chromosome, two causal variants per locus
# (the second at a fixed distance drawn from U(0, 100) kb), one sign per locus,
# |effect| scale ~ U(0.1, 1.9).
# Phenotypes: y = 10 scale(Z)/sqrt(n) + scale(G_causal) beta + N(0, 1) with
# beta = amplitude * Sign * Scale / sqrt(n) and Z = 2 - sex, as in the toy data.
# Amplitude 8 is the authors' value; amplitude 4 is a weaker setting. Replicate r
# redraws the noise with seed 1000 + r; its null partner has the same noise and
# no genetic effect. Writes design.rds and phenotypes.tab.

suppressMessages({library(data.table)})
args <- commandArgs(TRUE)
dir <- args[1]
R_REP <- 20
AMPLITUDES <- c(8, 4)

samples <- fread(file.path(dir, 'samples.tsv'))
n <- nrow(samples)
read_chr <- function(ch) {
  leg <- fread(file.path(dir, sprintf('hap_chr%d.legend', ch)))
  H <- as.matrix(fread(file.path(dir, sprintf('hap_chr%d.haps', ch)), header = FALSE))   # variants x haplotypes
  stopifnot(ncol(H) == 2 * n, nrow(H) == nrow(leg))
  X <- t(H[, seq(1, 2 * n, 2)] + H[, seq(2, 2 * n, 2)])                                 # n x variants
  colnames(X) <- leg$id
  list(X = X, map = data.frame(CHR = ch, SNP = leg$id, BP = leg$pos, cM = leg$cm, A0 = leg$a0, A1 = leg$a1))
}
parts <- lapply(c(21, 22), read_chr)
X <- do.call(cbind, lapply(parts, `[[`, 'X')); storage.mode(X) <- 'integer'
map <- do.call(rbind, lapply(parts, `[[`, 'map'))
rownames(X) <- samples$sample

# ---- the authors' choose_causal(), verbatim apart from the input table ----
suppressMessages(library(dplyr))
Variants <- tibble(CHR = map$CHR, SNP = map$SNP, BP = map$BP)
chr.list <- c(21, 22)
set.seed(2019)
choose_causal <- function(n.loci, n.signals) {
    locus.size <- n.signals/n.loci
    Variants.num <- Variants %>% group_by(CHR) %>% summarise(N=n())
    Variants.num <- Variants.num %>% mutate(Share = N/sum(Variants.num$N), Loci = round(Share*n.loci))
    if(sum(Variants.num$Loci) != n.loci) {
        Variants.num$Loci[1] = Variants.num$Loci[1]+(n.loci-sum(Variants.num$Loci))
    }
    n.loci.chr <- Variants.num$Loci
    names(n.loci.chr) <- chr.list
    causal.dist <- round(sort(runif(locus.size, 0, 100)),2)
    Variants.causal <- tibble()
    for(chr in chr.list) {
        chr.idx <- as.character(chr)
        Variants.chr <- Variants %>% filter(CHR==chr)
        causal.idx <- round(seq(round(0.02*nrow(Variants.chr)),
                                round(0.98*nrow(Variants.chr)), length.out=n.loci.chr[chr.idx]))
        bp.first <- Variants.chr$BP[causal.idx]
        bp.next <- as.numeric(sapply(causal.dist[-1], function(dkb) {bp.first + 1000 * dkb}))
        bp.next <- matrix(bp.next, ncol=length(causal.dist)-1)
        locus.idx <- seq(length(causal.idx))
        for(j in 1:length(bp.first)) {
            if(length(causal.dist)>1) {
                for(k in 1:ncol(bp.next)) {
                    candidates <- order(abs(Variants.chr$BP-bp.next[j,k]))
                    candidates <- candidates[which(candidates>causal.idx[j])]
                    candidates <- setdiff(candidates, causal.idx)
                    causal.idx <- c(causal.idx, candidates[1])
                    locus.idx <- c(locus.idx, j)
                }
            }
        }
        Variants.causal.chr <- Variants.chr[causal.idx,]
        Variants.causal.chr$Locus <- locus.idx
        Variants.causal.chr <- Variants.causal.chr %>% arrange(BP) %>% distinct()
        Variants.causal <- rbind(Variants.causal, Variants.causal.chr)
    }
    Variants.causal <- Variants.causal %>%
        mutate(Locus=rep(seq(1,n.loci),each=locus.size)) %>%
        mutate(Sign=rep(2*rbinom(nrow(Variants.causal)/locus.size,1,0.5)-1,each=locus.size)) %>%
        mutate(Scale=round(runif(nrow(Variants.causal),0.1,1.9),3))
    return(Variants.causal)
}
causal <- as.data.frame(choose_causal(20, 40))
stopifnot(nrow(causal) == 40, length(unique(causal$Locus)) == 20)
# ---------------------------------------------------------------------------

G <- X[, causal$SNP]
Gs <- scale(G)
Z <- 2 - samples$sex
sex_term <- as.numeric(10 * scale(Z) / sqrt(n))
pheno <- data.table(FID = samples$sample, IID = samples$sample, sex = samples$sex)
genetic <- list()
for (a in AMPLITUDES) genetic[[as.character(a)]] <- as.numeric(Gs %*% (a * causal$Sign * causal$Scale / sqrt(n)))
for (r in seq_len(R_REP)) {
  set.seed(1000 + r)
  e <- rnorm(n)
  for (a in AMPLITUDES) pheno[[sprintf('y_a%d_r%02d', a, r)]] <- round(sex_term + genetic[[as.character(a)]] + e, 6)
  pheno[[sprintf('y_null_r%02d', r)]] <- round(sex_term + e, 6)
}
fwrite(pheno, file.path(dir, 'phenotypes.tab'), sep = '\t')
h2 <- sapply(genetic, function(g) var(g) / (var(g) + 1 + var(sex_term)))
saveRDS(list(X = X, map = map, samples = as.data.frame(samples), causal = causal, n_loci = 20,
             amplitudes = AMPLITUDES, replicates = R_REP, h2 = h2),
        file.path(dir, 'design.rds'))
cat(sprintf('%s: n = %d, p = %d, 40 causal variants in 20 loci, h2 (amplitude 8, 4) = %s\n',
            basename(dir), n, ncol(X), paste(round(h2, 3), collapse = ', ')))
print(causal)
