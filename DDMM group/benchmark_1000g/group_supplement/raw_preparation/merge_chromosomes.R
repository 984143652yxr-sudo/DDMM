# Replacement for utils/merge_chromosomes.sh (which calls plink --merge-list):
# column-binds the snpknock2 outputs of two chromosomes into one PLINK file set
# and writes the group file without its header, as the authors' script does.
#   Rscript merge_chromosomes.R PREFIX_CHR21 PREFIX_CHR22 OUT_PREFIX
suppressMessages(library(snpStats))
a <- commandArgs(TRUE)
read1 <- function(p) read.plink(paste0(p, '.bed'), paste0(p, '.bim'), paste0(p, '.fam'))
x <- read1(a[1]); y <- read1(a[2])
stopifnot(identical(rownames(x$genotypes), rownames(y$genotypes)))
G <- cbind(x$genotypes, y$genotypes)
map <- rbind(x$map, y$map)
write.plink(a[3], snp.major = TRUE, snps = G, pedigree = x$fam$pedigree, id = x$fam$member,
            father = x$fam$father, mother = x$fam$mother, sex = x$fam$sex, phenotype = x$fam$affected,
            chromosome = map$chromosome, genetic.distance = map$cM, position = map$position,
            allele.1 = map$allele.1, allele.2 = map$allele.2)
grp <- rbind(read.table(paste0(a[1], '_grp.txt'), header = TRUE), read.table(paste0(a[2], '_grp.txt'), header = TRUE))
write.table(grp, paste0(a[3], '_grp.txt'), quote = FALSE, row.names = FALSE, col.names = FALSE)
