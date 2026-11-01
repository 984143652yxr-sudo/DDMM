# Data sources

- **Toy example:** KnockoffGWAS example genotypes, phenotypes, partitions, and
  published discovery lists (Sesia et al., 2021, *PNAS* 118, e2105841118).
- **1000G benchmark:** 2,000 SNPs from chromosomes 21 and 22 of the 1000 Genomes
  Phase 3 panel. EUR has 503 individuals; MIX has 504 (168 each EUR/AFR/EAS);
  ALL has 2,504. Outcomes are simulated: 40 causal variants in 20 loci,
  amplitudes 8 and 4, and a global null, with 20 replicates per setting.
- **Saved selections:** DDMM, KnockoffGWAS, and SuSiE outputs retained as inputs
  for rescoring. SuSiE uses L = 40, 95%/90% credible sets, and PIP > 0.9
  (Wang et al., 2020, *JRSS B* 82, 1273–1300).

The compact `design.rds` files contain genotypes, variant/sample metadata,
and simulated causal truth. Raw VCFs, phased haplotypes, and knockoff matrices
are excluded. Upstream data and software retain their own attribution and terms.

## Download original 1000G data

The bundled genotype subset is sufficient for `run_1000g.R`. For the source
VCFs, use the [1000 Genomes Phase 3 release (GRCh37)](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)
or the [IGSR data portal](https://www.internationalgenome.org/data/).

Download chromosomes 21 and 22, their indexes, and sample information:

```sh
mkdir -p raw_1000g
cd raw_1000g
base=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
for chr in 21 22; do
  file="ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  curl -fLO "$base/$file"
  curl -fLO "$base/$file.tbi"
done
curl -fLO "$base/integrated_call_samples_v3.20130502.ALL.panel"
```

The full VCFs are about 400 MB combined. The bundled panels contain only the
selected samples and variants, identified by `samples.txt` and `variants_qc.txt`;
`design.rds$X` stores their unphased 0/1/2 allele dosages. These dosage matrices
support DDMM and SuSiE. Regenerating HMM knockoffs also needs phased haplotypes,
a genetic map, and the [KnockoffGWAS software](https://github.com/msesia/knockoffgwas).
