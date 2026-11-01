"""Real-genotype version of the KnockoffGWAS toy example, from 1000 Genomes Phase 3.

For each cohort (EUR: the 503 EUR samples; MIX: 168 EUR, 168 AFR and 168 EAS
samples drawn at random with seed 2019, n = 504, the size of EUR; ALL: all 2504
samples) and each of
chromosomes 21 and 22 (GRCh37), keep biallelic SNPs with MAF >= 0.05 in the
cohort, starting at 30.0 Mb, thinned so that consecutive kept SNPs are at least
2 kb apart (array-like density), and stop at 1000 SNPs per chromosome, as in the
toy data (2 x 1000 variants). Genetic positions are interpolated from the
HapMap-based hg19 map shipped with BOLT-LMM.

Writes, per cohort, into OUT/<cohort>/:
  hap_chr{21,22}.haps/.legend/.sample   phased haplotypes for snpknock2 (--haps)
  geno_chr{21,22}.bim                   variant list (tab-delimited, as in KnockoffGWAS)
  map_chr{21,22}.txt                    genetic map at every variant, KnockoffGWAS format
  variants_qc.txt, samples.tsv          variant ids; samples with population and sex
  dosage.npy is not written; R reads the haps directly.

usage: python3 prepare_1000g.py VCF_DIR PANEL MAP_GZ OUT [COHORT ...]   (default: EUR MIX ALL)
"""
import gzip, os, subprocess, sys
import numpy as np

vcf_dir, panel_file, map_gz, out = sys.argv[1:5]
START, MIN_GAP, P_CHR, MAF = 30_000_000, 2000, 1000, 0.05
CHRS = [21, 22]

panel = [l.rstrip('\n').split('\t') for l in open(panel_file)][1:]
rng = np.random.default_rng(2019)
mix = set()
for sp in ('EUR', 'AFR', 'EAS'):
    ids_sp = [r[0] for r in panel if r[2] == sp]
    mix |= set(rng.choice(ids_sp, size=168, replace=False).tolist())
cohorts = {'EUR': [r for r in panel if r[2] == 'EUR'],
           'MIX': [r for r in panel if r[0] in mix],          # panel order
           'ALL': panel}
wanted = sys.argv[5:] or list(cohorts)
cohorts = {k: v for k, v in cohorts.items() if k in wanted}

gmap = {c: [] for c in CHRS}
with gzip.open(map_gz, 'rt') as f:
    next(f)
    for line in f:
        c, pos, rate, cm = line.split()
        if c.isdigit() and int(c) in gmap:
            gmap[int(c)].append((int(pos), float(rate), float(cm)))
for c in CHRS:
    gmap[c] = np.array(gmap[c])

for name, rows in cohorts.items():
    d = os.path.join(out, name); os.makedirs(d, exist_ok=True)
    ids = [r[0] for r in rows]
    open(os.path.join(d, 'samples.txt'), 'w').write('\n'.join(ids) + '\n')
    with open(os.path.join(d, 'samples.tsv'), 'w') as f:
        f.write('sample\tpop\tsuper_pop\tsex\n')
        for r in rows:
            f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{1 if r[3] == 'male' else 2}\n")
    qc = []
    for c in CHRS:
        vcf = os.path.join(vcf_dir, f'ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
        cmd = (f"bcftools view -S {d}/samples.txt --force-samples -r {c}:{START}-{START + 20_000_000} "
               f"-m2 -M2 -v snps {vcf} -Ou | bcftools view -q {MAF}:minor -Ou | "
               f"bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT[\\t%GT]\\n'")
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
        kept, last, seen = [], -10**9, set()
        g = gmap[c]
        for line in proc.stdout:
            t = line.rstrip('\n').split('\t')
            pos = int(t[1])
            if pos < last + MIN_GAP or pos < g[0, 0] or pos > g[-1, 0]:
                continue
            vid = t[2] if t[2] != '.' else f'{c}:{pos}'
            if vid in seen or ';' in vid:
                continue
            gts = t[5:]
            if any(len(x) != 3 or x[1] != '|' for x in gts):
                continue
            h = np.array([[int(x[0]), int(x[2])] for x in gts], dtype=np.int8).reshape(-1)
            kept.append((pos, vid, t[3], t[4], h)); seen.add(vid); last = pos
            if len(kept) == P_CHR:
                break
        proc.kill()
        assert len(kept) == P_CHR, (name, c, len(kept))
        pos = np.array([k[0] for k in kept])
        cm = np.interp(pos, g[:, 0], g[:, 2])
        rate = np.interp(pos, g[:, 0], g[:, 1])
        with open(os.path.join(d, f'hap_chr{c}.haps'), 'w') as f:
            for k in kept:
                f.write(' '.join(map(str, k[4])) + '\n')
        with open(os.path.join(d, f'hap_chr{c}.legend'), 'w') as f:
            f.write('chr id cm pos a0 a1\n')
            for k, m in zip(kept, cm):
                f.write(f'{c} {k[1]} {m:.6f} {k[0]} {k[2]} {k[3]}\n')
        with open(os.path.join(d, f'hap_chr{c}.sample'), 'w') as f:
            f.write('ID_1 ID_2 missing sex\n0 0 0 D\n')
            for r in rows:
                f.write(f"{r[0]} {r[0]} 0 {1 if r[3] == 'male' else 2}\n")
        with open(os.path.join(d, f'geno_chr{c}.bim'), 'w') as f:
            for k in kept:
                f.write(f'{c}\t{k[1]}\t0\t{k[0]}\t{k[2]}\t{k[3]}\n')
        with open(os.path.join(d, f'map_chr{c}.txt'), 'w') as f:
            f.write('Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)\n')
            for p_, r_, m_ in zip(pos, rate, cm):
                f.write(f'chr{c}\t{p_}\t{r_:.6f}\t{m_:.8f}\n')
        qc += [k[1] for k in kept]
        print(f'{name} chr{c}: {len(kept)} SNPs, {pos[0]}-{pos[-1]} bp ({(pos[-1]-pos[0])/1e6:.2f} Mb, '
              f'{cm[-1]-cm[0]:.2f} cM), n={len(ids)}', flush=True)
    open(os.path.join(d, 'variants_qc.txt'), 'w').write('\n'.join(qc) + '\n')
