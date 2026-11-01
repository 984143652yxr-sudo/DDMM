# Release validation

Validated on 2026-11-01 with R 4.5.0, data.table 1.18.4, and jsonlite 2.0.0.

- `Rscript reproduce.R --check-only`: all 1,361 input checksums and all three
  panels passed; 60 phenotype selections per method and 420 KnockoffGWAS
  discovery files per panel were present.
- `Rscript reproduce.R`, invoked by its full path from outside the repository:
  completed successfully. All six regenerated CSV files agreed with the
  reference values at tolerance `1e-8`; all four Markdown/LaTeX tables matched
  line by line. Both figure formats were regenerated.
- `Rscript report.R knockoffgwas`, also invoked from outside the repository:
  regenerated the toy-example reports from the bundled saved fits.
- `Rscript run.R --help`: the original fitting entry point loaded successfully.
- The release contains no filesystem symbolic links or raw haplotype/VCF files.

Model fits and knockoff generation were not rerun for this packaging check.
The installed data.table package reported that it was built under R 4.5.2;
all checks above passed on R 4.5.0.
