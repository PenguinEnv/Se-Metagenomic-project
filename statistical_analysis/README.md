# Statistical analysis (MAG and genus level)

Community statistics for the Se bioelectrochemical reactor metagenomes, written in R.

- Design: condition **OC** (open circuit) vs **BES** (bioelectrochemical) x
  location **electrode** (BF) vs **tube** (TB), 3 biological replicates per combination (n = 12)
- Input: CoverM quantification (276 MAGs x 12 samples), GTDB taxonomy
- Analyses: two-way PERMANOVA, CLR-PCA / Bray-Curtis PCoA, ALDEx2 differential abundance,
  genus-level OC vs BES comparison

**Results and interpretation: [`results/REPORT.md`](results/REPORT.md)**

## Layout

```
scripts/se_bes_metagenome_stats.R   Full analysis (one run produces all tables and figures)
prompts/analysis_prompt.md          The prompt used to request the analysis + analysis decisions
env/environment.yml                 R / package versions used (conda)
data/                               Put the input TSVs here (not committed; see .gitignore)
results/
  REPORT.md                         Summary and interpretation
  tables/                           PERMANOVA, betadisper, ALDEx2, genus summary tables (TSV)
  figures/                          Ordination, ALDEx2 effect, genus comparison figures (PDF + PNG)
  matrices/                         MAG/genus x sample matrices (same info as raw data; not committed)
  sessionInfo.txt                   R session info, seed, parameters
```

## Running

Run all commands from this `statistical_analysis/` directory.

```bash
# 1) Environment (conda/mamba)
micromamba create -f env/environment.yml
# or in R: install.packages(c("vegan", "ggplot2", "ggrepel", "patchwork"));
#          BiocManager::install("ALDEx2")

# 2) Put the input files in data/ and run (~15 min because of 1000 ALDEx2 MC instances)
Rscript scripts/se_bes_metagenome_stats.R \
    data/coverm_all.tsv data/MAG_abundance_taxonomy.tsv results

# Quick check: fewer ALDEx2 Monte-Carlo instances
ALDEX_MC=128 Rscript scripts/se_bes_metagenome_stats.R
```

Input format
- `coverm_all.tsv`: `Genome` + four columns per sample:
  `"<sample>_R1.clean.fq.gz Relative Abundance (%)"`, `Mean`, `Trimmed Mean`, `Read Count`.
  The script removes the `unmapped` row.
- `MAG_abundance_taxonomy.tsv`: `MAG`, `domain` ... `species`, then per-sample abundance.
- Sample names must match `^(OC|BES)(BF|TB)\d+$` (e.g. `OCBF01`, `BESTB03`).
