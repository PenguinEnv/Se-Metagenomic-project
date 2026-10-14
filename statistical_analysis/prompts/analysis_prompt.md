# Analysis prompt (Se bioelectrochemical metagenome statistics)

`scripts/se_bes_metagenome_stats.R` and `results/` in this folder were produced by giving the
prompt below to Claude Code (an AI coding assistant). It is kept here for reproducibility and
documentation. The original prompt was written in Korean; this is an English translation
(uploaded file paths replaced by file names).

- Date: 2026-10-14
- Input files: `coverm_all.tsv`, `MAG_abundance_taxonomy.tsv` (place them in `data/`)
- Run (from `statistical_analysis/`):
  `Rscript scripts/se_bes_metagenome_stats.R data/coverm_all.tsv data/MAG_abundance_taxonomy.tsv results`

---

## Prompt (English translation)

```text
@coverm_all.tsv @MAG_abundance_taxonomy.tsv Please run a statistical analysis of the Se
bioelectrochemical metagenome in R. I uploaded two files.

## Data
- coverm_all.tsv: quantification of 276 MAGs x 12 samples. Four columns per sample
  ("<sample> Relative Abundance (%)", "<sample> Mean", "<sample> Trimmed Mean",
  "<sample> Read Count"). There is an "unmapped" row. Sample name format:
  OCBF01_R1.clean.fq.gz etc.
- MAG_abundance_taxonomy.tsv: GTDB classification per MAG (domain to species) + relative
  abundance. First column MAG, then domain/phylum/class/order/family/genus/species, then
  per-sample abundance.

## Experimental design (n=12)
- condition: OC (open circuit) vs BES (bioelectrochemical)  <- sample name prefix
- location: electrode (BF) vs tube (TB)  <- BF/TB in the sample name
- 3 biological replicates per combination
- Samples: OCBF01-03, OCTB01-03, BESBF01-03, BESTB01-03

## Do everything at both the MAG level and the genus level
For the genus level, group by the genus column of MAG_abundance_taxonomy.tsv and sum
abundance/count of MAGs in the same genus to build a genus x sample matrix.
MAGs whose genus is unclassified/empty should not be dropped: label them like
'unclassified_<higher taxon>' or fall back to the family level.

## Analyses (each at MAG level and genus level)
1. Data parsing: split Relative Abundance and Read Count from coverm_all.tsv into separate
   matrices. Remove unmapped. Clean sample names (OCBF01 etc.). Build metadata
   (condition, location).
2. Two-way PERMANOVA: Bray-Curtis distance, adonis2(~ condition*location, 9999 perm)
   -> significance of condition / location / interaction. Both MAG and genus level.
3. Ordination (both MAG and genus level):
   - PCA after CLR transform (count + pseudocount)
   - Bray-Curtis PCoA
   Point colour = condition, shape = location, sample name labels.
4. ALDEx2 differential abundance (both MAG and genus level):
   - BES vs OC (integer Read Count, denom="all", test="t", effect=TRUE)
   - Extract significant features (wi.eBH < 0.05) and sort by effect size
   - Attach taxonomy (species name) to the MAG-level results
   - If possible, do location (electrode vs tube) the same way
5. Genus-level summary: bar chart or table comparing mean relative abundance of each genus
   in OC vs BES (so it is easy to see which genera dominate in BES)

## Outputs
- PERMANOVA result tables (MAG level, genus level)
- Ordination figures (2 at MAG level + 2 at genus level, PDF or PNG)
- ALDEx2 significant-result tables (MAG level with taxonomy, genus level separately,
  sorted by effect size)
- Genus-level OC vs BES comparison figure/table
- Interpretation of each result:
  * Are OC and BES communities significantly different? (Same conclusion at MAG vs genus level?)
  * Which species/genera increased significantly under BES (current supply)?
  * What is the difference between electrode and tube?

n=12 is small, so also report effect sizes, and comment on whether the MAG-level and
genus-level results agree or differ.

Save the prompt used for the computation separately so it can be uploaded to GitHub.
```

---

## Analysis decisions not specified in the prompt (made in the script)

| Item | Decision | Reason |
|---|---|---|
| Bray-Curtis input | CoverM Relative Abundance (%) as reported (unmapped row removed, not re-closed) | As instructed by the prompt. Re-closed (MAGs = 100%) results reported separately as a sensitivity analysis |
| Unclassified genus label | `unclassified_<deepest classified higher rank>` (e.g. `unclassified_Ancalomicrobiaceae`) | First option in the prompt; different lineages get different labels and are not merged |
| GTDB suffixed genera (`Pseudomonas_E` etc.) | Kept separate from `Pseudomonas` | They are distinct genera in GTDB |
| CLR pseudocount | read count + 1 | No zero counts (minimum 270 reads), so negligible impact |
| ALDEx2 Monte-Carlo instances | `mc.samples = 1000` (default 128) | Stable effect estimates with small n |
| ALDEx2 effect sign | + = higher in the target group (BES vs OC: + = BES higher; electrode vs tube: + = electrode higher) | ALDEx2 orders groups alphabetically, so the script checks and unifies the sign |
| ALDEx2 Wilcoxon p | Recomputed with an exact two-sided Wilcoxon test on the same MC instances (ALDEx2's original values kept in `*_ALDEx2_raw` columns) | ALDEx2 1.38 computes the lower tail as P(W < w), so only increases toward the second group (by input order) become falsely significant (see `results/REPORT.md`, section 4.1) |
| Random seed | 20261014 (PERMANOVA, betadisper, ALDEx2) | Reproducibility |
| Additional analyses | betadisper (homogeneity of dispersion), Aitchison PERMANOVA, simple-effect PERMANOVA within location/condition, `aldex.glm(~ condition + location)` | Help interpret a two-way design with n=12 |
