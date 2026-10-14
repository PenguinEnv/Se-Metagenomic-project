#!/usr/bin/env Rscript
# =============================================================================
# Se bioelectrochemical (BES) metagenome statistics - MAG level + genus level
#
# Design (n = 12): condition (OC vs BES) x location (electrode=BF vs tube=TB),
#                  3 biological replicates per combination
#
# 1. Parse coverm_all.tsv (Relative Abundance / Read Count matrices, drop unmapped)
#    + build genus-level matrices (summed) from GTDB taxonomy
# 2. two-way PERMANOVA (Bray-Curtis, adonis2 ~ condition * location, 9999 perm)
#    + betadisper (homogeneity of dispersion) + sensitivity (re-closed Bray-Curtis, Aitchison)
# 3. Ordination: CLR(count + pseudocount) PCA, Bray-Curtis PCoA
# 4. ALDEx2 (denom = "all", test = "t", effect = TRUE): BES vs OC, electrode vs tube
#    + supplementary: aldex.glm (~ condition + location)
# 5. Genus-level OC vs BES mean relative abundance comparison (table + figures)
#
# Usage (run from the statistical_analysis/ directory):
#   Rscript scripts/se_bes_metagenome_stats.R \
#       [data/coverm_all.tsv] [data/MAG_abundance_taxonomy.tsv] [results]
# =============================================================================

suppressPackageStartupMessages({
  library(vegan)
  library(ALDEx2)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

args        <- commandArgs(trailingOnly = TRUE)
coverm_file <- if (length(args) >= 1) args[1] else "data/coverm_all.tsv"
tax_file    <- if (length(args) >= 2) args[2] else "data/MAG_abundance_taxonomy.tsv"
out_dir     <- if (length(args) >= 3) args[3] else "results"

SEED       <- 20261014
N_PERM     <- 9999   # PERMANOVA / betadisper permutations
MC_SAMPLES <- as.integer(Sys.getenv("ALDEX_MC", "1000"))  # ALDEx2 Monte-Carlo instances
PSEUDO     <- 1      # CLR-PCA pseudocount (read count + 1)
ALPHA      <- 0.05   # wi.eBH threshold
TOP_GENERA <- 20     # number of top genera shown in the comparison figure

tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
mat_dir <- file.path(out_dir, "matrices")
for (d in c(tab_dir, fig_dir, mat_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Result tables use 4 significant digits (integer columns and digits = NA matrices are left as is)
write_tsv <- function(x, file, digits = 4) {
  if (!is.na(digits)) {
    rnd <- vapply(x, function(v) is.numeric(v) && any(v != round(v), na.rm = TRUE), logical(1))
    x[rnd] <- lapply(x[rnd], signif, digits)
  }
  write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
msg <- function(...) cat(sprintf(...), "\n", sep = "")

# -----------------------------------------------------------------------------
# Figure style (light surface, thin marks, recessive grid)
# -----------------------------------------------------------------------------
COL_SURFACE <- "#fcfcfb"
COL_INK     <- "#0b0b0b"
COL_INK2    <- "#52514e"
COL_MUTED   <- "#898781"
COL_GRID    <- "#e1e0d9"
COL_AXIS    <- "#c3c2b7"
COL_COND    <- c(OC = "#2a78d6", BES = "#eb6834")
LAB_COND    <- c(OC = "OC (open circuit)", BES = "BES (bioelectrochemical)")
SHP_LOC     <- c(electrode = 21, tube = 24)
LAB_LOC     <- c(electrode = "Electrode (BF)", tube = "Tube (TB)")
# 8 categorical colors for stacked bars (fixed order) + Other
COL_CAT8    <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                 "#e87ba4", "#008300", "#6250d6", "#e34948")
COL_OTHER   <- "#c3c2b7"

theme_se <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background   = element_rect(fill = COL_SURFACE, colour = NA),
      panel.background  = element_rect(fill = COL_SURFACE, colour = NA),
      panel.grid.major  = element_line(colour = COL_GRID, linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      axis.text         = element_text(colour = COL_INK2),
      axis.title        = element_text(colour = COL_INK2),
      plot.title        = element_text(colour = COL_INK, face = "bold"),
      plot.subtitle     = element_text(colour = COL_INK2, size = rel(0.85)),
      plot.caption      = element_text(colour = COL_MUTED, size = rel(0.75), hjust = 0),
      legend.position   = "bottom",
      legend.title      = element_blank(),
      legend.text       = element_text(colour = COL_INK2),
      strip.text        = element_text(colour = COL_INK, face = "bold")
    )
}

save_fig <- function(p, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = width, height = height,
         bg = COL_SURFACE, device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = width, height = height,
         dpi = 300, bg = COL_SURFACE, type = "cairo")
}

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))
p_txt <- function(p) ifelse(p < 0.001, "p<0.001", paste0("p=", formatC(p, format = "f", digits = 3)))

# =============================================================================
# 1. Data parsing
# =============================================================================
msg("== 1. Parsing input ==")
cov <- read.delim(coverm_file, check.names = FALSE, stringsAsFactors = FALSE)
cov <- cov[cov$Genome != "unmapped", , drop = FALSE]
rownames(cov) <- cov$Genome

# "OCBF01_R1.clean.fq.gz Relative Abundance (%)" -> "OCBF01"
clean_sample <- function(x) sub("_R[12].*$", "", sub(" .*$", "", x))
extract_metric <- function(suffix) {
  cols <- names(cov)[endsWith(names(cov), paste0(" ", suffix))]
  m <- as.matrix(cov[, cols, drop = FALSE])
  colnames(m) <- clean_sample(cols)
  storage.mode(m) <- "double"
  m
}
ra_mag <- extract_metric("Relative Abundance (%)")   # MAG x sample, % of all reads
rc_mag <- extract_metric("Read Count")               # MAG x sample, mapped reads
stopifnot(identical(colnames(ra_mag), colnames(rc_mag)), !anyNA(rc_mag))
rc_mag <- round(rc_mag)
samples <- colnames(ra_mag)

parts <- regmatches(samples, regexec("^(OC|BES)(BF|TB)(\\d+)$", samples))
if (any(lengths(parts) != 4)) stop("Unexpected sample name: ",
                                   paste(samples[lengths(parts) != 4], collapse = ", "))
meta <- data.frame(
  sample    = samples,
  condition = factor(sapply(parts, `[`, 2), levels = c("OC", "BES")),
  location  = factor(ifelse(sapply(parts, `[`, 3) == "BF", "electrode", "tube"),
                     levels = c("electrode", "tube")),
  replicate = sapply(parts, `[`, 4),
  stringsAsFactors = FALSE
)
meta$group <- interaction(meta$condition, meta$location, sep = "_")
rownames(meta) <- meta$sample
msg("Samples: %d | MAGs: %d", ncol(ra_mag), nrow(ra_mag))
print(table(meta$condition, meta$location))
write_tsv(cbind(meta,
                MAG_relabund_total_pct = round(colSums(ra_mag), 3),
                unmapped_pct = round(100 - colSums(ra_mag), 3),
                MAG_read_count_total = colSums(rc_mag)),
          file.path(tab_dir, "00_sample_metadata.tsv"))

# --- taxonomy ---------------------------------------------------------------
tax <- read.delim(tax_file, check.names = FALSE, stringsAsFactors = FALSE)
names(tax)[1] <- "MAG"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
stopifnot(all(ranks %in% names(tax)))
tax <- tax[, c("MAG", ranks)]
rownames(tax) <- tax$MAG
missing <- setdiff(rownames(ra_mag), tax$MAG)
if (length(missing)) stop("MAGs without taxonomy: ", paste(missing, collapse = ", "))
tax <- tax[rownames(ra_mag), ]

is_uncl <- function(x) is.na(x) | trimws(x) == "" | tolower(trimws(x)) %in% c("unclassified", "na")

# Genus label: if unclassified, use the deepest classified rank -> "unclassified_<higher taxon>"
tax$genus_label <- vapply(seq_len(nrow(tax)), function(i) {
  if (!is_uncl(tax$genus[i])) return(tax$genus[i])
  for (r in c("family", "order", "class", "phylum", "domain")) {
    if (!is_uncl(tax[i, r])) return(paste0("unclassified_", tax[i, r]))
  }
  "unclassified"
}, character(1))
# MAG display name: species > "Genus sp." > genus_label
tax$taxon_label <- ifelse(!is_uncl(tax$species), tax$species,
                          ifelse(!is_uncl(tax$genus), paste(tax$genus, "sp."), tax$genus_label))
write_tsv(tax, file.path(tab_dir, "00_MAG_taxonomy_genus_map.tsv"))

# --- genus-level matrices: sum MAGs sharing a genus label --------------------
ra_gen <- rowsum(ra_mag, tax$genus_label)
rc_gen <- rowsum(rc_mag, tax$genus_label)
gen_tax <- do.call(rbind, lapply(split(tax, tax$genus_label), function(d) {
  # higher rank of a genus label: the unique value, or "mixed" if not unique
  one <- function(v) if (length(unique(v)) == 1) v[1] else "mixed"
  data.frame(genus = d$genus_label[1], domain = one(d$domain), phylum = one(d$phylum),
             class = one(d$class), order = one(d$order), family = one(d$family),
             n_MAGs = nrow(d), MAGs = paste(d$MAG, collapse = ";"),
             stringsAsFactors = FALSE)
}))
rownames(gen_tax) <- gen_tax$genus
gen_tax <- gen_tax[rownames(ra_gen), ]
msg("Genus-level features: %d (of which %d 'unclassified_*' labels)",
    nrow(ra_gen), sum(grepl("^unclassified", rownames(ra_gen))))

write_tsv(data.frame(MAG = rownames(ra_mag), ra_mag, check.names = FALSE),
          file.path(mat_dir, "MAG_relabund.tsv"), digits = NA)
write_tsv(data.frame(MAG = rownames(rc_mag), rc_mag, check.names = FALSE),
          file.path(mat_dir, "MAG_readcount.tsv"), digits = NA)
write_tsv(data.frame(genus = rownames(ra_gen), ra_gen, check.names = FALSE),
          file.path(mat_dir, "genus_relabund.tsv"), digits = NA)
write_tsv(data.frame(genus = rownames(rc_gen), rc_gen, check.names = FALSE),
          file.path(mat_dir, "genus_readcount.tsv"), digits = NA)
write_tsv(gen_tax, file.path(tab_dir, "00_genus_taxonomy.tsv"))

levels_list <- list(
  MAG   = list(ra = ra_mag, rc = rc_mag),
  genus = list(ra = ra_gen, rc = rc_gen)
)

clr_transform <- function(counts, pseudo = PSEUDO) {
  l <- log(counts + pseudo)
  sweep(l, 2, colMeans(l))          # feature x sample, centred on each sample's log geometric mean
}

# =============================================================================
# 2. two-way PERMANOVA
# =============================================================================
msg("== 2. PERMANOVA ==")
tidy_adonis <- function(a, level, distance, data_used) {
  out <- data.frame(level = level, distance = distance, data = data_used,
                    term = rownames(a), Df = a$Df, SumOfSqs = a$SumOfSqs,
                    R2 = a$R2, F = a$F, p_value = a$`Pr(>F)`, row.names = NULL)
  out$sig <- ifelse(is.na(out$p_value), "",
                    ifelse(out$p_value < 0.001, "***", ifelse(out$p_value < 0.01, "**",
                    ifelse(out$p_value < 0.05, "*", ifelse(out$p_value < 0.1, ".", "ns")))))
  out
}
run_adonis <- function(d, level, distance, data_used) {
  set.seed(SEED)
  a <- adonis2(d ~ condition * location, data = meta, permutations = N_PERM, by = "terms")
  tidy_adonis(as.data.frame(a), level, distance, data_used)
}

dists <- list()
perm_main <- list(); perm_sens <- list(); disp <- list(); simple <- list()
for (lv in names(levels_list)) {
  ra <- levels_list[[lv]]$ra; rc <- levels_list[[lv]]$rc

  bc       <- vegdist(t(ra), method = "bray")                       # main analysis
  bc_close <- vegdist(t(sweep(ra, 2, colSums(ra), "/")), "bray")      # re-closed so MAGs sum to 100%
  ait      <- dist(t(clr_transform(rc)))                              # Aitchison (Euclidean on CLR)
  dists[[lv]] <- list(bray = bc, aitchison = ait)

  perm_main[[lv]] <- run_adonis(bc, lv, "Bray-Curtis", "CoverM relative abundance (%)")
  perm_sens[[lv]] <- rbind(
    perm_main[[lv]],
    run_adonis(bc_close, lv, "Bray-Curtis", "relative abundance re-closed to mapped MAGs = 100%"),
    run_adonis(ait, lv, "Aitchison", sprintf("CLR(read count + %g)", PSEUDO)))

  # Homogeneity of dispersion (is PERMANOVA significance a location or a dispersion effect?)
  for (f in c("condition", "location", "group")) {
    set.seed(SEED)
    bd <- betadisper(bc, meta[[f]])
    pt <- permutest(bd, permutations = N_PERM)
    disp[[paste(lv, f)]] <- data.frame(
      level = lv, distance = "Bray-Curtis", factor = f,
      mean_distance_to_centroid = paste(sprintf("%s=%.3f", names(bd$group.distances),
                                                bd$group.distances), collapse = "; "),
      F = pt$tab$F[1], p_value = pt$tab$`Pr(>F)`[1])
  }

  # Simple effects (exploratory): 3 vs 3 has few distinct permutations, minimum p = 0.1
  for (loc in levels(meta$location)) {
    idx <- meta$location == loc
    set.seed(SEED)
    a <- adonis2(as.dist(as.matrix(bc)[idx, idx]) ~ condition, data = meta[idx, ],
                 permutations = N_PERM)
    simple[[paste(lv, loc)]] <- data.frame(level = lv, test = paste0("BES vs OC within ", loc),
                                           n = sum(idx), R2 = a$R2[1], F = a$F[1],
                                           p_value = a$`Pr(>F)`[1])
  }
  for (cn in levels(meta$condition)) {
    idx <- meta$condition == cn
    set.seed(SEED)
    a <- adonis2(as.dist(as.matrix(bc)[idx, idx]) ~ location, data = meta[idx, ],
                 permutations = N_PERM)
    simple[[paste(lv, cn)]] <- data.frame(level = lv, test = paste0("electrode vs tube within ", cn),
                                          n = sum(idx), R2 = a$R2[1], F = a$F[1],
                                          p_value = a$`Pr(>F)`[1])
  }
}
perm_main_df <- do.call(rbind, perm_main)
write_tsv(perm_main[["MAG"]],   file.path(tab_dir, "01_PERMANOVA_MAG.tsv"))
write_tsv(perm_main[["genus"]], file.path(tab_dir, "01_PERMANOVA_genus.tsv"))
write_tsv(perm_main_df,         file.path(tab_dir, "01_PERMANOVA_MAG_vs_genus.tsv"))
write_tsv(do.call(rbind, perm_sens), file.path(tab_dir, "01b_PERMANOVA_sensitivity.tsv"))
write_tsv(do.call(rbind, disp),      file.path(tab_dir, "01c_betadisper.tsv"))
write_tsv(do.call(rbind, simple),    file.path(tab_dir, "01d_PERMANOVA_simple_effects.tsv"))
print(perm_main_df[, c("level", "term", "Df", "R2", "F", "p_value", "sig")], row.names = FALSE)

# =============================================================================
# 3. Ordination
# =============================================================================
msg("== 3. Ordination ==")
perm_subtitle <- function(tab) {
  tab <- tab[tab$term %in% c("condition", "location", "condition:location"), ]
  nm <- c(condition = "Condition", location = "Location", "condition:location" = "Interaction")
  paste0("PERMANOVA (", tab$distance[1], ", ", N_PERM, " perm)\n",
         paste(sprintf("%s R2=%.2f, %s", nm[tab$term], tab$R2, p_txt(tab$p_value)),
               collapse = " | "))
}
plot_ord <- function(scores, xlab, ylab, title, subtitle) {
  df <- cbind(meta, x = scores[meta$sample, 1], y = scores[meta$sample, 2])
  ggplot(df, aes(x, y)) +
    geom_hline(yintercept = 0, colour = COL_AXIS, linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = COL_AXIS, linewidth = 0.3) +
    geom_point(aes(fill = condition, shape = location), size = 3.8,
               colour = COL_SURFACE, stroke = 0.9) +
    geom_text_repel(aes(label = sample), size = 3, colour = COL_INK2,
                    box.padding = 0.45, point.padding = 0.3, min.segment.length = 0.2,
                    segment.colour = COL_AXIS, seed = 1, max.overlaps = Inf) +
    scale_fill_manual(values = COL_COND, labels = LAB_COND) +
    scale_shape_manual(values = SHP_LOC, labels = LAB_LOC) +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)),
           shape = guide_legend(override.aes = list(fill = COL_MUTED, size = 4))) +
    labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    theme_se()
}

ord_plots <- list(); ord_scores <- list()
for (lv in names(levels_list)) {
  rc <- levels_list[[lv]]$rc
  lv_name <- if (lv == "MAG") "MAG level" else "Genus level"

  # CLR-PCA
  clr <- clr_transform(rc)
  pca <- prcomp(t(clr), center = TRUE, scale. = FALSE)
  ve  <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  p_pca <- plot_ord(pca$x[, 1:2],
                    sprintf("PC1 (%.1f%%)", ve[1]), sprintf("PC2 (%.1f%%)", ve[2]),
                    sprintf("%s - CLR PCA (read count + %g)", lv_name, PSEUDO),
                    perm_subtitle(perm_sens[[lv]][perm_sens[[lv]]$distance == "Aitchison", ]))

  # Bray-Curtis PCoA
  pco <- cmdscale(dists[[lv]]$bray, k = nrow(meta) - 1, eig = TRUE)
  pe  <- 100 * pco$eig / sum(pco$eig[pco$eig > 0])
  p_pco <- plot_ord(pco$points[, 1:2],
                    sprintf("PCoA1 (%.1f%%)", pe[1]), sprintf("PCoA2 (%.1f%%)", pe[2]),
                    sprintf("%s - Bray-Curtis PCoA (relative abundance)", lv_name),
                    perm_subtitle(perm_main[[lv]]))

  tag <- if (lv == "MAG") "MAG" else "genus"
  save_fig(p_pca, sprintf("03_%s_CLR_PCA", tag), 7, 6)
  save_fig(p_pco, sprintf("03_%s_BrayCurtis_PCoA", tag), 7, 6)
  ord_plots[[lv]] <- list(pca = p_pca, pco = p_pco)
  ord_scores[[lv]] <- rbind(
    data.frame(level = lv, method = "CLR_PCA", sample = rownames(pca$x),
               axis1 = pca$x[, 1], axis2 = pca$x[, 2],
               axis1_pct = ve[1], axis2_pct = ve[2]),
    data.frame(level = lv, method = "BrayCurtis_PCoA", sample = rownames(pco$points),
               axis1 = pco$points[, 1], axis2 = pco$points[, 2],
               axis1_pct = pe[1], axis2_pct = pe[2]))
}
p_all <- (ord_plots$MAG$pca + ord_plots$MAG$pco) / (ord_plots$genus$pca + ord_plots$genus$pco) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_fig(p_all, "03_ordination_panel", 14, 12)
write_tsv(do.call(rbind, ord_scores), file.path(tab_dir, "03_ordination_scores.tsv"))

# =============================================================================
# 4. ALDEx2 differential abundance
# =============================================================================
msg("== 4. ALDEx2 ==")
# Positive effect = higher in the target group
comparisons <- list(
  BES_vs_OC         = list(factor = "condition", ref = "OC",   target = "BES"),
  electrode_vs_tube = list(factor = "location",  ref = "tube", target = "electrode")
)

# Correction of ALDEx2 1.38 expected Wilcoxon p-values
# ALDEx2::wilcox.fast computes only the upper tail P(W >= w), and aldex.ttest derives
# the other tail as 1 - P(W >= w) = P(W < w). For a discrete statistic the correct
# tail is P(W <= w), so that tail is too small by P(W = w). Example: in 6 vs 6, if the
# second group is entirely higher the two-sided p = 0 (exact minimum is 2/924 = 0.0022),
# so one direction (depending on group order) becomes falsely significant. Here both
# tails are computed from the exact Wilcoxon distribution on the same Monte-Carlo CLR
# instances, and the rest of the aggregation matches ALDEx2
# (per-direction 2*p, BH per MC instance, mean over MC, then the smaller direction).
wilcox_exact_mc <- function(clr_obj, conds) {
  mc   <- getMonteCarloInstances(clr_obj)            # per-sample feature x MC matrices
  n_mc <- numMCInstances(clr_obj)
  g1 <- conds == unique(conds)[1]
  n1 <- sum(g1); n2 <- sum(!g1)
  nf <- length(getFeatureNames(clr_obj))
  P_g <- P_l <- B_g <- B_l <- matrix(NA_real_, nf, n_mc)
  for (k in seq_len(n_mc)) {
    v <- vapply(mc, function(m) m[, k], numeric(nf))   # feature x sample
    r <- t(apply(v, 1, rank))
    W <- rowSums(r[, g1, drop = FALSE]) - n1 * (n1 + 1) / 2
    ties <- apply(v, 1, function(z) anyDuplicated(z) > 0)
    pg <- pwilcox(W - 1, n1, n2, lower.tail = FALSE)  # P(W >= w)
    pl <- pwilcox(W, n1, n2)                          # P(W <= w)
    for (i in which(ties)) {                          # normal approximation if ties
      pg[i] <- wilcox.test(v[i, g1], v[i, !g1], alternative = "greater", exact = FALSE)$p.value
      pl[i] <- wilcox.test(v[i, g1], v[i, !g1], alternative = "less", exact = FALSE)$p.value
    }
    P_g[, k] <- pmin(1, 2 * pg); P_l[, k] <- pmin(1, 2 * pl)
    B_g[, k] <- p.adjust(P_g[, k], "BH"); B_l[, k] <- p.adjust(P_l[, k], "BH")
  }
  data.frame(wi.ep  = pmin(rowMeans(P_g), rowMeans(P_l)),
             wi.eBH = pmin(rowMeans(B_g), rowMeans(B_l)))
}

run_aldex_t <- function(counts, cmp) {
  # ALDEx2 takes a character vector and orders levels alphabetically -> sign checked below
  conds <- as.character(meta[colnames(counts), cmp$factor])
  # Same steps as aldex(test = "t", effect = TRUE, denom = "all"), run separately
  # because the Wilcoxon correction needs the MC instances
  set.seed(SEED)
  x  <- aldex.clr(counts, conds, mc.samples = MC_SAMPLES, denom = "all", verbose = FALSE)
  tt <- aldex.ttest(x, paired.test = FALSE, verbose = FALSE)
  ef <- aldex.effect(x, include.sample.summary = FALSE, verbose = FALSE)
  wi <- wilcox_exact_mc(x, conds)
  res <- data.frame(ef, tt[rownames(ef), ], check.names = FALSE)
  rab_t <- res[[paste0("rab.win.", cmp$target)]]
  rab_r <- res[[paste0("rab.win.", cmp$ref)]]
  # Check effect sign: make positive = higher in target
  agree <- mean(sign(res$diff.btw) == sign(rab_t - rab_r), na.rm = TRUE)
  msg("    sign agreement diff.btw vs (rab.win.%s - rab.win.%s): %.2f", cmp$target, cmp$ref, agree)
  if (agree < 0.5) {
    msg("    -> ALDEx2 effect sign flipped so that + = higher in %s", cmp$target)
    res$diff.btw <- -res$diff.btw
    res$effect   <- -res$effect
  }
  msg("    wi.eBH < %.2f: corrected %d | ALDEx2 raw %d", ALPHA,
      sum(wi$wi.eBH < ALPHA), sum(res$wi.eBH < ALPHA))
  data.frame(feature = rownames(res),
             rab.all = res$rab.all, rab.win.ref = rab_r, rab.win.target = rab_t,
             diff.btw = res$diff.btw, diff.win = res$diff.win,
             effect = res$effect, overlap = res$overlap,
             we.ep = res$we.ep, we.eBH = res$we.eBH,
             wi.ep = wi$wi.ep, wi.eBH = wi$wi.eBH,
             wi.ep_ALDEx2_raw = res$wi.ep, wi.eBH_ALDEx2_raw = res$wi.eBH,
             stringsAsFactors = FALSE)
}

aldex_res <- list()
aldex_summary <- list()
for (lv in names(levels_list)) {
  rc  <- levels_list[[lv]]$rc
  tag <- if (lv == "MAG") "MAG" else "genus"
  ra  <- levels_list[[lv]]$ra
  for (cn in names(comparisons)) {
    cmp <- comparisons[[cn]]
    msg("  %s: %s", lv, cn)
    r <- run_aldex_t(rc, cmp)
    names(r)[names(r) == "rab.win.ref"]    <- paste0("rab.win.", cmp$ref)
    names(r)[names(r) == "rab.win.target"] <- paste0("rab.win.", cmp$target)

    # Add group mean relative abundance (%)
    g <- meta[colnames(ra), cmp$factor]
    r[[paste0("mean_relabund_", cmp$ref)]]    <- rowMeans(ra[r$feature, g == cmp$ref, drop = FALSE])
    r[[paste0("mean_relabund_", cmp$target)]] <- rowMeans(ra[r$feature, g == cmp$target, drop = FALSE])
    r$enriched_in <- ifelse(r$effect > 0, cmp$target, cmp$ref)
    r$significant_wi.eBH <- r$wi.eBH < ALPHA
    r$abs_effect_ge_1    <- abs(r$effect) >= 1

    if (lv == "MAG") {
      tx <- tax[r$feature, c("taxon_label", "species", "genus", "genus_label",
                             "family", "order", "class", "phylum")]
      r <- cbind(MAG = r$feature, tx, r[, -1], row.names = NULL)
    } else {
      gt <- gen_tax[r$feature, c("family", "order", "class", "phylum", "n_MAGs")]
      r <- cbind(genus = r$feature, gt, r[, -1], row.names = NULL)
    }
    r <- r[order(-r$effect), ]
    aldex_res[[paste(lv, cn)]] <- r

    write_tsv(r, file.path(tab_dir, sprintf("04_ALDEx2_%s_%s_all.tsv", tag, cn)))
    sig <- r[r$significant_wi.eBH, ]
    write_tsv(sig, file.path(tab_dir, sprintf("04_ALDEx2_%s_%s_significant.tsv", tag, cn)))

    aldex_summary[[paste(lv, cn)]] <- data.frame(
      level = lv, comparison = cn, n_features = nrow(r),
      n_sig_wi.eBH = nrow(sig),
      n_sig_up_target = sum(sig$effect > 0), n_sig_up_ref = sum(sig$effect < 0),
      n_sig_we.eBH = sum(r$we.eBH < ALPHA),
      n_abs_effect_ge_1 = sum(r$abs_effect_ge_1),
      n_sig_and_abs_effect_ge_1 = sum(r$significant_wi.eBH & r$abs_effect_ge_1),
      target = cmp$target, ref = cmp$ref)
  }
}
aldex_summary_df <- do.call(rbind, aldex_summary)
write_tsv(aldex_summary_df, file.path(tab_dir, "04_ALDEx2_summary.tsv"))
print(aldex_summary_df, row.names = FALSE)

# --- Supplementary: aldex.glm (~ condition + location), each factor adjusted for the other
msg("== 4b. ALDEx2 glm (~ condition + location) ==")
mm <- model.matrix(~ condition + location, data = meta)
glm_res <- list()
for (lv in names(levels_list)) {
  set.seed(SEED)
  x  <- aldex.clr(levels_list[[lv]]$rc, mm, mc.samples = MC_SAMPLES, denom = "all",
                  verbose = FALSE)
  gl <- aldex.glm(x, fdr.method = "BH")   # model matrix is passed via aldex.clr(conds = mm)
  pick <- function(col) {
    if (!col %in% names(gl)) stop("aldex.glm column not found: ", col)
    gl[[col]]
  }
  out <- data.frame(
    feature = rownames(gl),
    conditionBES_estimate = pick("conditionBES:Est"),
    conditionBES_p        = pick("conditionBES:pval"),
    conditionBES_pBH      = pick("conditionBES:pval.padj"),
    locationtube_estimate = pick("locationtube:Est"),
    locationtube_p        = pick("locationtube:pval"),
    locationtube_pBH      = pick("locationtube:pval.padj"),
    stringsAsFactors = FALSE)
  names(out)[1] <- if (lv == "MAG") "MAG" else "genus"
  if (lv == "MAG") out <- cbind(out[1], taxon_label = tax[out$MAG, "taxon_label"], out[-1])
  out <- out[order(out$conditionBES_pBH, -abs(out$conditionBES_estimate)), ]
  glm_res[[lv]] <- out
  write_tsv(out, file.path(tab_dir, sprintf("04b_ALDEx2_glm_%s.tsv",
                                            if (lv == "MAG") "MAG" else "genus")))
  msg("  %s: condition (adj. for location) BH<0.05: %d | location (adj. for condition) BH<0.05: %d",
      lv, sum(out$conditionBES_pBH < ALPHA), sum(out$locationtube_pBH < ALPHA))
}

# --- ALDEx2 effect plot (BES vs OC), MAG & genus ----------------------------
MARK_CAPTION <- paste("* t-test model: wi.eBH < 0.05 (Wilcoxon, exact).",
                      "+ glm ~ condition + location: condition BH < 0.05.")
plot_effect <- function(r, id_col, label_col, title, glm_tab, n_show = 30) {
  glm_sig <- glm_tab[[id_col]][glm_tab$conditionBES_pBH < ALPHA]
  r$glm_sig <- r[[id_col]] %in% glm_sig
  # Show significant features (t-test or glm) first, then fill by largest |effect|
  r <- r[order(-(r$significant_wi.eBH | r$glm_sig), -abs(r$effect)), ]
  show <- r[seq_len(min(max(n_show, sum(r$significant_wi.eBH | r$glm_sig)), nrow(r))), ]
  show$lab <- paste0(show[[label_col]],
                     if (label_col != id_col) paste0("  [", show[[id_col]], "]") else "",
                     ifelse(show$significant_wi.eBH, "  *", ""),
                     ifelse(show$glm_sig, "  +", ""))
  show <- show[order(show$effect), ]
  show$lab <- factor(show$lab, levels = show$lab)
  ggplot(show, aes(effect, lab, fill = enriched_in)) +
    geom_vline(xintercept = 0, colour = COL_AXIS, linewidth = 0.4) +
    geom_vline(xintercept = c(-1, 1), colour = COL_GRID, linewidth = 0.4, linetype = "dashed") +
    geom_col(width = 0.7) +
    scale_fill_manual(values = COL_COND, labels = c(OC = "Higher in OC", BES = "Higher in BES")) +
    labs(x = "ALDEx2 effect size (BES vs OC; + = higher in BES)", y = NULL, title = title,
         subtitle = sprintf("Shown: significant features plus the largest |effect| (%d in total)", nrow(show)),
         caption = paste(MARK_CAPTION, "Dashed lines: |effect| = 1.")) +
    theme_se(10) + theme(panel.grid.major.y = element_blank())
}
p_eff_mag <- plot_effect(aldex_res[["MAG BES_vs_OC"]], "MAG", "taxon_label",
                         "MAG level - ALDEx2 BES vs OC", glm_res[["MAG"]])
p_eff_gen <- plot_effect(aldex_res[["genus BES_vs_OC"]], "genus", "genus",
                         "Genus level - ALDEx2 BES vs OC", glm_res[["genus"]])
save_fig(p_eff_mag, "04_MAG_ALDEx2_BES_vs_OC_effect", 10.5, 9.5)
save_fig(p_eff_gen, "04_genus_ALDEx2_BES_vs_OC_effect", 9, 8)

# =============================================================================
# 5. Genus-level OC vs BES summary
# =============================================================================
msg("== 5. Genus-level OC vs BES summary ==")
grp_mean <- function(m, idx) rowMeans(m[, idx, drop = FALSE])
grp_sd   <- function(m, idx) apply(m[, idx, drop = FALSE], 1, sd)
is_oc  <- meta$condition == "OC";  is_bes <- meta$condition == "BES"
is_bf  <- meta$location == "electrode"; is_tb <- meta$location == "tube"
ag  <- aldex_res[["genus BES_vs_OC"]];  rownames(ag) <- ag$genus
agl <- glm_res[["genus"]];              rownames(agl) <- agl$genus

gsum <- data.frame(
  genus = rownames(ra_gen),
  gen_tax[rownames(ra_gen), c("phylum", "class", "family", "n_MAGs")],
  mean_all = grp_mean(ra_gen, TRUE),
  mean_OC = grp_mean(ra_gen, is_oc),   sd_OC = grp_sd(ra_gen, is_oc),
  mean_BES = grp_mean(ra_gen, is_bes), sd_BES = grp_sd(ra_gen, is_bes),
  mean_OC_electrode  = grp_mean(ra_gen, is_oc & is_bf),
  mean_BES_electrode = grp_mean(ra_gen, is_bes & is_bf),
  mean_OC_tube       = grp_mean(ra_gen, is_oc & is_tb),
  mean_BES_tube      = grp_mean(ra_gen, is_bes & is_tb),
  stringsAsFactors = FALSE, row.names = NULL)
gsum$diff_BES_minus_OC <- gsum$mean_BES - gsum$mean_OC
gsum$log2FC_BES_OC     <- log2((gsum$mean_BES + 0.01) / (gsum$mean_OC + 0.01))  # +0.01% offset
gsum$ALDEx2_effect     <- ag[gsum$genus, "effect"]
gsum$ALDEx2_wi.eBH     <- ag[gsum$genus, "wi.eBH"]
gsum$ALDEx2_we.eBH     <- ag[gsum$genus, "we.eBH"]
gsum$ALDEx2_sig        <- gsum$ALDEx2_wi.eBH < ALPHA
gsum$glm_condition_pBH <- agl[gsum$genus, "conditionBES_pBH"]
gsum$dominant_in       <- ifelse(gsum$diff_BES_minus_OC > 0, "BES", "OC")
gsum <- gsum[order(-gsum$diff_BES_minus_OC), ]
write_tsv(gsum, file.path(tab_dir, "05_genus_OC_vs_BES_summary.tsv"))

# --- Figure 5a: top genera mean abundance OC vs BES (+ replicate points, ALDEx2 effect)
top <- head(gsum[order(-gsum$mean_all), ], TOP_GENERA)
top <- top[order(top$diff_BES_minus_OC), ]
top$lab <- paste0(top$genus, ifelse(top$ALDEx2_sig, " *", ""),
                  ifelse(top$glm_condition_pBH < ALPHA, " +", ""))
lab_levels <- top$lab

SCOPES <- c(all = "All samples (6 vs 6)", electrode = "Electrode, BF (3 vs 3)",
            tube = "Tube, TB (3 vs 3)")
scope_idx <- function(sc) if (sc == "all") rep(TRUE, nrow(meta)) else meta$location == sc
bar_df <- do.call(rbind, lapply(names(SCOPES), function(sc) {
  do.call(rbind, lapply(c("OC", "BES"), function(cn) {
    j <- scope_idx(sc) & meta$condition == cn
    data.frame(lab = top$lab, scope = SCOPES[[sc]], condition = cn,
               mean = rowMeans(ra_gen[top$genus, j, drop = FALSE]))
  }))
}))
pt_df <- do.call(rbind, lapply(names(SCOPES), function(sc) {
  idx <- scope_idx(sc)
  do.call(rbind, lapply(seq_len(nrow(top)), function(i) {
    data.frame(lab = top$lab[i], scope = SCOPES[[sc]], sample = meta$sample[idx],
               condition = as.character(meta$condition[idx]), location = meta$location[idx],
               value = ra_gen[top$genus[i], meta$sample[idx]])
  }))
}))
fix_levels <- function(d) {
  d$lab <- factor(d$lab, levels = lab_levels)
  d$scope <- factor(d$scope, levels = SCOPES)
  d$condition <- factor(d$condition, levels = c("OC", "BES"))
  d
}
bar_df <- fix_levels(bar_df); pt_df <- fix_levels(pt_df)

dodge <- position_dodge(width = 0.8)
p_bar <- ggplot(bar_df, aes(mean, lab, fill = condition)) +
  geom_col(position = dodge, width = 0.75, colour = COL_SURFACE, linewidth = 0.4) +
  geom_point(data = pt_df, aes(value, lab, group = condition, shape = location),
             position = position_jitterdodge(jitter.height = 0.1, jitter.width = 0,
                                             dodge.width = 0.8, seed = 1),
             size = 1.4, fill = COL_SURFACE, colour = COL_INK2, stroke = 0.35,
             inherit.aes = FALSE) +
  facet_grid(~ scope) +
  scale_fill_manual(values = COL_COND, labels = LAB_COND) +
  scale_shape_manual(values = SHP_LOC, labels = LAB_LOC) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Relative abundance (% of all reads)", y = NULL,
       title = sprintf("Top %d genera by mean abundance - OC vs BES", TOP_GENERA),
       subtitle = "Bars = group mean, points = replicates. Rows ordered by overall BES - OC difference (BES-dominant at top).") +
  theme_se(10) + theme(panel.grid.major.y = element_blank(),
                       panel.spacing.x = unit(1.2, "lines"))

eff_df <- data.frame(lab = factor(top$lab, levels = lab_levels), effect = top$ALDEx2_effect,
                     dir = factor(ifelse(top$ALDEx2_effect > 0, "BES", "OC"), levels = c("OC", "BES")),
                     scope = "ALDEx2, all samples")
p_eff <- ggplot(eff_df, aes(effect, lab, fill = dir)) +
  geom_vline(xintercept = 0, colour = COL_AXIS, linewidth = 0.4) +
  geom_vline(xintercept = c(-1, 1), colour = COL_GRID, linewidth = 0.4, linetype = "dashed") +
  geom_col(width = 0.6) +
  facet_grid(~ scope) +
  scale_fill_manual(values = COL_COND, guide = "none") +
  labs(x = "Effect (+ = higher in BES)", y = NULL) +
  theme_se(10) + theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())

p_genus <- (p_bar + p_eff + plot_layout(widths = c(3, 0.6), guides = "collect") +
              plot_annotation(caption = paste("Genus label marks:", MARK_CAPTION),
                              theme = theme_se(10))) &
  theme(legend.position = "bottom")
save_fig(p_genus, "05_genus_OC_vs_BES_barplot", 16, 9)

# --- Figure 5b: per-sample genus composition (top 8 + Other) ---------------------
top8 <- head(gsum$genus[order(-gsum$mean_all)], 8)
comp <- do.call(rbind, lapply(meta$sample, function(s) {
  v <- ra_gen[, s]
  data.frame(sample = s, genus = c(top8, "Other MAGs"),
             value = c(v[top8], sum(v[!names(v) %in% top8])))
}))
comp$genus  <- factor(comp$genus, levels = rev(c(top8, "Other MAGs")))
comp$sample <- factor(comp$sample, levels = meta$sample)
comp <- merge(comp, meta[, c("sample", "condition", "location")], by = "sample")
comp$panel <- factor(paste(comp$condition, LAB_LOC[as.character(comp$location)], sep = " / "),
                     levels = c(outer(c("OC", "BES"), LAB_LOC, paste, sep = " / ")))
p_comp <- ggplot(comp, aes(sample, value, fill = genus)) +
  geom_col(width = 0.8, colour = COL_SURFACE, linewidth = 0.4) +
  facet_grid(~ panel, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = setNames(c(COL_CAT8, COL_OTHER), c(top8, "Other MAGs")),
                    breaks = c(top8, "Other MAGs")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Relative abundance (% of all reads)",
       title = "Genus composition per sample (top 8 genera)",
       caption = "Bars sum to the MAG-mapped fraction; the remainder up to 100% is unmapped reads.") +
  theme_se(10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank()) +
  guides(fill = guide_legend(nrow = 3))
save_fig(p_comp, "05_genus_composition_stacked", 11, 6.5)

# =============================================================================
# Session info
# =============================================================================
writeLines(c(capture.output(sessionInfo()),
             "", sprintf("SEED=%d N_PERM=%d MC_SAMPLES=%d PSEUDO=%g ALPHA=%g",
                         SEED, N_PERM, MC_SAMPLES, PSEUDO, ALPHA)),
           file.path(out_dir, "sessionInfo.txt"))
msg("Done. Outputs in: %s", normalizePath(out_dir))
