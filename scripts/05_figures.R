# Step 5 — Summary figures: PCA, volcano plot and heatmap of significant DMPs
# (with sepsis samples shown alongside, to see whether DMPs also change with
# systemic inflammation).

source("scripts/_common.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
})

d     <- readRDS(results_path("rds", "02_normalized.rds"))
dmp   <- readRDS(results_path("rds", "04_dmp_all.rds"))$dmp
pheno <- d$pheno

group_levels  <- c("Healthy", "Cancer", "Sepsis")
group_colours <- c(Healthy = "#2E86AB", Cancer = "#E4572E", Sepsis = "#7D5BA6")
pheno$Sample_Group <- factor(pheno$Sample_Group, levels = group_levels)
pheno$Batch <- ifelse(pheno$Shared_Chip, "Chip shared by healthy + cancer", "Single-group chip")

# ---- PCA on the most variable CpGs --------------------------------------------
top_var <- order(apply(d$M, 1, var), decreasing = TRUE)[1:min(10000, nrow(d$M))]
pca     <- prcomp(t(d$M[top_var, ]), scale. = FALSE)
var_pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_df <- data.frame(pca$x[, 1:2], Group = pheno$Sample_Group, Batch = pheno$Batch)
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = Group, shape = Batch)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_colour_manual(values = group_colours) +
  labs(x = paste0("PC1 (", var_pct[1], "%)"), y = paste0("PC2 (", var_pct[2], "%)"),
       colour = NULL, shape = NULL, title = "PCA — 10,000 most variable CpGs (M-values)") +
  theme_minimal() +
  theme(legend.position = "right")
ggsave(results_path("figures", "05_pca.png"), p_pca, width = 8, height = 5, dpi = 150)

# ---- Volcano plot --------------------------------------------------------------
dmp$status <- ifelse(dmp$significant, dmp$direction, "not significant")
p_volcano <- ggplot(dmp, aes(delta_beta, -log10(p_value), colour = status)) +
  geom_point(size = 0.4, alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1) * cfg$dmp$min_abs_delta_beta, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(hyper = "#E4572E", hypo = "#2E86AB", "not significant" = "grey75")) +
  labs(x = paste0("Delta beta (", cfg$dmp$contrast, ")"), y = expression(-log[10](p)), colour = NULL,
       title = paste0(sum(dmp$significant), " DMPs (FDR < ", cfg$dmp$fdr_max,
                      ", |delta beta| >= ", cfg$dmp$min_abs_delta_beta, ")")) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_minimal()
ggsave(results_path("figures", "05_volcano.png"), p_volcano, width = 7, height = 5, dpi = 150)

# ---- Heatmap of the top significant DMPs, with sepsis samples alongside ----------
# Clustering all DMPs needs an n x n distance matrix (44k CpGs ~ 15 GB), so the
# heatmap shows the most significant ones; `dmp` is already sorted by p-value.
max_heatmap_rows <- 2000
sig_cpgs <- head(dmp$CpG[dmp$significant], max_heatmap_rows)
if (length(sig_cpgs) >= 2) {
  ord <- order(pheno$Sample_Group)
  ht <- Heatmap(
    d$beta[sig_cpgs, ord],
    name = "Beta",
    col = circlize::colorRamp2(c(0, 0.5, 1), c("#2E86AB", "white", "#E4572E")),
    top_annotation = HeatmapAnnotation(Group = pheno$Sample_Group[ord], Batch = pheno$Batch[ord],
                                       col = list(Group = group_colours,
                                                  Batch = c("Chip shared by healthy + cancer" = "grey30",
                                                            "Single-group chip" = "grey80"))),
    column_split = pheno$Sample_Group[ord],
    cluster_column_slices = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    show_column_names = FALSE,
    column_title_gp = grid::gpar(fontsize = 10)
  )
  png(results_path("figures", "05_heatmap.png"), width = 9, height = 7, units = "in", res = 150)
  draw(ht, column_title = paste0("Top ", length(sig_cpgs), " of ", sum(dmp$significant),
                                 " DMPs (Cancer vs Healthy)"))
  invisible(dev.off())
}

message("Figures written to ", file.path(cfg$paths$results_dir, "figures"))
