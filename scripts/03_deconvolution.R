# Step 3 — Tissue-of-origin deconvolution with the Moss et al. (2018) atlas.
# Each cfDNA methylome is modelled as a mixture of 25 reference cell-type
# methylomes; non-negative least squares (NNLS) estimates each contribution,
# which is then rescaled to sum to 1 (same approach as the original method).

source("scripts/_common.R")
suppressPackageStartupMessages({
  library(nnls)
  library(ggplot2)
})

d <- readRDS(results_path("rds", "02_normalized.rds"))

# The atlas has 7,890 rows but 6,105 distinct CpGs: a CpG that marks several cell
# types is listed once per type (identical copies). The original implementation
# (deconvolve.py) keeps those copies, so we do too, to reproduce the method.
atlas_df <- read.csv(cfg$paths$atlas, check.names = FALSE)
atlas_df <- atlas_df[atlas_df[[1]] %in% rownames(d$beta_unfiltered) & stats::complete.cases(atlas_df), ]
cpgs  <- atlas_df[[1]]
atlas <- as.matrix(atlas_df[, -1])
colnames(atlas) <- sub("_EPIC$", "", colnames(atlas))
message(length(unique(cpgs)), " distinct atlas CpGs (", length(cpgs), " rows) are present in the data")

deconvolve <- function(sample_beta) {
  ok <- !is.na(sample_beta)
  x  <- nnls(atlas[ok, ], sample_beta[ok])$x
  x / sum(x)
}
fractions <- sapply(colnames(d$beta_unfiltered), function(s) deconvolve(d$beta_unfiltered[cpgs, s]))
rownames(fractions) <- colnames(atlas)

blood_types <- c("Monocytes", "B-cells", "CD4T-cells", "NK-cells", "CD8T-cells",
                 "Neutrophils", "Erythrocyte_progenitors")
group_levels <- c("Healthy", "Cancer", "Sepsis")

frac_df <- data.frame(
  Unit_ID      = colnames(fractions),
  Sample_Group = factor(d$pheno$Sample_Group, levels = group_levels),
  Disease      = d$pheno$Disease,
  t(fractions),
  blood_total  = colSums(fractions[blood_types, ]),
  check.names  = FALSE
)
write.csv(frac_df, results_path("tables", "03_cell_type_fractions.csv"), row.names = FALSE)
saveRDS(frac_df, results_path("rds", "03_deconvolution.rds"))

# ---- Stacked bar chart: cell types contributing >= 1% on average, rest as "Other"
long <- data.frame(
  Unit_ID      = rep(colnames(fractions), each = nrow(fractions)),
  Sample_Group = rep(frac_df$Sample_Group, each = nrow(fractions)),
  cell_type    = rep(rownames(fractions), times = ncol(fractions)),
  fraction     = as.vector(fractions)
)
major <- names(which(rowMeans(fractions) >= 0.01))
long$cell_type <- ifelse(long$cell_type %in% major, long$cell_type, "Other")
long <- aggregate(fraction ~ Unit_ID + Sample_Group + cell_type, long, sum)
long$cell_type <- factor(long$cell_type, levels = c(major[order(-rowMeans(fractions)[major])], "Other"))

fill_colours <- c(setNames(scales::hue_pal()(length(major)), levels(long$cell_type)[seq_along(major)]),
                  Other = "grey80")
p_bar <- ggplot(long, aes(Unit_ID, fraction, fill = cell_type)) +
  geom_col() +
  facet_grid(~Sample_Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = fill_colours) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Estimated contribution to cfDNA", fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
ggsave(results_path("figures", "03_tissue_of_origin.png"), p_bar, width = 12, height = 5.5, dpi = 150)

# ---- Group comparison for the cell types of interest (descriptive)
of_interest <- c("blood_total", "Neutrophils", "Colon_epithelial_cells", "Breast", "Lung_cells",
                 "Hepatocytes", "Vascular_endothelial_cells")
cmp <- do.call(rbind, lapply(of_interest, function(ct) {
  data.frame(cell_type = ct, Sample_Group = frac_df$Sample_Group, fraction = frac_df[[ct]])
}))
cmp$cell_type <- factor(cmp$cell_type, levels = of_interest)

p_cmp <- ggplot(cmp, aes(Sample_Group, fraction, colour = Sample_Group)) +
  geom_boxplot(outlier.shape = NA, colour = "grey60") +
  geom_jitter(width = 0.15, height = 0, size = 1.8) +  # height = 0: never move values vertically
  facet_wrap(~cell_type, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = c(Healthy = "#2E86AB", Cancer = "#E4572E", Sepsis = "#7D5BA6")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Estimated contribution", colour = NULL) +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(results_path("figures", "03_fractions_by_group.png"), p_cmp, width = 11, height = 5.5, dpi = 150)

# Each disease group against healthy donors (Wilcoxon rank-sum, two-sided)
wilcox <- do.call(rbind, lapply(of_interest, function(ct) {
  healthy <- frac_df[[ct]][frac_df$Sample_Group == "Healthy"]
  do.call(rbind, lapply(c("Cancer", "Sepsis"), function(g) {
    other <- frac_df[[ct]][frac_df$Sample_Group == g]
    data.frame(cell_type = ct, comparison = paste(g, "vs Healthy"),
               median_group = median(other), median_healthy = median(healthy),
               p_wilcoxon = wilcox.test(other, healthy, exact = FALSE)$p.value)
  }))
}))
write.csv(wilcox, results_path("tables", "03_fractions_by_group.csv"), row.names = FALSE)
print(wilcox, row.names = FALSE, digits = 3)
