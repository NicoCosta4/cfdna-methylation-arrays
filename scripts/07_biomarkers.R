# Step 7 — Biomarker funnel: from cfDNA DMPs to a tumor-derived signature with low
# background in blood. Each filter answers one question:
#   1. Does the CpG differ between cancer and healthy cfDNA?            (step 4)
#   2. Is the change NOT reproduced by systemic inflammation (sepsis)?  (step 4)
#   3. Is the opposite state found in every leukocyte type?             (FlowSorted)
#   4. Is the CpG altered the same way in the tumors themselves?        (TCGA)
# Thresholds live in config.yml (biomarker section).

source("scripts/_common.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
})

b      <- cfg$biomarker
d      <- readRDS(results_path("rds", "02_normalized.rds"))
dmps   <- readRDS(results_path("rds", "04_dmp_all.rds"))
blood  <- readRDS(file.path(dirname(cfg$paths$atlas), "flowsorted_blood_means.rds"))
tcga   <- readRDS(file.path(dirname(cfg$paths$atlas), "tcga_tumor_normal.rds"))

funnel <- data.frame()
log_step <- function(cand, name) {
  funnel <<- rbind(funnel, data.frame(step = name, hyper = sum(cand$direction == "hyper"),
                                      hypo = sum(cand$direction == "hypo"), total = nrow(cand)))
  cand
}

# 1 — cfDNA DMPs (Cancer vs Healthy)
cand <- log_step(dmps$dmp[dmps$dmp$significant, ], "1. DMP in cancer vs healthy cfDNA")

# 2 — drop CpGs that shift in sepsis in the same direction as in cancer
cand$sepsis_delta_beta <- dmps$sepsis$delta_beta[match(cand$CpG, dmps$sepsis$CpG)]
shifted_in_sepsis <- sign(cand$sepsis_delta_beta) == sign(cand$delta_beta) &
  abs(cand$sepsis_delta_beta) >= b$sepsis_max_same_direction_delta
cand <- log_step(cand[!shifted_in_sepsis, ], "2. Not shifted in sepsis")

# 3 — low background in blood: opposite methylation state in every leukocyte type
bm <- blood[match(cand$CpG, rownames(blood)), , drop = FALSE]
cand$blood_max_beta <- round(apply(bm, 1, max), 3)
cand$blood_min_beta <- round(apply(bm, 1, min), 3)
low_background <- ifelse(cand$direction == "hyper",
                         cand$blood_max_beta <= b$blood_max_beta_for_hyper,
                         cand$blood_min_beta >= b$blood_min_beta_for_hypo)
cand <- log_step(cand[low_background %in% TRUE, ], "3. Low background in leukocytes")

# 4 — altered in TCGA tumors vs normal tissue, same direction, in >= N tumor types.
# CpGs absent from the 450K array cannot be checked and therefore do not pass.
tcga_delta <- sapply(b$tcga_projects, function(p) {
  x <- tcga[tcga$project == p, ]
  x$delta_beta[match(cand$CpG, x$CpG)]
})
tcga_delta <- matrix(tcga_delta, nrow = nrow(cand), dimnames = list(cand$CpG, b$tcga_projects))
agrees <- sign(tcga_delta) == sign(cand$delta_beta) & abs(tcga_delta) >= b$tcga_min_abs_delta_beta
cand$tcga_types_agreeing <- rowSums(agrees, na.rm = TRUE)
cand <- cbind(cand, setNames(as.data.frame(round(tcga_delta, 3)), paste0("delta_", sub("TCGA-", "", b$tcga_projects))))
signature <- log_step(cand[cand$tcga_types_agreeing >= b$tcga_min_projects, ], "4. Altered in TCGA tumors")

write.csv(funnel, results_path("tables", "07_biomarker_funnel.csv"), row.names = FALSE)
write.csv(signature, results_path("tables", "07_signature.csv"), row.names = FALSE)
print(funnel, row.names = FALSE)

# ---- Funnel figure -------------------------------------------------------------------
funnel_long <- rbind(data.frame(step = funnel$step, direction = "hyper", n = funnel$hyper),
                     data.frame(step = funnel$step, direction = "hypo", n = funnel$hypo))
funnel_long$step <- factor(funnel_long$step, levels = rev(funnel$step))
# Side-by-side bars: stacking on a log scale would distort the proportions
p_funnel <- ggplot(funnel_long, aes(n, step, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = format(n, big.mark = ",")), position = position_dodge(width = 0.8),
            hjust = -0.15, size = 3) +
  scale_x_log10(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = c(hyper = "#E4572E", hypo = "#2E86AB")) +
  labs(x = "CpGs (log scale)", y = NULL, fill = NULL, title = "Biomarker selection funnel") +
  theme_minimal()
ggsave(results_path("figures", "07_funnel.png"), p_funnel, width = 8, height = 3.5, dpi = 150)

# ---- Core panel: compact, assay-sized subset of the signature --------------------------
# Hypermethylated only (more robust in cfDNA than loss of methylation), altered in every
# TCGA tumor type, near-zero blood background, in CpG islands of promoters (biologically
# interpretable, assay-friendly); ranked by cfDNA significance.
eligible <- signature[signature$direction == "hyper" &
                        signature$tcga_types_agreeing == length(b$tcga_projects) &
                        signature$blood_max_beta <= b$panel_blood_max_beta &
                        signature$island_relation %in% b$panel_island_relations &
                        signature$gene_region %in% b$panel_gene_regions, ]
eligible <- eligible[order(eligible$p_value), ]

# Regional support: several eligible CpGs in the same promoter island point to a
# regional methylation change rather than a single noisy probe.
eligible$cpgs_in_gene <- as.vector(table(eligible$gene)[eligible$gene])
if (isTRUE(b$panel_one_cpg_per_gene)) eligible <- eligible[!duplicated(eligible$gene), ]

panel <- head(eligible, b$panel_size)
write.csv(panel, results_path("tables", "07_core_panel.csv"), row.names = FALSE)
message("Core panel: ", nrow(panel), " CpGs in ", length(unique(panel$gene)), " genes; ",
        sum(panel$cpgs_in_gene > 1), " genes supported by >1 eligible CpG")

# ---- Heatmap and per-sample score for a set of CpGs -----------------------------------
pheno <- d$pheno
pheno$Sample_Group <- factor(pheno$Sample_Group, levels = c("Healthy", "Cancer", "Sepsis"))
ord <- order(pheno$Sample_Group, pheno$Disease)
group_colours <- c(Healthy = "#2E86AB", Cancer = "#E4572E", Sepsis = "#7D5BA6")

plot_cpg_set <- function(set, label, file_prefix) {
  if (nrow(set) < 2) return(invisible(NULL))
  row_labels <- ifelse(set$gene == "", set$CpG, paste0(set$CpG, " (", set$gene, ")"))
  if ("cpgs_in_gene" %in% names(set)) {
    row_labels <- paste0(row_labels, ifelse(set$cpgs_in_gene > 1, paste0(" [", set$cpgs_in_gene, " CpGs]"), ""))
  }
  ht <- Heatmap(
    d$beta[set$CpG, ord],
    name = "Beta",
    col = circlize::colorRamp2(c(0, 0.5, 1), c("#2E86AB", "white", "#E4572E")),
    top_annotation = HeatmapAnnotation(Group = pheno$Sample_Group[ord], col = list(Group = group_colours)),
    left_annotation = rowAnnotation(Direction = set$direction,
                                    col = list(Direction = c(hyper = "#E4572E", hypo = "#2E86AB"))),
    column_split = pheno$Sample_Group[ord],
    cluster_column_slices = FALSE,
    row_labels = row_labels,
    show_row_names = nrow(set) <= 60,
    row_names_gp = grid::gpar(fontsize = 7),
    column_labels = pheno$Unit_ID[ord],
    column_names_gp = grid::gpar(fontsize = 6)
  )
  png(results_path("figures", paste0(file_prefix, "_heatmap.png")), width = 10, height = 7, units = "in", res = 150)
  draw(ht, column_title = paste(nrow(set), "CpG", label))
  invisible(dev.off())

  # Tumor-like methylation: beta for hypermethylated markers, 1 - beta for hypomethylated
  tumor_like <- ifelse(set$direction == "hyper", 1, -1) * (d$beta[set$CpG, , drop = FALSE] - 0.5) + 0.5
  scores <- data.frame(pheno[, c("Unit_ID", "Sample_Group", "Disease")], score = colMeans(tumor_like))
  write.csv(scores, results_path("tables", paste0(file_prefix, "_scores.csv")), row.names = FALSE)

  p_score <- ggplot(scores, aes(Sample_Group, score, colour = Sample_Group)) +
    geom_boxplot(outlier.shape = NA, colour = "grey60") +
    geom_jitter(width = 0.15, height = 0, size = 2) +
    scale_colour_manual(values = group_colours) +
    labs(x = NULL, y = "Score (mean tumor-like methylation)",
         title = paste0(tools::toTitleCase(label), " score (", nrow(set), " CpGs)")) +
    theme_minimal() +
    theme(legend.position = "none")
  ggsave(results_path("figures", paste0(file_prefix, "_score.png")), p_score, width = 5.5, height = 4.5, dpi = 150)
}

plot_cpg_set(signature, "signature", "07_signature")
plot_cpg_set(panel, "core panel", "07_core_panel")
