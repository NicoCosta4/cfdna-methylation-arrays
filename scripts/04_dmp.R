# Step 4 — Differentially methylated positions (DMPs) with limma.
# Statistics on M-values, effect sizes as delta-beta. One model with the three
# groups (Healthy, Cancer, Sepsis) gives:
#   - the primary contrast (Cancer - Healthy)
#   - Sepsis - Healthy, used later to discard inflammation-driven CpGs
# A sensitivity analysis refits Cancer vs Healthy using only the chips that hold
# both groups, with chip as a covariate, to check robustness to batch effects.

source("scripts/_common.R")
suppressPackageStartupMessages(library(limma))

d     <- readRDS(results_path("rds", "02_normalized.rds"))
pheno <- d$pheno
group <- factor(pheno[[cfg$dmp$group_column]], levels = c("Healthy", "Cancer", "Sepsis"))

# Promoter / gene body / intergenic, from the UCSC RefGene annotation
gene_region <- function(refgene_group) {
  parts <- strsplit(ifelse(is.na(refgene_group), "", refgene_group), ";")
  promoter_terms <- c("TSS1500", "TSS200", "5'UTR", "5UTR", "1stExon", "exon_1")
  vapply(parts, function(p) {
    if (length(p) == 0) return("intergenic")
    if (any(p %in% promoter_terms)) "promoter" else "body"
  }, character(1))
}
# RefGene lists one gene name per transcript ("RBP1;RBP1"); keep unique names
unique_genes <- function(x) {
  vapply(strsplit(ifelse(is.na(x), "", x), ";"),
         function(g) paste(unique(g), collapse = ";"), character(1))
}

# Fit `design` on the given units and return an annotated table for `contrast`
run_contrast <- function(units, design, contrast) {
  cm  <- makeContrasts(contrasts = contrast, levels = design)
  fit <- eBayes(contrasts.fit(lmFit(d$M[, units], design), cm))
  res <- topTable(fit, number = Inf, sort.by = "none")

  g <- factor(pheno[units, cfg$dmp$group_column])
  means <- sapply(levels(g), function(l) rowMeans(d$beta[, units][, g == l, drop = FALSE]))
  delta <- as.vector(means[, colnames(means) %in% rownames(cm)] %*%
                       cm[colnames(means)[colnames(means) %in% rownames(cm)], 1])

  ann <- d$annotation[rownames(res), ]
  out <- data.frame(
    CpG             = rownames(res),
    chr             = ann$chr,
    pos             = ann$pos,
    gene            = unique_genes(ann$UCSC_RefGene_Name),
    gene_region     = gene_region(ann$UCSC_RefGene_Group),
    island_relation = ann$Relation_to_Island,
    round(means, 3),
    delta_beta      = round(delta, 3),
    logFC_M         = round(res$logFC, 3),
    t               = round(res$t, 3),
    p_value         = signif(res$P.Value, 3),
    fdr             = signif(res$adj.P.Val, 3)
  )
  out$significant <- out$fdr < cfg$dmp$fdr_max & abs(out$delta_beta) >= cfg$dmp$min_abs_delta_beta
  out$direction   <- ifelse(out$delta_beta > 0, "hyper", "hypo")
  out[order(out$p_value), ]
}

# ---- Main model: all biological units, three groups ---------------------------
all_units <- rownames(pheno)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

dmp    <- run_contrast(all_units, design, cfg$dmp$contrast)
sepsis <- run_contrast(all_units, design, cfg$dmp$sepsis_contrast)

# ---- Sensitivity: Healthy and Cancer units on shared chips, chip as covariate ---
shared <- rownames(pheno)[pheno$Shared_Chip & group %in% c("Healthy", "Cancer")]
g_s    <- droplevels(group[match(shared, rownames(pheno))])
chip   <- factor(pheno[shared, "Sentrix_ID"])
design_s <- model.matrix(~0 + g_s + chip)
colnames(design_s) <- c(levels(g_s), paste0("chip", levels(chip)[-1]))
sens <- run_contrast(shared, design_s, cfg$dmp$contrast)

sig <- dmp[dmp$significant, ]
sens_sig <- sens[match(sig$CpG, sens$CpG), ]
sensitivity_summary <- data.frame(
  primary_dmps                = nrow(sig),
  same_direction_in_shared    = sum(sign(sens_sig$delta_beta) == sign(sig$delta_beta)),
  nominal_p05_same_direction  = sum(sens_sig$p_value < 0.05 & sign(sens_sig$delta_beta) == sign(sig$delta_beta)),
  delta_beta_correlation      = round(cor(sig$delta_beta, sens_sig$delta_beta), 3),
  units_in_sensitivity        = length(shared),
  chips_in_sensitivity        = nlevels(chip)
)

# ---- Outputs -----------------------------------------------------------------------
summary_table <- as.data.frame.matrix(table(factor(sig$gene_region, levels = c("promoter", "body", "intergenic")),
                                            factor(sig$direction, levels = c("hyper", "hypo"))))
write.csv(summary_table, results_path("tables", "04_dmp_summary.csv"))
write.csv(sig, results_path("tables", "04_dmp_significant.csv"), row.names = FALSE)
write.csv(sensitivity_summary, results_path("tables", "04_sensitivity_shared_chips.csv"), row.names = FALSE)
saveRDS(list(dmp = dmp, sepsis = sepsis, sensitivity = sens), results_path("rds", "04_dmp_all.rds"))

message(nrow(sig), " DMPs (FDR < ", cfg$dmp$fdr_max, ", |delta-beta| >= ", cfg$dmp$min_abs_delta_beta,
        "): ", sum(sig$direction == "hyper"), " hyper, ", sum(sig$direction == "hypo"),
        " hypo in ", cfg$dmp$contrast)
print(sensitivity_summary, row.names = FALSE)
