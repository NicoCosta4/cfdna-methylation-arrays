# Step 2 — Normalization (stratified quantile normalization, minfi::preprocessQuantile),
# probe filtering, and averaging of technical replicates into biological units.
# Every filter is logged so the README can report how many probes each removed.

source("scripts/_common.R")
suppressPackageStartupMessages({
  library(minfi)
  library(limma)
  library(ggplot2)
})

qc     <- readRDS(results_path("rds", "01_qc.rds"))
rgSet  <- qc$rgSet
detP   <- qc$detP
groups <- rgSet$Sample_Group

mSetRaw <- preprocessRaw(rgSet)
# preprocessQuantile estimates each sample's sex to normalize chrX/chrY separately.
# With cfDNA and pooled samples this estimate can be "inconsistent" (a warning);
# harmless here, as sex-chromosome probes are removed below.
gset <- preprocessQuantile(rgSet)

# ---- Beta-value densities before and after normalization ---------------------
# A random subset of probes keeps the plot fast; the shape is identical.
set.seed(1)
probes <- sample(intersect(featureNames(gset), featureNames(mSetRaw)), 50000)

density_data <- function(betas, label) {
  data.frame(
    beta   = as.vector(betas[probes, ]),
    Sample = rep(colnames(betas), each = length(probes)),
    Group  = rep(groups, each = length(probes)),
    Data   = label
  )
}
# offset = 100 (Illumina's default) avoids 0/0 betas for probes with no signal
dens <- rbind(density_data(getBeta(mSetRaw, offset = 100), "Raw"),
              density_data(getBeta(gset), "Normalized"))
dens$Data <- factor(dens$Data, levels = c("Raw", "Normalized"))

p_dens <- ggplot(dens, aes(beta, group = Sample, colour = Group)) +
  geom_density(linewidth = 0.3, alpha = 0.6) +
  facet_wrap(~Data) +
  labs(x = "Beta value", y = "Density", colour = NULL) +
  theme_minimal()
ggsave(results_path("figures", "02_beta_densities.png"), p_dens, width = 9, height = 4, dpi = 150)

# Unfiltered normalized betas are kept for deconvolution: reference atlases use
# their own CpG sets, which our filters would otherwise shrink unnecessarily.
beta_unfiltered <- getBeta(gset)

# ---- Probe filtering ----------------------------------------------------------
filter_log <- data.frame(filter = "Probes after normalization", removed = NA_integer_,
                         remaining = nrow(gset))

apply_filter <- function(gset, keep, name) {
  filter_log <<- rbind(filter_log, data.frame(filter = name, removed = sum(!keep),
                                              remaining = sum(keep)))
  gset[keep, ]
}

detP <- detP[match(featureNames(gset), rownames(detP)), ]
gset <- apply_filter(gset, rowSums(detP < cfg$qc$probe_detp_max) == ncol(gset),
                     sprintf("Detection p >= %s in any sample", cfg$qc$probe_detp_max))

gset <- apply_filter(gset, startsWith(featureNames(gset), "cg"), "Non-CpG probes (ch.*)")

if (cfg$filtering$remove_sex_chromosomes) {
  gset <- apply_filter(gset, !as.character(seqnames(gset)) %in% c("chrX", "chrY"),
                       "Sex chromosomes (chrX, chrY)")
}

if (cfg$filtering$remove_snp_probes) {
  no_snps <- featureNames(dropLociWithSnps(gset))
  gset <- apply_filter(gset, featureNames(gset) %in% no_snps,
                       "SNP at CpG or single-base extension site")
}

if (cfg$filtering$remove_cross_reactive) {
  crosshyb <- readLines("reference/crosshyb_probes.txt")
  gset <- apply_filter(gset, !featureNames(gset) %in% crosshyb, "Cross-reactive probes")
}

write.csv(filter_log, results_path("tables", "02_probe_filtering.csv"), row.names = FALSE)
print(filter_log, row.names = FALSE)

# ---- Average technical replicates into biological units ------------------------
# Arrays from the same person or pool (different DNA input, repeated arrays) are
# not independent; averaging them makes every unit count once in the statistics.
pheno <- as.data.frame(colData(gset))
unit  <- pheno[[cfg$dmp$unit_column]]

M_unit    <- avereps(t(getM(gset)), ID = unit)
beta_unit <- avereps(t(getBeta(gset)), ID = unit)
beta_unfiltered_unit <- avereps(t(beta_unfiltered[, colnames(gset)]), ID = unit)

pheno_unit <- pheno[!duplicated(unit), c(cfg$dmp$unit_column, "Sample_Group", "Disease", "Sex",
                                         "Sentrix_ID", "Shared_Chip")]
pheno_unit$n_arrays <- as.vector(table(unit)[pheno_unit[[cfg$dmp$unit_column]]])
pheno_unit$Age      <- tapply(pheno$Age, unit, median, na.rm = TRUE)[pheno_unit[[cfg$dmp$unit_column]]]
rownames(pheno_unit) <- pheno_unit[[cfg$dmp$unit_column]]

ann <- as.data.frame(getAnnotation(gset))
ann <- ann[, c("chr", "pos", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Group")]

saveRDS(list(beta = t(beta_unit)[, rownames(pheno_unit)],
             M = t(M_unit)[, rownames(pheno_unit)],
             beta_unfiltered = t(beta_unfiltered_unit)[, rownames(pheno_unit)],
             pheno = pheno_unit, annotation = ann),
        results_path("rds", "02_normalized.rds"))
message("Normalization done: ", nrow(gset), " CpGs x ", ncol(gset), " arrays -> ",
        nrow(pheno_unit), " biological units")
