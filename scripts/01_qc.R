# Step 1 — Quality control: read IDATs, compute detection p-values, drop samples
# whose mean detection p-value is above the threshold in config.yml, and record
# how groups are distributed across chips (batch structure).

source("scripts/_common.R")
suppressPackageStartupMessages({
  library(minfi)
  library(ggplot2)
})

ss <- read.csv(cfg$paths$samplesheet, stringsAsFactors = FALSE)
ss$Basename <- file.path(cfg$paths$idat_dir, ss$Basename)

# The dataset spans two EPIC revisions (1,052,641 vs 1,051,815 probes, a
# 826-probe difference); force = TRUE combines them on their shared probes.
rgSet <- read.metharray.exp(targets = ss, force = TRUE)
rgSet@annotation <- array_annotation(cfg$array)
sampleNames(rgSet) <- ss$Sample_Name

detP <- detectionP(rgSet)

qc <- data.frame(
  Sample_Name       = colnames(detP),
  Sample_Group      = ss$Sample_Group,
  Sentrix_ID        = ss$Sentrix_ID,
  mean_detP         = colMeans(detP),
  failed_probes_pct = 100 * colMeans(detP > cfg$qc$probe_detp_max)
)
qc$pass <- qc$mean_detP < cfg$qc$sample_detp_max
write.csv(qc, results_path("tables", "01_qc_samples.csv"), row.names = FALSE)

chips <- as.data.frame.matrix(table(ss$Sentrix_ID, ss$Sample_Group))
write.csv(chips, results_path("tables", "01_chips_by_group.csv"))

qc$Sample_Name <- factor(qc$Sample_Name, levels = qc$Sample_Name)
p_failed <- ggplot(qc, aes(Sample_Name, failed_probes_pct, fill = Sample_Group)) +
  geom_col() +
  facet_grid(~Sample_Group, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = paste0("% probes with detection p > ", cfg$qc$probe_detp_max), fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), legend.position = "none")
ggsave(results_path("figures", "01_failed_probes.png"), p_failed, width = 11, height = 4.5, dpi = 150)

# minfi intensity QC: samples below the dashed line ((mMed + uMed) / 2 < 10.5)
# have low overall signal and are usually failed hybridizations.
intensity <- as.data.frame(getQC(preprocessRaw(rgSet)))
intensity$Sample_Group <- ss$Sample_Group

p_int <- ggplot(intensity, aes(mMed, uMed, colour = Sample_Group)) +
  geom_abline(intercept = 21, slope = -1, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(x = "Median methylated intensity (log2)", y = "Median unmethylated intensity (log2)", colour = NULL) +
  theme_minimal()
ggsave(results_path("figures", "01_median_intensities.png"), p_int, width = 6, height = 5, dpi = 150)

if (any(!qc$pass)) {
  message("Removing samples that fail QC: ", paste(qc$Sample_Name[!qc$pass], collapse = ", "))
}
rgSet <- rgSet[, qc$pass]
detP  <- detP[, qc$pass]

saveRDS(list(rgSet = rgSet, detP = detP), results_path("rds", "01_qc.rds"))
message("QC done: ", sum(qc$pass), "/", nrow(qc), " samples kept")
