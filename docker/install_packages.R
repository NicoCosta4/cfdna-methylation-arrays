# Installs every R package the pipeline needs, then fails loudly if any is missing.
# BiocManager::install() only warns when a package fails, which would let a broken
# image build "successfully" — the check at the end prevents that.

pkgs <- c(
  # core analysis
  "minfi", "limma",
  # one-off download of the cross-reactive probe list (scripts/fetch_reference_lists.R)
  "ExperimentHub",
  # EPIC v1 (hg19) manifest + annotation
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  # tissue-of-origin deconvolution (non-negative least squares)
  "nnls",
  # biomarker filters: purified leukocyte methylomes and TCGA breast cancer data
  "FlowSorted.Blood.EPIC", "TCGAbiolinks",
  # fast reading of TCGA files, config and figures
  "yaml", "data.table", "ggplot2", "ComplexHeatmap", "circlize"
)

BiocManager::install(pkgs, ask = FALSE, update = FALSE, Ncpus = parallel::detectCores())

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Failed to install: ", paste(missing, collapse = ", "))
}
message("All ", length(pkgs), " packages installed.")
