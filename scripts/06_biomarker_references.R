# Step 6 — Reference methylomes for the biomarker filters (heavy downloads, run once;
# summaries are cached in data/reference/):
#   - Purified leukocytes (FlowSorted.Blood.EPIC, GSE110554): mean beta per cell type,
#     to measure each CpG's background methylation in blood.
#   - TCGA primary tumors vs normal tissue (450K) for the tumor types in the cohort
#     (colon, breast, lung), to check that a CpG is altered in the tumors themselves.
# The TCGA samples drawn at random are written to reference/tcga_samples.csv the first
# time; later runs reuse that list, so the selection does not change if GDC reorders.

source("scripts/_common.R")
suppressPackageStartupMessages({
  library(minfi)
  library(ExperimentHub)
  library(TCGAbiolinks)
  library(data.table)
})

ref_dir <- dirname(cfg$paths$atlas)
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
setExperimentHubOption("CACHE", file.path("data", "cache", "ExperimentHub"))
setExperimentHubOption("ASK", FALSE)

# ---- Purified leukocytes --------------------------------------------------------
blood_out <- file.path(ref_dir, "flowsorted_blood_means.rds")
if (!file.exists(blood_out)) {
  rg <- ExperimentHub()[["EH1136"]]
  rg <- rg[, rg$CellType != "MIX"]          # keep purified populations, drop mixtures
  beta <- getBeta(preprocessNoob(rg))
  blood_means <- sapply(split(seq_len(ncol(beta)), rg$CellType),
                        function(i) rowMeans(beta[, i, drop = FALSE], na.rm = TRUE))
  saveRDS(blood_means, blood_out)
  message("Leukocyte references: ", paste(colnames(blood_means), table(rg$CellType)[colnames(blood_means)],
                                         sep = "=", collapse = ", "))
}

# ---- TCGA tumors vs normal tissue ----------------------------------------------------
tcga_out     <- file.path(ref_dir, "tcga_tumor_normal.rds")
samples_file <- "reference/tcga_samples.csv"
projects     <- cfg$biomarker$tcga_projects
n_per_group  <- cfg$biomarker$tcga_samples_per_group

if (!file.exists(tcga_out)) {
  query_project <- function(project) {
    GDCquery(project = project, data.category = "DNA Methylation",
             data.type = "Methylation Beta Value", platform = "Illumina Human Methylation 450")
  }

  if (!file.exists(samples_file)) {
    set.seed(1)
    chosen <- do.call(rbind, lapply(projects, function(p) {
      r <- getResults(query_project(p))
      pick <- function(type) {
        idx <- which(r$sample_type == type)
        idx[sample.int(length(idx), min(n_per_group, length(idx)))]
      }
      sel <- r[c(pick("Primary Tumor"), pick("Solid Tissue Normal")), ]
      data.frame(project = p, barcode = sel$cases, sample_type = sel$sample_type,
                 file_id = sel$file_id, file_name = sel$file_name)
    }))
    write.csv(chosen, samples_file, row.names = FALSE)
  }
  chosen <- read.csv(samples_file, stringsAsFactors = FALSE)

  summaries <- lapply(projects, function(p) {
    q <- query_project(p)
    q$results[[1]] <- q$results[[1]][q$results[[1]]$file_id %in% chosen$file_id[chosen$project == p], ]
    GDCdownload(q, directory = file.path("data", "tcga"), files.per.chunk = 10)

    files <- list.files(file.path("data", "tcga", p), pattern = "level3betas\\.txt$",
                        recursive = TRUE, full.names = TRUE)
    files <- files[basename(files) %in% chosen$file_name[chosen$project == p]]
    type  <- chosen$sample_type[match(basename(files), chosen$file_name)]

    beta <- do.call(cbind, lapply(files, function(f) {
      x <- fread(f, header = FALSE, col.names = c("probe", "beta"))
      setNames(x$beta, x$probe)
    }))
    tumor  <- rowMeans(beta[, type == "Primary Tumor", drop = FALSE], na.rm = TRUE)
    normal <- rowMeans(beta[, type == "Solid Tissue Normal", drop = FALSE], na.rm = TRUE)
    message(p, ": ", sum(type == "Primary Tumor"), " tumors, ", sum(type == "Solid Tissue Normal"), " normals")
    data.frame(CpG = rownames(beta), project = p, tumor = tumor, normal = normal, delta_beta = tumor - normal)
  })
  saveRDS(do.call(rbind, summaries), tcga_out)
}

message("References ready in ", ref_dir)
