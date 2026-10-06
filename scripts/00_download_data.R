# Download the IDAT files of the selected GSE122126 samples (see
# config/samplesheet.csv) and the Moss et al. (2018) reference atlas.
# Only the 49 cfDNA arrays used here are fetched, not the whole series.

source("scripts/_common.R")

ss        <- read.csv(cfg$paths$samplesheet, stringsAsFactors = FALSE)
idat_dir  <- cfg$paths$idat_dir
max_tries <- 3
dir.create(idat_dir, recursive = TRUE, showWarnings = FALSE)

# A download interrupted by a dropped connection or Ctrl+C leaves a truncated
# file behind; `gzip -t` fails on it, so we treat it as missing.
gz_is_complete <- function(path) {
  file.exists(path) && system2("gzip", c("-t", shQuote(path)), stdout = FALSE, stderr = FALSE) == 0
}

download_with_retries <- function(url, dest) {
  for (attempt in seq_len(max_tries)) {
    if (gz_is_complete(dest)) return(invisible(TRUE))
    unlink(dest)
    tryCatch(download.file(url, dest, mode = "wb", quiet = TRUE),
             error = function(e) message("  attempt ", attempt, " failed: ", conditionMessage(e)),
             warning = function(w) message("  attempt ", attempt, " failed: ", conditionMessage(w)))
  }
  if (!gz_is_complete(dest)) {
    unlink(dest)
    stop("Could not download ", basename(dest), " after ", max_tries, " attempts. ",
         "The NCBI server refused or dropped the connection; retry later.")
  }
}

# GEO stores each sample's files in .../samples/GSMnnnnnnn/GSMxxxxxxx/suppl/
gsm_dir <- function(gsm) {
  paste0("https://ftp.ncbi.nlm.nih.gov/geo/samples/", sub("[0-9]{3}$", "nnn", gsm), "/", gsm, "/suppl/")
}

for (i in seq_len(nrow(ss))) {
  for (channel in c("Grn", "Red")) {
    file <- paste0(ss$Basename[i], "_", channel, ".idat.gz")
    dest <- file.path(idat_dir, file)
    if (!gz_is_complete(dest)) {
      message("[", i, "/", nrow(ss), "] ", file)
      download_with_retries(paste0(gsm_dir(ss$GSM[i]), file), dest)
    }
  }
}

n_idats  <- length(list.files(idat_dir, pattern = "\\.idat\\.gz$"))
expected <- 2 * nrow(ss)
if (n_idats < expected) stop("Expected ", expected, " IDAT files, found ", n_idats)
message("Done: ", n_idats, " IDAT files in ", idat_dir)

# ---- Tissue-of-origin reference atlas (Moss et al. 2018) ----------------------
# Distributed under a research-only license that does not allow redistribution,
# so it is downloaded here (pinned to a fixed commit + checksum) instead of being
# versioned in this repository.
atlas_url <- paste0("https://raw.githubusercontent.com/nloyfer/meth_atlas/",
                    "fd4bc448d14b7f66300fa9496acd1789f1267817/reference_atlas.csv")
atlas_md5 <- "c64af0b9d2c5f8648e61897067b534f7"
atlas_path <- cfg$paths$atlas

if (!file.exists(atlas_path) || tools::md5sum(atlas_path) != atlas_md5) {
  dir.create(dirname(atlas_path), recursive = TRUE, showWarnings = FALSE)
  download.file(atlas_url, atlas_path, quiet = TRUE)
}
if (tools::md5sum(atlas_path) != atlas_md5) {
  unlink(atlas_path)
  stop("Checksum mismatch for the reference atlas downloaded from ", atlas_url)
}
message("Done: ", atlas_path)
