# Shared helpers sourced by every pipeline step. Run scripts from the repo root.

suppressPackageStartupMessages(library(yaml))

cfg <- yaml::read_yaml("config/config.yml")

# minfi annotation slot for each supported array. EPIC v1 IDATs are detected
# automatically, but EPICv2 must be set by hand, so we always set it explicitly.
array_annotation <- function(array) {
  switch(array,
    EPIC   = c(array = "IlluminaHumanMethylationEPIC",   annotation = "ilm10b4.hg19"),
    EPICv2 = c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38"),
    stop("Unknown array type in config.yml: ", array)
  )
}

# Build a path inside the results directory, creating parent folders as needed.
results_path <- function(...) {
  path <- file.path(cfg$paths$results_dir, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}
