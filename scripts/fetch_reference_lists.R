# One-off: download the list of cross-reactive probes and store it as plain text
# in reference/, so the pipeline runs offline and the exact list is versioned.
# The list is the one DMRcate uses (ExperimentHub EH3129): probes whose sequence
# aligns to non-target regions of the genome (>= 47 bp match), compiled from
# Chen et al. 2013 (450K) and Pidsley et al. 2016 (EPIC).
#
# You only need to rerun this to refresh the list; the result is committed.

suppressPackageStartupMessages(library(ExperimentHub))

out_file <- "reference/crosshyb_probes.txt"

crosshyb <- ExperimentHub()[["EH3129"]]
crosshyb <- sort(unique(as.character(crosshyb)))

dir.create(dirname(out_file), showWarnings = FALSE)
writeLines(crosshyb, out_file)
message("Wrote ", length(crosshyb), " probe IDs to ", out_file)
