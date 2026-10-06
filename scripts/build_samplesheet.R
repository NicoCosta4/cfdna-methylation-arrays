# One-off: build config/samplesheet.csv from the GEO metadata of GSE122126.
# Selects the plasma cfDNA samples profiled on EPIC (GPL21145) from healthy
# donors, solid-tumor patients and sepsis patients, and annotates:
#   Unit_ID     biological unit: technical replicates (same person or pool at a
#               different DNA input, or repeated arrays) share it and are averaged
#   Shared_Chip TRUE for chips that hold both healthy and cancer samples, used by
#               the batch-robustness sensitivity analysis
# The result is committed; rerun only to regenerate it.

geo_url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122126&targ=gsm&form=text&view=brief"
lines   <- readLines(geo_url, warn = FALSE)

# Split the GEO "brief" text into one block per sample
starts  <- grep("^\\^SAMPLE", lines)
blocks  <- split(lines, cumsum(seq_along(lines) %in% starts))[-1]

field <- function(block, key) sub(paste0("^!", key, " = "), "", grep(paste0("^!", key, " = "), block, value = TRUE))
characteristic <- function(block, name) {
  ch  <- field(block, "Sample_characteristics_ch1")
  hit <- grep(paste0("^", name, ":"), ch, value = TRUE)
  if (length(hit) == 0) NA_character_ else trimws(sub("^[^:]+:", "", hit))
}

meta <- do.call(rbind, lapply(blocks, function(b) {
  grn <- grep("Grn\\.idat", field(b, "Sample_supplementary_file"), value = TRUE)
  data.frame(
    GSM         = sub("^\\^SAMPLE = ", "", b[1]),
    title       = field(b, "Sample_title"),
    platform    = field(b, "Sample_platform_id"),
    sample_type = characteristic(b, "sample type"),
    Disease     = characteristic(b, "disease state"),
    Sex         = characteristic(b, "gender"),
    Age         = suppressWarnings(as.numeric(characteristic(b, "age"))),
    Basename    = if (length(grn)) sub("_Grn\\.idat\\.gz$", "", basename(grn)) else NA_character_
  )
}))

cancers <- c("Colon adenocarcinoma", "Breast adenocarcinoma", "Non-small-cell lung carcinoma",
             "Small-cell lung carcinoma", "Cancer of Unkown Primary")
ss <- meta[meta$platform == "GPL21145" & meta$sample_type == "cfDNA" &
             meta$Disease %in% c("Normal", "Sepsis", cancers), ]

ss$Sample_Name  <- sub("^cfDNA_", "", ss$title)
ss$Unit_ID      <- gsub("_(rep[0-9]+|[0-9]+ng)", "", ss$Sample_Name)
ss$Sample_Group <- ifelse(ss$Disease == "Normal", "Healthy",
                          ifelse(ss$Disease == "Sepsis", "Sepsis", "Cancer"))
ss$Sentrix_ID       <- sub("^GSM[0-9]+_([0-9]+)_.*$", "\\1", ss$Basename)
ss$Sentrix_Position <- sub("^.*_(R[0-9]{2}C[0-9]{2})$", "\\1", ss$Basename)
ss$Sex[ss$Sex == "NA"] <- NA

groups_per_chip <- tapply(ss$Sample_Group, ss$Sentrix_ID, function(g) all(c("Healthy", "Cancer") %in% g))
ss$Shared_Chip  <- groups_per_chip[ss$Sentrix_ID]

ss <- ss[order(factor(ss$Sample_Group, levels = c("Healthy", "Cancer", "Sepsis")), ss$Disease, ss$Sample_Name),
         c("Sample_Name", "GSM", "Unit_ID", "Sample_Group", "Disease", "Sex", "Age",
           "Sentrix_ID", "Sentrix_Position", "Shared_Chip", "Basename")]

write.csv(ss, "config/samplesheet.csv", row.names = FALSE)
message("Wrote ", nrow(ss), " samples (", length(unique(ss$Unit_ID)), " biological units) to config/samplesheet.csv")
print(table(Group = ss$Sample_Group[!duplicated(ss$Unit_ID)]))
