# ==============================================================================
# STAGE 1 · STEP 5 — Batch annotation across cohorts
# ------------------------------------------------------------------------------
# Applies the identification logic to every experimental cohort in one loop and
# writes one dual-sheet workbook per cohort:
#   Sheet 1 "Clean Data"            : original feature table + Metabolite Name
#   Sheet 2 "Identification Details": identified features joined to IDEOM metadata
#
# A cohort = one (treatment x matrix) feature table. Add your own cohort files to
# INPUT_FILES below — labels are generic on purpose.
#
# INPUTS : 03_identified.xlsx, ideom_reference.xlsx, per-cohort feature tables
# OUTPUT : results/annotated/Annotated_<cohort>.xlsx
# ==============================================================================

library(readxl); library(writexl); library(dplyr)

# ---- CONFIG ------------------------------------------------------------------
RESULTS_FILE  <- "results/03_identified.xlsx"
IDEOM_FILE    <- "data/ideom_reference.xlsx"
OUTPUT_FOLDER <- "results/annotated"

# One entry per cohort (treatment x matrix). Edit to match your study.
INPUT_FILES <- c(
  "data/cohort_matrix1_pathogen1.xlsx",
  "data/cohort_matrix2_pathogen1.xlsx",
  "data/cohort_matrix1_pathogen2.xlsx",
  "data/cohort_matrix2_pathogen2.xlsx"
)

METADATA_COLS <- c("Alignment_ID", "Metabolite", "Formula", "Map", "DB",
  "Pathway", "KEGGpath", "KEGGid#", "Lmid#", "HMDBID#", "MetacycID#",
  "CAS#", "CHEBI#", "METLIN#", "InchiKEY", "SMILES", "Synonyms",
  "source_id", "HMDB_taxonomy", "CATEGORY", "species",
  "alt_HMDBID#", "alt_CAS#", "alt_CHEBI#", "alt_METLIN#")
# ------------------------------------------------------------------------------

if (!dir.exists(OUTPUT_FOLDER)) dir.create(OUTPUT_FOLDER, recursive = TRUE)

results_file <- read_excel(RESULTS_FILE)
ideom_db     <- read_excel(IDEOM_FILE)

for (file_name in INPUT_FILES) {
  if (!file.exists(file_name)) { message("Skipping missing: ", file_name); next }
  data_file <- read_excel(file_name)

  # Sheet 1 — Clean Data with Metabolite Name merged in
  sheet1 <- data_file %>%
    left_join(results_file %>% select(Alignment_ID, Metabolite = Best_Name),
              by = c("Alignment ID" = "Alignment_ID")) %>%
    rename(`Metabolite Name` = Metabolite)

  # Sheet 2 — full IDEOM metadata for identified features only
  identified_ids <- sheet1 %>%
    filter(!is.na(`Metabolite Name`), `Metabolite Name` != "Not Identified") %>%
    pull(`Alignment ID`)

  sheet2 <- results_file %>%
    filter(Alignment_ID %in% identified_ids) %>%
    left_join(ideom_db, by = c("Best_Name" = "Metabolite")) %>%
    select(any_of(METADATA_COLS))

  out_name <- file.path(OUTPUT_FOLDER, paste0("Annotated_", basename(file_name)))
  write_xlsx(list("Clean Data" = sheet1, "Identification Details" = sheet2), out_name)
  message("Written: ", out_name)
}
