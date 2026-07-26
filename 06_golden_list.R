# ==============================================================================
# STAGE 1 · STEP 6 — Golden-List extraction
# ------------------------------------------------------------------------------
# After downstream statistics produce a short list of significant Alignment IDs,
# this step pulls the full metadata for exactly those IDs so the results table,
# pathway map and supplementary files can be built from one curated workbook.
#
# INPUTS : 03_identified.xlsx (or 04_best_hits.xlsx) + target_ids.xlsx (col: Alignment)
# OUTPUTS: 06_golden_list.xlsx            (full metadata for the target set)
#          06_golden_list_ids.xlsx        (slim ID-only, for pathway tools)
# ==============================================================================

library(readxl); library(writexl); library(dplyr)

# ---- CONFIG ------------------------------------------------------------------
MAIN_FILE    <- "results/03_identified.xlsx"
TARGET_FILE  <- "data/target_ids.xlsx"          # one column: Alignment
OUTPUT_FULL  <- "results/06_golden_list.xlsx"
OUTPUT_IDS   <- "results/06_golden_list_ids.xlsx"
# ------------------------------------------------------------------------------

main_data  <- read_excel(MAIN_FILE)
target_ids <- read_excel(TARGET_FILE)$Alignment

# Full metadata for the Golden List
extracted_full <- main_data %>% filter(Alignment_ID %in% target_ids)
write_xlsx(extracted_full, OUTPUT_FULL)

# Slim ID-only version for pathway tools (e.g., MetaboAnalyst, Cytoscape)
extracted_ids <- main_data %>%
  filter(Alignment_ID %in% target_ids) %>%
  select(any_of(c("Alignment_ID", "KEGGid#", "HMDBID#")))
write_xlsx(extracted_ids, OUTPUT_IDS)

cat("Golden List size:", nrow(extracted_full), "\n")
