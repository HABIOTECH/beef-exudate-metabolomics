# ==============================================================================
# STAGE 1 · STEP 4 — Best-hit selection and isomer handling
# ------------------------------------------------------------------------------
# Within the 10 ppm window a feature can match several DB entries:
#   * Ranked candidates (different ppm) -> keep the lowest-ppm hit as best.
#   * True isobaric isomers (identical exact mass) -> MS1 cannot resolve them,
#     so their names are concatenated as "name A ; name B" (nothing discarded).
# Retention time is intentionally NOT used as a global match filter (see README,
# Design notes): RT drifts between batches and would silently drop true positives.
#
# INPUT  : 03_identified.xlsx
# OUTPUT : 04_best_hits.xlsx   (one curated row per feature; PPM_Error retained)
# ==============================================================================

library(readxl); library(writexl); library(dplyr)

# ---- CONFIG ------------------------------------------------------------------
INPUT_FILE  <- "results/03_identified.xlsx"
OUTPUT_FILE <- "results/04_best_hits.xlsx"
# ------------------------------------------------------------------------------

dat <- read_excel(INPUT_FILE)

best <- dat %>%
  group_by(Alignment_ID) %>%
  arrange(IDEOM_PPM, .by_group = TRUE) %>%        # lowest ppm error wins
  summarise(
    Exp_Mz        = first(Exp_Mz),
    Exp_Neutral   = first(Exp_Neutral),
    Best_Name     = first(IDEOM_Name),
    PPM_Error     = first(IDEOM_PPM),
    # merge any isobaric isomer names into a single auditable cell
    Isomer_Names  = paste(unique(na.omit(IDEOM_Name)), collapse = " ; "),
    HMDB_IDs      = first(HMDB_IDs),
    HMDB_Names    = first(HMDB_Names),
    .groups = "drop"
  )

write_xlsx(best, OUTPUT_FILE)
cat("Best hits per feature:", nrow(best), "\n")
