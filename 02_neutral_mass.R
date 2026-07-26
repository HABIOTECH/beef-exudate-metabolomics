# ==============================================================================
# STAGE 1 · STEP 2 — Neutral-mass recovery (adduct -> monoisotopic neutral mass)
# ------------------------------------------------------------------------------
# Reverses ionization chemistry per feature so masses can be searched against
# metabolite databases (which store neutral monoisotopic masses).
#
# INPUT  : 01_features_cleaned.xlsx   (convert the Step 1 CSV to XLSX, or read CSV)
# OUTPUT : 02_neutral_mass.xlsx       (adds a Neutral_Mass column)
# ==============================================================================

library(dplyr); library(readxl); library(writexl)

# ---- CONFIG ------------------------------------------------------------------
INPUT_FILE  <- "results/01_features_cleaned.xlsx"
OUTPUT_FILE <- "results/02_neutral_mass.xlsx"
# ------------------------------------------------------------------------------

my_data <- read_excel(INPUT_FILE)

calculate_neutral_mass <- function(mz, adduct) {
  p <- 1.007276; na <- 22.989769; k <- 38.963707
  nh4 <- 18.033823; h2o <- 18.010565
  if (is.na(adduct) || trimws(adduct) == "") return(NA_real_)
  switch(adduct,
    "[M+H]+"     = mz - p,
    "[M+Na]+"    = mz - na,
    "[M+K]+"     = mz - k,
    "[M+NH4]+"   = mz - nh4,
    "[M-H]-"     = mz + p,
    "[M]+"       = mz,
    "[M+H-H2O]+" = mz - p + h2o,
    "[M-H2O+H]+" = mz - p + h2o,
    "[M+2H]2+"   = (mz * 2) - (p * 2),
    "[M+H]2+"    = (mz * 2) - p,
    "[2M+H]+"    = (mz - p) / 2,
    NA_real_)   # unrecognised adduct -> flagged, skipped in Step 3
}

my_data$Neutral_Mass <- mapply(calculate_neutral_mass,
                               my_data$`Average Mz`,
                               my_data$`Adduct type`)

write_xlsx(my_data, OUTPUT_FILE)
cat("Unrecognised adducts (Neutral_Mass = NA):",
    sum(is.na(my_data$Neutral_Mass)), "\n")
