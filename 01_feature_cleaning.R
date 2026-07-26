# ==============================================================================
# STAGE 1 · STEP 1 — Feature cleaning
# ------------------------------------------------------------------------------
# Removes chemical noise from a raw positive-mode LC-MS feature table:
#   (1) contaminants  (2) redundant Na/K/NH4 adducts + 13C isotopes  (3) polymers
# All filters use a 10 ppm mass tolerance; adduct grouping adds a 0.1 min RT window.
#
# INPUT  : feature_table.xlsx   (peak-picker export: MS-DIAL / MZmine)
# OUTPUT : 01_features_cleaned.csv   (every feature tagged with a Status column)
# ==============================================================================

library(dplyr); library(readxl); library(readr)

# ---- CONFIG ------------------------------------------------------------------
INPUT_FILE  <- "data/feature_table.xlsx"
OUTPUT_FILE <- "results/01_features_cleaned.csv"
PPM_TOL     <- 10     # mass-accuracy tolerance (ppm)
RT_TOL      <- 0.1    # retention-time tolerance (min) for adduct grouping
# ------------------------------------------------------------------------------

my_data <- read_excel(INPUT_FILE) %>% arrange(`Average Mz`)

# --- Filter 1: contaminants ---------------------------------------------------
contaminants <- read_tsv(
  "https://raw.githubusercontent.com/stanstrup/commonMZ/master/inst/extdata/contaminants_%2B.tsv",
  show_col_types = FALSE)
my_data$Is_Contaminant <- sapply(my_data$`Average Mz`, function(mz)
  any(abs(contaminants$mz - mz) / mz * 1e6 <= PPM_TOL, na.rm = TRUE))

# --- Filter 2: redundant adducts and 13C isotopes -----------------------------
adduct_diffs <- c(Na = 21.9819, K = 37.9559, NH4 = 17.0265, C13 = 1.0033)
my_data$Is_Adduct <- FALSE
for (i in seq_len(nrow(my_data) - 1)) {
  if (my_data$Is_Adduct[i]) next
  mz1 <- my_data$`Average Mz`[i]; rt1 <- my_data$`Average Rt(min)`[i]
  coel <- which(abs(my_data$`Average Rt(min)` - rt1) <= RT_TOL &
                  seq_len(nrow(my_data)) > i)
  for (idx in coel) {
    mz2 <- my_data$`Average Mz`[idx]
    if (any(abs((mz2 - mz1) - adduct_diffs) <= PPM_TOL / 1e6 * mz2))
      my_data$Is_Adduct[idx] <- TRUE
  }
}

# --- Filter 3: polymers / repeating units -------------------------------------
polymers <- read_tsv(
  "https://raw.githubusercontent.com/stanstrup/commonMZ/master/inst/extdata/repeating_units_%2B.tsv",
  show_col_types = FALSE)
my_data$Is_Polymer <- FALSE
for (i in seq_len(nrow(my_data))) {
  if (my_data$Is_Polymer[i] || my_data$Is_Contaminant[i]) next
  mz1 <- my_data$`Average Mz`[i]
  for (p_mass in polymers$mass) {
    hit <- which(abs(my_data$`Average Mz` - (mz1 + p_mass)) /
                   (mz1 + p_mass) * 1e6 <= PPM_TOL)
    if (length(hit)) my_data$Is_Polymer[hit] <- TRUE
  }
}

# --- Unified Status column and export -----------------------------------------
my_data <- my_data %>% mutate(Status = case_when(
  Is_Contaminant ~ "Contaminant",
  Is_Polymer     ~ "Polymer",
  Is_Adduct      ~ "Redundant Adduct/Isotope",
  TRUE           ~ "Clean Biological"))

write_csv(my_data, OUTPUT_FILE)
cat("Status counts:\n"); print(table(my_data$Status))
