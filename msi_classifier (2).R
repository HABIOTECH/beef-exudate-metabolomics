#!/usr/bin/env Rscript

# ==============================================================================
# msi_classifier.R
#
# Assigns Metabolomics Standards Initiative (MSI) confidence levels (2-4) to
# annotated LC-MS/MS features based on retention time, m/z, and MS/MS spectral
# match evidence.
#
# MSI Level 1 is deliberately not assigned: it requires comparison against an
# authentic standard run under identical analytical conditions, which cannot be
# inferred from an annotation table alone.
#
# Usage:
#   Rscript msi_classifier.R <input.xlsx> [output_dir] [sheet]
#
# Example:
#   Rscript msi_classifier.R data/QC_Annotation_merged.xlsx results/ 1
#
# Dependencies: R >= 4.0, readxl, writexl, dplyr
#   install.packages(c("readxl", "writexl", "dplyr"))
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# Classification thresholds
#
# Adjust these to match the acceptance criteria of your study. They are the only
# tunable parameters in the rule set; everything downstream is deterministic.
# ------------------------------------------------------------------------------
DOT_PRODUCT_MIN       <- 80   # minimum spectral dot product for "good" MS/MS
FRAGMENT_PRESENCE_MIN <- 50   # minimum fragment presence (%) for "good" MS/MS

# Values (case-insensitive) treated as TRUE when coercing Excel match flags.
TRUE_TOKENS <- c("TRUE", "T", "1", "YES", "Y")

# Columns the input workbook must contain.
REQUIRED_COLUMNS <- c(
  "Metabolite name", "Formula", "INCHIKEY", "SMILES",
  "RT matched", "m/z matched", "MS/MS matched",
  "Dot product", "Fragment presence %"
)

MSI_LEVELS <- c("MSI Level 2", "MSI Level 3", "MSI Level 4")

# ------------------------------------------------------------------------------
# Arguments
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1L) {
  stop(
    "Usage: Rscript msi_classifier.R <input.xlsx> [output_dir] [sheet]",
    call. = FALSE
  )
}

input_file <- args[1]
output_dir <- if (length(args) >= 2L) args[2] else "."
sheet_name <- if (length(args) >= 3L) args[3] else 1

# Allow a sheet to be given either by index or by name.
if (!is.na(suppressWarnings(as.integer(sheet_name)))) {
  sheet_name <- as.integer(sheet_name)
}

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

#' Coerce a column to logical.
#'
#' Excel round-trips boolean flags inconsistently: the same column may arrive as
#' logical, as numeric 0/1, or as the strings "TRUE" / "T" / "YES". Missing
#' values are treated as FALSE (absence of evidence is not evidence of a match).
to_logical <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  toupper(trimws(as.character(x))) %in% TRUE_TOKENS
}

#' TRUE where a character column carries usable content.
is_present <- function(x, exclude = character(0)) {
  x <- trimws(as.character(x))
  !is.na(x) & x != "" & !(x %in% exclude)
}

#' Abort with a readable message if the input is missing required columns.
check_columns <- function(dat, required) {
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop(
      "Input is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# 1. Read merged annotation table
# ------------------------------------------------------------------------------
message("Reading: ", input_file)
dat <- readxl::read_excel(input_file, sheet = sheet_name)
check_columns(dat, REQUIRED_COLUMNS)
message("Features read: ", nrow(dat))

# ------------------------------------------------------------------------------
# 2. Derive evidence flags
# ------------------------------------------------------------------------------
dat <- dat %>%
  mutate(
    # Normalise match flags to real logicals.
    `RT matched`    = to_logical(`RT matched`),
    `m/z matched`   = to_logical(`m/z matched`),
    `MS/MS matched` = to_logical(`MS/MS matched`),

    # Numeric scores; non-numeric entries become NA and fail the threshold test.
    Dot_product_num       = suppressWarnings(as.numeric(`Dot product`)),
    Fragment_presence_num = suppressWarnings(as.numeric(`Fragment presence %`)),

    # 2a. Any structural information present.
    Has_ID =
      is_present(`Metabolite name`, exclude = "Unknown") |
      is_present(Formula) |
      is_present(INCHIKEY) |
      is_present(SMILES),

    # 2b. Good-quality MS/MS evidence.
    Good_MS2 =
      `MS/MS matched` &
      !is.na(Dot_product_num) & Dot_product_num >= DOT_PRODUCT_MIN &
      !is.na(Fragment_presence_num) &
      Fragment_presence_num >= FRAGMENT_PRESENCE_MIN,

    # 2c. Good RT and m/z match.
    Good_RT_mz = `RT matched` & `m/z matched`
  )

# ------------------------------------------------------------------------------
# 3. MSI level classification
#
#   Level 2  structural ID + RT/m-z match + high-quality MS/MS
#   Level 3  structural ID + RT *or* m/z match, MS/MS evidence insufficient
#   Level 4  everything else (unknown)
# ------------------------------------------------------------------------------
dat <- dat %>%
  mutate(
    MSI_Level = case_when(
      Has_ID & Good_RT_mz & Good_MS2            ~ MSI_LEVELS[1],
      Has_ID & (`m/z matched` | `RT matched`)   ~ MSI_LEVELS[2],
      TRUE                                      ~ MSI_LEVELS[3]
    ),
    MSI_Level = factor(MSI_Level, levels = MSI_LEVELS, ordered = TRUE)
  )

# ------------------------------------------------------------------------------
# 4. Write flat annotated table
# ------------------------------------------------------------------------------
output_flat <- file.path(output_dir, "merged_QC_annotation_MSIclassified.xlsx")
writexl::write_xlsx(dat, output_flat)
message("Wrote: ", normalizePath(output_flat))

# ------------------------------------------------------------------------------
# 5. Write multi-sheet workbook split by MSI level
# ------------------------------------------------------------------------------

# Sort each level by strength of evidence, then retention time, where available.
sort_by_evidence <- function(x) {
  if (all(c("Dot_product_num", "RT") %in% names(x))) {
    dplyr::arrange(x, dplyr::desc(Dot_product_num), RT)
  } else if ("Dot_product_num" %in% names(x)) {
    dplyr::arrange(x, dplyr::desc(Dot_product_num))
  } else {
    x
  }
}

split_list <- lapply(split(dat, dat$MSI_Level, drop = TRUE), sort_by_evidence)
names(split_list) <- gsub(" ", "_", names(split_list))   # -> "MSI_Level_2"

sheets_to_write <- c(list(All = dat), split_list)

output_levels <- file.path(
  output_dir, "merged_QC_annotation_MSIclassified_byLevel.xlsx"
)
writexl::write_xlsx(sheets_to_write, path = output_levels)
message("Wrote: ", normalizePath(output_levels))

# ------------------------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------------------------
counts <- table(dat$MSI_Level, useNA = "ifany")

cat("\nMSI level counts:\n")
print(counts)
cat("\nTotal features: ", nrow(dat), "\n", sep = "")
cat(sprintf(
  "Annotated (Level 2 or 3): %d (%.1f%%)\n",
  sum(counts[MSI_LEVELS[1:2]]),
  100 * sum(counts[MSI_LEVELS[1:2]]) / nrow(dat)
))
