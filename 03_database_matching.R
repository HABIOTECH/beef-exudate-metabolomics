# ==============================================================================
# STAGE 1 · STEP 3 — Database matching (HMDB + IDEOM, 10 ppm)
# ------------------------------------------------------------------------------
# Part A : one-time, memory-safe streaming parse of the ~6 GB HMDB XML into a
#          compact (~15 MB) CSV that is cached and reused on every later run.
# Part B : match each Clean Biological feature's neutral mass against HMDB and
#          IDEOM at 10 ppm; IDEOM is primary (richer metadata), HMDB is back-up.
#
# INPUTS : hmdb_metabolites.xml (download once from hmdb.ca), ideom_reference.xlsx,
#          02_neutral_mass.xlsx
# OUTPUT : 03_identified.xlsx   (best hits + ppm errors per feature)
# ==============================================================================

# ---- CONFIG ------------------------------------------------------------------
HMDB_XML     <- "data/hmdb_metabolites.xml"
HMDB_SLIM    <- "data/HMDB_Slim.csv"          # cached after first run
IDEOM_FILE   <- "data/ideom_reference.xlsx"   # in-house IDEOM reference
INPUT_FILE   <- "results/02_neutral_mass.xlsx"
OUTPUT_FILE  <- "results/03_identified.xlsx"
PPM_TOL      <- 10
# ------------------------------------------------------------------------------

# ============================ PART A — HMDB stream ============================
if (!file.exists(HMDB_SLIM)) {
  con_in  <- file(HMDB_XML, "r")
  con_out <- file(HMDB_SLIM, "w")
  writeLines("Accession,Name,Neutral_Mass", con_out)

  acc <- ""; metab_name <- ""; mass <- ""; count <- 0
  repeat {
    lines <- readLines(con_in, n = 50000, warn = FALSE)
    if (length(lines) == 0) break
    for (line in lines) {
      if (grepl("^  <metabolite>", line)) {
        acc <- ""; metab_name <- ""; mass <- ""
      } else if (grepl("^    <accession>", line)) {
        acc <- sub(".*<accession>(.*)</accession>.*", "\\1", line)
      } else if (grepl("^    <name>", line)) {
        raw <- sub(".*<name>(.*)</name>.*", "\\1", line)
        metab_name <- paste0('"', gsub('"', '""', raw), '"')
      } else if (grepl("^    <monoisotopic_molecular_weight>", line)) {
        mass <- sub(".*<monoisotopic_molecular_weight>(.*)</monoisotopic_molecular_weight>.*",
                    "\\1", line)
      } else if (grepl("^  </metabolite>", line)) {
        if (nzchar(acc) && nzchar(mass)) {
          writeLines(paste(acc, metab_name, mass, sep = ","), con_out)
          count <- count + 1
        }
      }
    }
  }
  close(con_in); close(con_out)
  message("Parsed ", count, " HMDB metabolites to ", HMDB_SLIM)
}

# ============================ PART B — matching ==============================
library(readxl); library(readr); library(writexl); library(dplyr)

exp_data <- read_excel(INPUT_FILE)
hmdb_db  <- read_csv(HMDB_SLIM, show_col_types = FALSE) %>%
            filter(!is.na(Neutral_Mass))
ideom_db <- read_excel(IDEOM_FILE)

results <- vector("list", nrow(exp_data))

for (i in seq_len(nrow(exp_data))) {
  if (exp_data$Status[i] != "Clean Biological") next
  mz <- exp_data$`Average Mz`[i]
  nm <- exp_data$Neutral_Mass[i]
  id <- exp_data$`Alignment ID`[i]
  if (is.na(mz) || is.na(nm)) next

  # IDEOM (primary)
  ideom_hits <- ideom_db %>%
    mutate(ppm_err = abs(searchmass - mz) / searchmass * 1e6) %>%
    filter(ppm_err <= PPM_TOL) %>% arrange(ppm_err)

  # HMDB (secondary)
  hmdb_hits <- hmdb_db %>%
    mutate(ppm_err = abs(Neutral_Mass - nm) / nm * 1e6) %>%
    filter(ppm_err <= PPM_TOL) %>% arrange(ppm_err)

  if (nrow(ideom_hits) > 0 || nrow(hmdb_hits) > 0) {
    results[[i]] <- tibble(
      Alignment_ID = id,
      Exp_Mz       = mz,
      Exp_Neutral  = nm,
      IDEOM_Name   = if (nrow(ideom_hits)) paste(unique(ideom_hits$Metabolite), collapse = " ; ") else NA,
      IDEOM_PPM    = if (nrow(ideom_hits)) round(min(ideom_hits$ppm_err), 2) else NA,
      HMDB_IDs     = if (nrow(hmdb_hits))  paste(hmdb_hits$Accession, collapse = " | ") else NA,
      HMDB_Names   = if (nrow(hmdb_hits))  paste(hmdb_hits$Name, collapse = " | ") else NA)
  }
}

write_xlsx(bind_rows(results), OUTPUT_FILE)
cat("Identified features written to", OUTPUT_FILE, "\n")
