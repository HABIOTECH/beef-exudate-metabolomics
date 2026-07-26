# LC-MS Untargeted Metabolomics Annotation Pipeline

**Standard Operating Procedure**

> This document is the prose SOP for the annotation pipeline. The R code that
> implements each step lives alongside it in the repository as separate scripts;
> each step below points to the file that runs it rather than embedding the code.
>
> **Script filenames referenced below are placeholders** — adjust them to match
> the files you uploaded to the repo.

---

## 1. Overview

This document describes the end-to-end computational pipeline used in the Homics
lab for annotating untargeted LC-MS metabolomics features acquired in positive
ionization mode on high-resolution mass spectrometers (Orbitrap, Q-TOF). The
pipeline takes a raw feature table exported from peak-picking software (e.g.
MS-DIAL, MZmine) and returns a curated, biologically meaningful identification
table linked to HMDB and IDEOM metadata.

It serves two purposes: as an internal Standard Operating Procedure to
standardize how projects are handled by the team, and as documentation for the R
code base that implements each step. Every stage specifies its inputs, outputs,
parameters, and scientific rationale.

### 1.1 Pipeline at a glance

| Step | Stage | Purpose | Primary output | Script |
|---|---|---|---|---|
| 1 | Feature cleaning | Remove contaminants, redundant adducts/isotopes, polymers | `Combined_Annotated_Metabolites.csv` | `R/01_feature_cleaning.R` |
| 2 | Neutral mass recovery | Convert observed *m/z* to monoisotopic neutral mass | `Data_with_Neutral_Mass.xlsx` | `R/02_neutral_mass.R` |
| 3 | Database matching | Match neutral mass against HMDB and IDEOM at 10 ppm | `Best_Identified_Metabolites_with_IDs.xlsx` | `R/03a_hmdb_streaming_parser.R`, `R/03b_database_matching.R` |
| 4 | Best-hit & isomer handling | Select lowest-ppm hit; merge isobaric isomers | *(same workbook)* | *(embedded in Step 3)* |
| 5 | Batch annotation | Apply identifications across all experimental cohorts | `Annotated_Results/` (folder) | `R/05_batch_annotation.R` |
| 6 | Golden List extraction | Filter to statistically significant features with metadata | `Extracted_Target_Metabolites.xlsx` | `R/06_golden_list.R` |

### 1.2 Input data requirements

The pipeline assumes a single feature table per cohort (one `.xlsx` file)
exported from the peak-picking tool. The following columns are required — names
must match exactly, as they are referenced literally in the code:

- **Alignment ID** — unique integer identifier for each feature.
- **Average Mz** — observed mass-to-charge ratio.
- **Average Rt(min)** — retention time in minutes.
- **Adduct type** — adduct assignment string, e.g. `[M+H]+`, `[M+Na]+`, `[M+NH4]+`.

Sample-intensity columns are retained unchanged throughout the pipeline and
re-merged with identifications at the reporting stage.

### 1.3 Software environment

| | |
|---|---|
| **Language** | R (≥ 4.2.0) |
| **CRAN packages** | `dplyr`, `tidyr`, `readr`, `readxl`, `writexl`, `xml2`, `stringr` |
| **External DBs** | HMDB (`hmdb_metabolites.xml`, ~6 GB uncompressed) and IDEOM reference file (`DB.xlsx`) |
| **Contaminant DBs** | `stanstrup/commonMZ` (`contaminants_+.tsv`, `repeating_units_+.tsv`), fetched live from GitHub |
| **RAM** | ≥ 8 GB recommended (Step 3 uses a streaming XML parser to keep memory below 1 GB) |

Install the CRAN dependencies in one call:

```r
install.packages(c("dplyr", "tidyr", "readr", "readxl", "writexl", "xml2", "stringr"))
```

---

## 2. Step 1 — Feature cleaning

The raw feature table is filtered in three sequential passes to remove chemical
noise that would otherwise waste downstream database queries and inflate the
false-positive rate. All three filters use the same 10 ppm mass-accuracy
tolerance; the adduct/isotope filter adds a 0.1-minute retention-time window to
require co-elution.

> **Implemented in `R/01_feature_cleaning.R`.**

### 2.1 Filters applied

| Filter | What it removes | Reference / mass differences |
|---|---|---|
| Contaminants | Plasticizers, solvents, common MS background ions | `stanstrup/commonMZ` `contaminants_+.tsv` |
| Adducts & isotopes | Na⁺ / K⁺ / NH₄⁺ adducts and ¹³C isotopologues of a feature already in the table | Δ*m/z* = 21.9819 (Na), 37.9559 (K), 17.0265 (NH₄), 1.0033 (¹³C) |
| Polymers | PEG/PPG and other repeating-unit series | `stanstrup/commonMZ` `repeating_units_+.tsv` |

### 2.2 Rationale

In positive-mode LC-MS a single true metabolite routinely generates several
features: the protonated parent `[M+H]+`, sodium/potassium/ammonium adducts, and
the first ¹³C isotopologue. Searching every one of these against a metabolite
database produces duplicated hits and inflates apparent coverage. By requiring
co-elution (RT within 0.1 min) before flagging a pair as an adduct relationship,
the filter avoids collapsing genuinely different metabolites that happen to
differ by a sodium mass.

### 2.3 Inputs and outputs

| | |
|---|---|
| **Input** | `Metabolites_database.xlsx` (feature table from peak picker) |
| **Output** | `Combined_Annotated_Metabolites.csv` (all features tagged with a `Status` column) |
| **Status values** | `Clean Biological` \| `Contaminant` \| `Polymer` \| `Redundant Adduct/Isotope` |

---

## 3. Step 2 — Neutral mass recovery

Mass spectrometers detect charged ions, but metabolite databases store neutral
monoisotopic masses. This step reverses the ionization chemistry for every
feature, producing the neutral mass that can legitimately be searched against
HMDB and IDEOM.

> **Implemented in `R/02_neutral_mass.R`.**

### 3.1 Supported adduct chemistries

| Adduct | Formula applied | Notes |
|---|---|---|
| `[M+H]+` | M = mz − 1.007276 | Default assumption if Adduct type is blank. |
| `[M+Na]+` | M = mz − 22.989769 | Sodium adduct. |
| `[M+K]+` | M = mz − 38.963707 | Potassium adduct. |
| `[M+NH4]+` | M = mz − 18.033823 | Ammonium adduct. |
| `[M-H]−` | M = mz + 1.007276 | Present for mixed-polarity files. |
| `[M]+` | M = mz | Radical cation (no proton transfer). |
| `[M+H-H2O]+` | M = mz − 1.007276 + 18.010565 | In-source water loss. |
| `[M+2H]2+` | M = 2·mz − 2·1.007276 | Doubly protonated. |
| `[M+H]2+` | M = 2·mz − 1.007276 | Doubly charged, single proton. |
| `[2M+H]+` | M = (mz − 1.007276) / 2 | Dimer. |
| Unknown / blank | `"Not Identified"` | Feature flagged and skipped in Step 3. |

### 3.2 Inputs and outputs

| | |
|---|---|
| **Input** | `Combined_Annotated_Metabolites.xlsx` (converted from the Step 1 CSV) |
| **Output** | `Data_with_Neutral_Mass.xlsx` (new column: `Neutral_Mass`) |

---

## 4. Step 3 — Database matching (HMDB + IDEOM)

Two orthogonal reference databases are queried at 10 ppm: HMDB for breadth of
human-relevant metabolites, and IDEOM (provided in-house as `DB.xlsx`), which
couples identifications to formulas, KEGG pathways, and cross-references.
Searching both increases coverage while preserving the richer metadata IDEOM
provides.

> **Implemented in two parts:**
> `R/03a_hmdb_streaming_parser.R` (one-time HMDB parse) and
> `R/03b_database_matching.R` (matching against HMDB + IDEOM).

### 4.1 HMDB preparation — memory-safe streaming parser

The full HMDB XML is ~6 GB unpacked. Loading it with `xml2::read_xml()` would
exceed typical laptop memory. The parser streams the file line-by-line, extracts
only the three fields the pipeline needs, and writes a ~15 MB CSV
(`HMDB_Slim.csv`) that is reused on all subsequent runs.

### 4.2 Parameters

| | |
|---|---|
| **Mass accuracy** | 10 ppm (standard for Orbitrap / Q-TOF in positive mode) |
| **RT filtering** | Bypassed at the global matching step. See Design Note §8.1. |
| **Scope** | Only features with `Status = "Clean Biological"` are matched; contaminants, adducts and polymers are skipped. |
| **Source fields (HMDB)** | `accession`, `name`, `monoisotopic_molecular_weight` |
| **Source fields (IDEOM)** | `searchmass`, `Metabolite`, `Formula`, `Map`, `DB`, `Pathway`, KEGG/HMDB/CAS/CHEBI/METLIN/LMID IDs, `InchiKEY`, `SMILES`, `Synonyms` |

### 4.3 Inputs and outputs

| | |
|---|---|
| **Inputs** | `hmdb_metabolites.xml` (downloaded once from hmdb.ca) + `DB.xlsx` (IDEOM) + `Data_with_Neutral_Mass.xlsx` |
| **Derived** | `HMDB_Slim.csv` (cached on first run) |
| **Output** | `Best_Identified_Metabolites_with_IDs.xlsx` (one best hit per feature) |

---

## 5. Step 4 — Best-hit selection and isomer handling

Within the 10 ppm window a single feature frequently matches several database
entries. Two cases must be distinguished:

- **Ranked candidates** — different metabolites with different ppm errors. The
  candidate with the lowest ppm error is retained as the best hit. The ppm error
  itself is written to a `PPM_Error` column for downstream auditing.
- **True isobaric isomers** — different metabolites with identical exact mass.
  Because MS1 alone cannot resolve them, their names are concatenated into a
  single cell as `"name A ; name B"` so no identification is silently discarded.

This logic is embedded in the Step 3 matching script (arrange by ppm error, then
paste unique names). No separate script is required.

---

## 6. Step 5 — Batch annotation across cohorts

Production projects contain multiple experimental groups. In the current
workflow these are: *E. coli* tissue, *E. coli* exudate, *Salmonella* tissue,
and a fourth exudate cohort (file `Data Exudate SA.xlsx`). The identification
logic from Steps 2–4 is applied to each cohort file in a single loop, producing
one dual-sheet workbook per cohort inside an `Annotated_Results/` folder.

> **⚠️ Verify cohort name:** the source document labels the fourth cohort
> "*Salmonella aureus* exudate," which is not a valid organism name — it appears
> to conflate *Salmonella* with *Staphylococcus aureus*. The filename `SA` is
> ambiguous. Confirm which organism this cohort is and correct the label before
> publishing.

> **Implemented in `R/05_batch_annotation.R`.**

### 6.1 Workbook structure

- **Sheet 1 — Clean Data:** original feature table plus a `Metabolite Name`
  column merged in on `Alignment ID`.
- **Sheet 2 — Identification Details:** identified features only, joined to
  IDEOM metadata (Formula, KEGG pathway, HMDB/KEGG/CAS/CHEBI/METLIN/LMID IDs,
  InChIKey, SMILES, synonyms, taxonomy, category, species, alt-IDs).

---

## 7. Step 6 — Golden List extraction

After downstream statistical analysis (e.g. MetaboAnalyst, limma, or PLS-DA VIP
filtering), the analyst produces a short list of Alignment IDs deemed
statistically significant. This final step pulls the full metadata for exactly
those IDs, so the paper's results table, pathway map, and supplementary
materials can be assembled from a single curated workbook.

> **Implemented in `R/06_golden_list.R`.**

### 7.1 Inputs and outputs

| | |
|---|---|
| **Inputs** | `Best_Identified_Metabolites_with_IDs.xlsx` + `Target_IDs.xlsx` (one column: `Alignment`) |
| **Output** | `Extracted_Target_Metabolites.xlsx` (for the results table) and optionally `Extracted_HMDB_IDs.xlsx` (slim, for pathway tools) |

---

## 8. Design notes and quality-control guidance

### 8.1 Why retention time is bypassed at the matching stage

RT is cohort- and column-dependent and drifts between batches. Enforcing an RT
window at the global MS1-matching step would silently discard true positives
every time a new LC method is used. The pipeline therefore matches on exact mass
only and relies on the ppm-ranking and isomer-merging logic of Step 4 to control
ambiguity. RT is still used defensively inside Step 1 (co-elution requirement for
adduct grouping), where it strengthens rather than restricts the logic.

### 8.2 Why two databases are queried

HMDB contributes breadth — tens of thousands of biologically observed compounds
across human biofluids and tissues. IDEOM contributes depth — each record ships
with the KEGG pathway, cross-references, InChIKey and SMILES required for pathway
enrichment. Running both in parallel and reporting both columns lets the
biologist pick whichever identifier their downstream tool needs, without
re-running the pipeline.

### 8.3 Quality-control checkpoints

- **After Step 1:** contaminants typically make up 5–15% of features. If the
  figure is > 25%, suspect solvent or vial carry-over.
- **After Step 2:** inspect any rows where `Neutral_Mass` is `NA` — these
  indicate adduct strings the calculator did not recognise, which should be added
  to the adduct-handling logic before re-running.
- **After Step 3:** identification rates of 30–60% of Clean Biological features
  are typical for untargeted positive-mode work. A sudden drop indicates an
  instrument calibration drift beyond 10 ppm.
- **After Step 5:** open each `Annotated_*.xlsx` and confirm the Identification
  Details sheet row count equals the number of Clean Biological IDs with non-NA
  names.

### 8.4 File conventions

All intermediate and final files live in the project working directory. CSV is
used where downstream tools (MetaboAnalyst, Python scripts) benefit from plain
text; XLSX is used wherever the team will inspect or edit the file manually. Raw
cohort filenames are prefixed with `Annotated_` rather than renamed, so the link
back to the peak-picker output is always obvious.

---

## 9. Change log

| Version | Date | Summary |
|---|---|---|
| 1.0 | 2026-04-19 | Initial consolidated pipeline. Supersedes the previously separate R scripts (contaminant-only filter, combined filters, combined filters with Status export, HMDB direct parser, HMDB streaming parser, and the batch-annotation and Golden-List scripts). |

*— End of document —*
