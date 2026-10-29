###############################################################################
# tools/make_minimal_dataset.R
# Creates the data files for sharing: same rows as the full study files, but
# only the columns used by the analysis code (listed in data/README.md). All
# other columns (e.g. province, job, study degree, other clinical or
# neuropsychological scores) are dropped.
# Run from the repository root, with the full study files in data/:
#   source("tools/make_minimal_dataset.R")
# Output: data_minimal/data_adult.xlsx, data_minimal/data_psy.xlsx
###############################################################################

suppressPackageStartupMessages({ library(readxl); library(writexl); library(dplyr) })

IN_DIR  <- Sys.getenv("MATRIKS_DATA_DIR", unset = "data")
OUT_DIR <- "data_minimal"
dir.create(OUT_DIR, showWarnings = FALSE)

CONTROL_ITEMS <- c(4, 13, 14, 15, 16, 27, 30, 33, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47)
PATIENT_ITEMS <- 34:51

item_columns <- function(version, items) {
  c(paste0("response_mat", version, "_item.", items), paste0("time_mat", version, "_item.", items))
}

ADULT_COLUMNS <- c("new_id", "age", "gender", "diagnostic", "anni_scolarita",
                   item_columns(50, CONTROL_ITEMS))
PSY_COLUMNS   <- c("external_code", "Diagnosi", "age", "gender", "anni_scolarita",
                   "BPRS_PRE", "BPRS_POST", "HONOS_PRE", "HONOS_POST",
                   "Matrici_Attentive", "SPAN_diretto", "SPAN_inverso",
                   item_columns(88, PATIENT_ITEMS))

make_minimal <- function(file, columns) {
  d <- read_xlsx(file.path(IN_DIR, file))
  missing <- setdiff(columns, names(d))
  if (length(missing) > 0) stop(file, ": missing columns ", paste(missing, collapse = ", "), call. = FALSE)
  message(file, ": ", ncol(d), " -> ", length(columns), " columns; dropped: ",
          paste(setdiff(names(d), columns)[!grepl("_item\\.", setdiff(names(d), columns))], collapse = ", "),
          " (and the item columns not used by the code)")
  write_xlsx(select(d, all_of(columns)), file.path(OUT_DIR, file))
}

make_minimal("data_adult.xlsx", ADULT_COLUMNS)
make_minimal("data_psy.xlsx",   PSY_COLUMNS)
