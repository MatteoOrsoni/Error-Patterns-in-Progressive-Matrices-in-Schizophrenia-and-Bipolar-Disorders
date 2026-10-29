###############################################################################
# 00_setup.R
# Shared configuration for all analysis scripts:
#   - package checks and loading
#   - file paths (relative to the repository root)
#   - sample definitions (exclusions, item mapping, response labels)
#   - MCMC settings and random seed
#   - plotting palette/theme
#   - helper functions used by more than one script
#
# Every script in R/ sources this file. Run all scripts from the repository
# root (the folder containing run_all.R).
###############################################################################

if (!file.exists(file.path("R", "00_setup.R"))) {
  stop("Run the scripts from the repository root (the folder that contains ",
       "'run_all.R'), e.g. setwd('<path>/matriks-error-patterns').", call. = FALSE)
}

# ---- Packages ----------------------------------------------------------------
require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         ". Run source('install_dependencies.R') first.", call. = FALSE)
  }
  invisible(TRUE)
}

CORE_PACKAGES  <- c("readxl", "writexl", "dplyr", "tidyr", "purrr", "tibble",
                    "stringr", "ggplot2", "patchwork", "emmeans", "broom", "clue")
BAYES_PACKAGES <- c("cmdstanr", "posterior", "bayesplot")   # checked in scripts 03-05

require_packages(CORE_PACKAGES)
suppressPackageStartupMessages({
  library(readxl);  library(writexl); library(dplyr);   library(tidyr)
  library(purrr);   library(tibble);  library(stringr); library(ggplot2)
  library(patchwork); library(emmeans); library(broom)
})

# ---- Paths -------------------------------------------------------------------
# Defaults can be overridden with the environment variables MATRIKS_DATA_DIR and
# MATRIKS_OUTPUT_DIR (e.g. to run the pipeline on the simulated example data
# written by tools/simulate_example_data.R).
DATA_DIR   <- Sys.getenv("MATRIKS_DATA_DIR",   unset = "data")
OUTPUT_DIR <- Sys.getenv("MATRIKS_OUTPUT_DIR", unset = "outputs")
ADULT_FILE <- file.path(DATA_DIR, "data_adult.xlsx")   # healthy controls
PSY_FILE   <- file.path(DATA_DIR, "data_psy.xlsx")     # psychiatric inpatients
STAN_FILE  <- file.path("stan", "multinomial_exchangeable.stan")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Folder for CmdStan's CSV output. It must contain ASCII characters only:
# CmdStan writes this path into the CSV header in the system encoding, and
# cmdstanr fails to read the files back if the path contains e.g. accented
# letters. Override with the environment variable MATRIKS_STAN_CSV_DIR.
STAN_CSV_DIR <- Sys.getenv("MATRIKS_STAN_CSV_DIR",
                           unset = file.path(tempdir(), "matriks_stan_csv"))

# ---- Reproducibility ----------------------------------------------------------
SEED <- 42
set.seed(SEED)

# MCMC settings used in the paper (4 chains x 2,000 warm-up + 4,000 sampling).
MCMC <- list(chains = 4, parallel_chains = 4, iter_warmup = 2000,
             iter_sampling = 4000, adapt_delta = 0.95, max_treedepth = 12)

# Quick smoke test (NOT the settings of the paper): set the environment
# variable MATRIKS_QUICK_TEST=1 to check that the pipeline runs end to end.
QUICK_TEST <- identical(Sys.getenv("MATRIKS_QUICK_TEST"), "1")
if (QUICK_TEST) {
  MCMC$chains <- 2; MCMC$parallel_chains <- 2
  MCMC$iter_warmup <- 150; MCMC$iter_sampling <- 150
  message("QUICK TEST MODE: reduced MCMC settings; results will NOT match the paper.")
}

# ---- Sample definitions ---------------------------------------------------------
# Psychiatric participants excluded from all analyses (N = 62 in the raw file,
# N = 59 analysed). The reason for exclusion is documented in README.md.
EXCLUDED_PATIENT_CODES <- c("BSPDC40", "BSPDC41", "BSPDC45")

# The 18 items shared by the two MatriKS versions, in corresponding order:
# item k of CONTROL_ITEMS (MatriKS 12+, columns response_mat50_item.<n>) is the
# same stimulus as item k of PATIENT_ITEMS (MatriKS88, response_mat88_item.<n>).
CONTROL_ITEMS <- c(4, 13, 14, 15, 16, 27, 30, 33, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47)
PATIENT_ITEMS <- 34:51
stopifnot(length(CONTROL_ITEMS) == 18, length(PATIENT_ITEMS) == 18)

# Sensitivity analysis (not the main analysis): MatriKS88 item numbers (34-51)
# to exclude from the Bayesian models, e.g. Sys.setenv(MATRIKS_EXCLUDE_ITEMS = "47").
# Empty (default) = main analysis on all 18 items.
MODEL_EXCLUDED_ITEMS <- suppressWarnings(as.integer(strsplit(Sys.getenv("MATRIKS_EXCLUDE_ITEMS", ""), ",")[[1]]))
if (anyNA(MODEL_EXCLUDED_ITEMS) || !all(MODEL_EXCLUDED_ITEMS %in% PATIENT_ITEMS)) {
  stop("MATRIKS_EXCLUDE_ITEMS must contain MatriKS88 item numbers between 34 and 51.", call. = FALSE)
}

CONTROL_ITEM_STUB <- "response_mat50_item."
PATIENT_ITEM_STUB <- "response_mat88_item."
item_regex <- function(stub, items) {
  paste0("^", gsub(".", "\\.", stub, fixed = TRUE), "(", paste(items, collapse = "|"), ")$")
}
CONTROL_ITEM_REGEX <- item_regex(CONTROL_ITEM_STUB, CONTROL_ITEMS)
PATIENT_ITEM_REGEX <- item_regex(PATIENT_ITEM_STUB, PATIENT_ITEMS)

# Response labels recorded by MatriKS for each item.
ERROR_LABELS <- list(
  incomplete_correlate = c("ic.flip", "ic.inc", "ic.neg", "ic.scale"),
  repetition           = c("r.diag", "r.left", "r.top"),
  wrong_principle      = c("wp.copy", "wp.matrix", "wp1"),
  difference           = c("d.union", "diff", "diff1", "diff2"),
  r_ic                 = c("r.ic"),
  correct              = c("correct")
)
NON_RESPONSE_LABELS <- c("skip")

# Error type assigned to "r.ic" responses (the label MatriKS records, in three
# items, for the choice of either repetition distractor; see README.md).
#   "R"    main analysis (paper): counted as repetition errors
#   "IC"   sensitivity analysis (Appendix I): counted as incomplete correlate errors
#   "drop" excluded from the error-type analyses
# Set with the environment variable MATRIKS_RIC_POLICY (default "R"); see
# run_sensitivity_analyses.R.
RIC_POLICY <- Sys.getenv("MATRIKS_RIC_POLICY", unset = "R")
if (!RIC_POLICY %in% c("R", "IC", "drop")) {
  stop("MATRIKS_RIC_POLICY must be \"R\", \"IC\" or \"drop\".", call. = FALSE)
}
if (RIC_POLICY != "R") message("Sensitivity analysis: \"r.ic\" responses classified as ", RIC_POLICY, ".")

# Performance matching (scripts 03 and 07): optimal matching in triplets on
# accuracy over the answered common items. Pairs further apart than
# MATCH_CALIPER (0.10, about two of the 18 items) are not allowed, and a triplet
# is kept only if the accuracies of its three members differ by at most the
# caliper. MATCH_AGE_WEIGHT only breaks ties between candidates with the same
# accuracy (closest age first): its largest possible contribution to the total
# distance is far below the smallest non-zero difference in squared accuracy.
MATCH_CALIPER    <- 0.10
MATCH_AGE_WEIGHT <- 1e-9

# ---- Plotting ------------------------------------------------------------------
DX_COLORS <- c("Bipolar" = "#84994F", "Control" = "#A72703", "Schizophrenic" = "#FCB53B")

DX_THEME <- theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 18),
    axis.title   = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold")
  )

# ---- Helper functions -----------------------------------------------------------

# Numeric conversion that also accepts a comma as decimal separator.
to_num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x), fixed = TRUE)))

# Years of education are stored as text, with "NULL" for missing values.
to_years <- function(x) as.numeric(na_if(as.character(x), "NULL"))

check_columns <- function(data, cols, file) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop(sprintf("%s: missing columns: %s", file, paste(missing, collapse = ", ")), call. = FALSE)
  }
}

check_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Data file not found: '%s'. See data/README.md.", path), call. = FALSE)
  }
}

# Healthy controls: adults (>= 18 years) without any recorded diagnosis.
# Rows with gender "B" (and rows with missing gender) are excluded.
load_controls <- function(path = ADULT_FILE) {
  check_file(path)
  d <- read_xlsx(path)
  check_columns(d, c("new_id", "age", "gender", "diagnostic", "anni_scolarita",
                     paste0(CONTROL_ITEM_STUB, CONTROL_ITEMS)), path)
  d %>%
    filter(age >= 18, is.na(diagnostic), gender != "B") %>%
    mutate(
      external_code  = as.character(new_id),
      Diagnosi       = "Control",
      anni_scolarita = to_years(anni_scolarita)
    )
}

# Psychiatric inpatients (Bipolar / Schizophrenic), after exclusions.
load_patients <- function(path = PSY_FILE) {
  check_file(path)
  d <- read_excel(path)
  check_columns(d, c("external_code", "Diagnosi", "age", "gender", "anni_scolarita",
                     "BPRS_PRE", "BPRS_POST", "HONOS_PRE", "HONOS_POST",
                     paste0(PATIENT_ITEM_STUB, PATIENT_ITEMS)), path)
  d %>%
    filter(!external_code %in% EXCLUDED_PATIENT_CODES) %>%
    mutate(
      external_code  = as.character(external_code),
      BPRS_PRE       = as.numeric(as.character(BPRS_PRE)),
      BPRS_POST      = as.numeric(as.character(BPRS_POST)),
      HONOS_PRE      = as.numeric(as.character(HONOS_PRE)),
      HONOS_POST     = as.numeric(as.character(HONOS_POST)),
      anni_scolarita = to_years(anni_scolarita)
    )
}

# Warn about response labels that are neither known error/correct labels nor
# non-responses: such responses would be silently ignored by the counts below.
check_response_labels <- function(data, item_regex, dataset_label) {
  known  <- c(unlist(ERROR_LABELS), NON_RESPONSE_LABELS)
  values <- unique(unlist(lapply(select(data, matches(item_regex)), as.character)))
  values <- values[!is.na(values)]
  unknown <- setdiff(values, known)
  if (length(unknown) > 0) {
    warning(sprintf("%s: unrecognised response labels (not counted): %s",
                    dataset_label, paste(unknown, collapse = ", ")), call. = FALSE)
  }
  invisible(unknown)
}

# Per-participant counts of correct responses and of each error type on the
# selected items. Skipped items count neither as correct nor as errors.
# "r.ic" responses are counted separately (r_ic) and included in total_errors.
add_error_scores <- function(df, item_regex) {
  count_labels <- function(row, labels) sum(row %in% labels, na.rm = TRUE)
  df %>%
    rowwise() %>%
    mutate(
      incomplete_correlate = count_labels(c_across(matches(item_regex)), ERROR_LABELS$incomplete_correlate),
      repetition           = count_labels(c_across(matches(item_regex)), ERROR_LABELS$repetition),
      wrong_principle      = count_labels(c_across(matches(item_regex)), ERROR_LABELS$wrong_principle),
      difference           = count_labels(c_across(matches(item_regex)), ERROR_LABELS$difference),
      r_ic                 = count_labels(c_across(matches(item_regex)), ERROR_LABELS$r_ic),
      total_matriks        = count_labels(c_across(matches(item_regex)), ERROR_LABELS$correct)
    ) %>%
    ungroup() %>%
    mutate(total_errors = incomplete_correlate + repetition + wrong_principle + difference + r_ic)
}

# Combined participant-level dataset (patients + controls) scored on the
# 18 common items; used by the ANOVAs, the error-type proportions and the
# symptom-severity regression.
build_scored_dataset <- function() {
  controls <- load_controls()
  patients <- load_patients()
  check_response_labels(controls, CONTROL_ITEM_REGEX, "Controls")
  check_response_labels(patients, PATIENT_ITEM_REGEX, "Patients")

  keep <- c("external_code", "gender", "age", "anni_scolarita", "Diagnosi",
            "total_errors", "total_matriks", "incomplete_correlate", "repetition",
            "wrong_principle", "difference", "r_ic")

  patients_scored <- add_error_scores(patients, PATIENT_ITEM_REGEX)
  controls_scored <- add_error_scores(controls, CONTROL_ITEM_REGEX)

  list(
    patients = patients_scored,
    general  = bind_rows(select(patients_scored, all_of(keep)),
                         select(controls_scored, all_of(keep))) %>%
      mutate(Diagnosi = factor(Diagnosi, levels = c("Bipolar", "Control", "Schizophrenic")))
  )
}

# p-value formatted as in the paper: "< .001", otherwise without leading zero.
fmt_p <- function(p) {
  vapply(p, function(v) {
    if (v < .001) return("< .001")
    sub("^0", "", formatC(v, format = "f", digits = if (v < .01) 3 else 2))
  }, character(1))
}

# Optimal matching in triplets (one Control, one Bipolar, one Schizophrenia
# participant). `data` has one row per participant with columns subject_key,
# diagnosis, age and the matching score. The smallest group is the anchor
# (Schizophrenia, n = 28). For each of the other two groups, the anchors are
# paired with distinct participants by solving the linear assignment problem
# (Hungarian algorithm, clue::solve_LSAP), which minimises the total distance
# over all anchors instead of pairing them one at a time. The distance is the
# squared accuracy difference (plus the age tie-break); pairs outside the
# caliper get a prohibitive distance, so the solution keeps as many pairs as
# possible within the caliper. Triplets whose accuracy range still exceeds the
# caliper are dropped. The squared distance makes the solution unique up to
# exact ties in accuracy and age; the procedure is deterministic.
match_three_groups <- function(data, score_col = "total_accuracy",
                               caliper = MATCH_CALIPER, age_weight = MATCH_AGE_WEIGHT) {
  groups <- c("Control", "Bipolar", "Schizophrenia")
  dat <- data %>%
    filter(diagnosis %in% groups) %>%
    mutate(diagnosis = as.character(diagnosis))
  stopifnot(all(groups %in% dat$diagnosis), !anyNA(dat[[score_col]]), !anyNA(dat$age))

  anchor_group <- dat %>% count(diagnosis) %>% arrange(n, diagnosis) %>% slice(1) %>% pull(diagnosis)
  anchors <- dat %>%
    filter(diagnosis == anchor_group) %>%
    arrange(.data[[score_col]], subject_key) %>%
    mutate(match_set = row_number())

  picks <- map(setdiff(groups, anchor_group), function(g) {
    pool <- dat %>% filter(diagnosis == g) %>% arrange(subject_key)
    diff <- abs(outer(anchors[[score_col]], pool[[score_col]], "-"))
    cost <- diff^2 + age_weight * abs(outer(anchors$age, pool$age, "-"))
    cost[diff > caliper + 1e-9] <- 1e3
    assignment <- as.integer(clue::solve_LSAP(cost))   # column chosen for each anchor
    pool[assignment, ] %>% mutate(match_set = anchors$match_set)
  })

  out <- bind_rows(anchors, picks) %>%
    group_by(match_set) %>%
    mutate(triplet_range = max(.data[[score_col]]) - min(.data[[score_col]])) %>%
    ungroup()
  dropped <- out %>% filter(triplet_range > caliper + 1e-9)
  if (nrow(dropped) > 0) {
    message(sprintf("Matching: %d triplet(s) outside the caliper (%.2f) dropped: %s",
                    n_distinct(dropped$match_set), caliper,
                    paste(dropped$subject_key, collapse = ", ")))
  }
  out <- out %>%
    filter(triplet_range <= caliper + 1e-9) %>%
    mutate(match_set = dense_rank(match_set),
           diagnosis = factor(diagnosis, levels = groups)) %>%
    arrange(match_set, diagnosis)
  if (nrow(out) == 0) stop("Matching produced no complete triplet.", call. = FALSE)
  out
}

write_session_info <- function(path = file.path(OUTPUT_DIR, "session_info.txt")) {
  lines <- capture.output(sessionInfo())
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    v <- tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) "not installed")
    lines <- c(lines, "", paste("CmdStan version:", v))
  }
  writeLines(lines, path)
  invisible(path)
}

MATRIKS_SETUP_DONE <- TRUE
