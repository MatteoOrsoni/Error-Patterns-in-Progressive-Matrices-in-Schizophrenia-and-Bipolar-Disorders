###############################################################################
# tools/simulate_example_data.R
# Writes SIMULATED data files with the same structure as the study files
# (columns, types and response labels documented in data/README.md), so that
# the pipeline can be run without access to the study data. The values are
# random: results obtained on these files say nothing about the study.
# Run from the repository root:
#   source("tools/simulate_example_data.R")
# Output: data_example/data_adult.xlsx, data_example/data_psy.xlsx
# Then run the pipeline on them, e.g.
#   Sys.setenv(MATRIKS_DATA_DIR = "data_example", MATRIKS_OUTPUT_DIR = "outputs_example",
#              MATRIKS_QUICK_TEST = "1")
#   source("run_all.R")
###############################################################################

suppressPackageStartupMessages({ library(writexl); library(dplyr) })

OUT_DIR <- "data_example"
dir.create(OUT_DIR, showWarnings = FALSE)
set.seed(2024)

CONTROL_ITEMS <- c(4, 13, 14, 15, 16, 27, 30, 33, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47)
PATIENT_ITEMS <- 34:51
RIC_ITEMS     <- c(15, 38, 40)   # MatriKS 12+ items whose repetition distractors are recorded as "r.ic"

# Distractor labels of each item (the same 18 stimuli in the two versions).
distractors <- read.csv(file.path("data", "distractor_distance.csv"), stringsAsFactors = FALSE) %>%
  mutate(distractor = ifelse(category == "R" & matriks12_item %in% RIC_ITEMS, "r.ic", distractor))
item_labels <- lapply(CONTROL_ITEMS, function(i) filter(distractors, matriks12_item == i))

n_items    <- length(CONTROL_ITEMS)
difficulty <- sort(rnorm(n_items, 0, 1))      # later items are harder
TYPE_PROBS <- list(                           # error-type probabilities (R, D, IC, WP)
  Control       = c(R = .30, D = .05, IC = .30, WP = .35),
  Bipolar       = c(R = .28, D = .08, IC = .30, WP = .34),
  Schizophrenic = c(R = .25, D = .10, IC = .30, WP = .35)
)

# Responses (labels) and response times (integer seconds) of one group.
simulate_items <- function(n, group, ability_mean, p_skip, p_fast) {
  ability <- rnorm(n, ability_mean, 1)
  fast    <- runif(n) < p_fast                # participants who respond very quickly
  resp <- matrix(NA_character_, n, n_items)
  time <- matrix(NA_real_, n, n_items)
  for (k in seq_len(n_items)) {
    labs    <- item_labels[[k]]
    correct <- runif(n) < plogis(ability - difficulty[k])
    type    <- sample(names(TYPE_PROBS[[group]]), n, replace = TRUE, prob = TYPE_PROBS[[group]])
    error   <- vapply(type, function(t) {
      pool <- unique(labs$distractor[labs$category == t])
      if (length(pool) == 0) pool <- unique(labs$distractor)
      pool[sample.int(length(pool), 1)]
    }, character(1))
    resp[, k] <- ifelse(correct, "correct", error)
    resp[runif(n) < p_skip, k] <- "skip"
    time[, k] <- ifelse(fast, rpois(n, 2), round(rlnorm(n, 3 + 0.05 * k, 0.6)))
  }
  list(resp = resp, time = time)
}

item_columns <- function(sim, version, items) {
  resp <- as.data.frame(sim$resp, stringsAsFactors = FALSE)
  time <- as.data.frame(sim$time)
  names(resp) <- paste0("response_mat", version, "_item.", items)
  names(time) <- paste0("time_mat", version, "_item.", items)
  bind_cols(resp, time)
}

years_text <- function(n, mean) {
  y <- sprintf("%.1f", pmin(pmax(round(rnorm(n, mean, 3)), 5), 21))
  y[runif(n) < .02] <- "NULL"                 # missing values are stored as "NULL"
  y
}

# ---- Healthy participants (MatriKS 12+) -------------------------------------------
n_adult <- 320
adult <- tibble(
  new_id         = sprintf("SIM-A%03d", seq_len(n_adult)),
  age            = round(runif(n_adult, 14, 75), 2),            # rows under 18 are excluded by the code
  gender         = sample(c("F", "M", "B"), n_adult, replace = TRUE, prob = c(.55, .44, .01)),
  diagnostic     = ifelse(runif(n_adult) < .05, "Example diagnosis", NA_character_),  # excluded
  anni_scolarita = years_text(n_adult, 13)
) %>%
  bind_cols(item_columns(simulate_items(n_adult, "Control", 0.3, .005, 0), 50, CONTROL_ITEMS))

# ---- Psychiatric inpatients (MatriKS88) ---------------------------------------------
n_bd <- 34; n_sz <- 28
n_psy <- n_bd + n_sz
dx <- c(rep("Bipolar", n_bd), rep("Schizophrenic", n_sz))
sim_bd <- simulate_items(n_bd, "Bipolar",       -0.6, .02, .10)
sim_sz <- simulate_items(n_sz, "Schizophrenic", -0.9, .02, .15)
psy_items <- item_columns(list(resp = rbind(sim_bd$resp, sim_sz$resp),
                               time = rbind(sim_bd$time, sim_sz$time)), 88, PATIENT_ITEMS)
bprs_pre  <- round(rnorm(n_psy, 50, 10))
honos_pre <- round(rnorm(n_psy, 35, 8))
psy <- tibble(
  external_code     = sprintf("SIM-P%02d", seq_len(n_psy)),
  Diagnosi          = dx,
  age               = round(runif(n_psy, 19, 75), 2),
  gender            = sample(c("F", "M"), n_psy, replace = TRUE),
  anni_scolarita    = years_text(n_psy, 11),
  BPRS_PRE          = bprs_pre,
  BPRS_POST         = bprs_pre - round(abs(rnorm(n_psy, 10, 6))),
  HONOS_PRE         = honos_pre,
  HONOS_POST        = honos_pre - round(abs(rnorm(n_psy, 6, 4))),
  Matrici_Attentive = sub(".", ",", sprintf("%.1f", rnorm(n_psy, 38, 9)), fixed = TRUE),  # comma decimals
  SPAN_diretto      = round(rnorm(n_psy, 5.5, 0.9), 2),
  SPAN_inverso      = round(rnorm(n_psy, 3.6, 0.8), 2)
) %>%
  bind_cols(psy_items)

# Missing cognitive scores (Model 2 compares included and excluded patients).
no_cognitive <- sample(n_psy, 8)
psy$Matrici_Attentive[no_cognitive] <- NA
psy$SPAN_diretto[no_cognitive]      <- NA
psy$SPAN_inverso[no_cognitive]      <- NA
psy$BPRS_POST[sample(n_psy, 3)]     <- NA
# Two patients with an incomplete administration (items not presented are empty).
for (i in sample(n_psy, 2)) {
  last_items <- tail(PATIENT_ITEMS, 4)
  psy[i, paste0("response_mat88_item.", last_items)] <- NA
  psy[i, paste0("time_mat88_item.", last_items)]     <- NA
}

write_xlsx(adult, file.path(OUT_DIR, "data_adult.xlsx"))
write_xlsx(psy,   file.path(OUT_DIR, "data_psy.xlsx"))
message("Simulated data written to ", OUT_DIR, "/ (", n_adult, " adults, ", n_psy, " patients).")
