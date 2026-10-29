###############################################################################
# 03_matching_and_models.R
# Bayesian multinomial mixed-effects models of error type (Stan via cmdstanr).
#   1. Item-level responses on the 18 common items
#   2. Performance matching: optimal matching in triplets (one Control, one
#      Bipolar, one Schizophrenia participant) on accuracy, with caliper
#      (match_three_groups() in 00_setup.R)  -> Table 4 (balance of the matched sample)
#   3. Model 1: diagnosis (matched three-group sample)
#   4. Model 2: attention and verbal memory x diagnosis (patients only);
#      patients included vs excluded for missing cognitive scores (Table E)
#   5. Posterior summaries, sampling diagnostics, diagnostic plots
# Outputs (in outputs/): data_and_matching_outputs.xlsx, matched_participants.csv (also read by
#   script 07), matched_sample_balance.xlsx (Table 4), model2_missing_covariates.xlsx,
#   multinomial_model_results.xlsx, model1_fit.rds, model2_fit.rds,
#   model_stan_data.rds (read by script 08), diagnostics/*.png
###############################################################################

if (!exists("MATRIKS_SETUP_DONE")) source(file.path("R", "00_setup.R"))
require_packages(BAYES_PACKAGES)
suppressPackageStartupMessages({ library(cmdstanr); library(posterior); library(bayesplot) })

OUTCOME_LEVELS      <- c("R", "D", "IC", "WP")          # R = reference category
NONREFERENCE_LEVELS <- c("D vs R", "IC vs R", "WP vs R")
DIAG_DIR <- file.path(OUTPUT_DIR, "diagnostics")
dir.create(DIAG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Helpers --------------------------------------------------------------------
harmonize_diagnosis <- function(x) {
  x_chr <- str_to_lower(trimws(as.character(x)))
  case_when(
    x_chr == "control"                                ~ "Control",
    x_chr == "bipolar"                                ~ "Bipolar",
    x_chr %in% c("schizophrenia", "schizophrenic")    ~ "Schizophrenia",
    TRUE                                              ~ as.character(x)
  )
}

# Classify a MatriKS response label into correct / R / D / IC / WP.
# Skipped items become NA (excluded). "r.ic" follows RIC_POLICY (see 00_setup.R).
classify_response <- function(x, ric_policy = c("R", "IC", "drop")) {
  ric_policy <- match.arg(ric_policy)
  x <- str_to_lower(trimws(as.character(x)))
  case_when(
    is.na(x) | x %in% c("", NON_RESPONSE_LABELS, "null", "na") ~ NA_character_,
    x == "r.ic" & ric_policy == "drop"                         ~ NA_character_,
    x == "r.ic" & ric_policy == "IC"                           ~ "IC",
    x == "r.ic" & ric_policy == "R"                            ~ "R",
    startsWith(x, "correct")                                   ~ "correct",
    startsWith(x, "wp")                                        ~ "WP",
    startsWith(x, "ic")                                        ~ "IC",
    startsWith(x, "r")                                         ~ "R",
    startsWith(x, "d")                                         ~ "D",
    TRUE                                                       ~ "UNKNOWN"
  )
}

# z-scores; a constant variable gets 0.
safe_z <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

# Long format: one row per participant x item. `item` is the harmonised
# index 1..18 (position in the item list), identical for the two versions.
extract_item_responses <- function(data, id_col, dx_col, item_stub, item_ids, dataset_label) {
  item_cols <- paste0(item_stub, item_ids)
  check_columns(data, item_cols, dataset_label)
  out <- data %>%
    transmute(id        = as.character(.data[[id_col]]),
              diagnosis = harmonize_diagnosis(.data[[dx_col]]),
              across(all_of(item_cols), as.character)) %>%
    pivot_longer(all_of(item_cols), names_to = "item_col", values_to = "response_raw") %>%
    mutate(item_raw       = as.integer(str_extract(item_col, "\\d+$")),
           item           = match(item_raw, item_ids),
           source_dataset = dataset_label,
           response_type  = classify_response(response_raw, RIC_POLICY))
  bad <- which(out$response_type == "UNKNOWN")
  if (length(bad) > 0) {
    warning(sprintf("%s: %d responses with unrecognised labels (discarded): %s", dataset_label,
                    length(bad), paste(sort(unique(out$response_raw[bad])), collapse = ", ")),
            call. = FALSE)
  }
  filter(out, !is.na(response_type), response_type != "UNKNOWN")
}

# Participant-level totals. total_accuracy = correct / answered items
# (skipped items excluded from the denominator); used for matching.
summarise_participant_items <- function(item_level_data) {
  item_level_data %>%
    group_by(id, diagnosis, source_dataset) %>%
    summarise(total_items   = n(),
              total_correct = sum(response_type == "correct"),
              total_errors  = sum(response_type != "correct"),
              R  = sum(response_type == "R"),  D  = sum(response_type == "D"),
              IC = sum(response_type == "IC"), WP = sum(response_type == "WP"),
              .groups = "drop") %>%
    mutate(total_accuracy = total_correct / total_items)
}

# Item administration on the 18 common items, from the raw labels: items not
# administered (empty cell or "NULL") and items skipped ("skip"), per participant.
count_item_administration <- function(data, item_stub, item_ids, dataset_label) {
  data %>%
    transmute(subject_key = paste(dataset_label, as.character(external_code), sep = "::"),
              across(all_of(paste0(item_stub, item_ids)), as.character)) %>%
    pivot_longer(-subject_key, values_to = "response_raw") %>%
    group_by(subject_key) %>%
    mutate(response_raw = str_to_lower(trimws(response_raw))) %>%
    summarise(n_not_administered = sum(is.na(response_raw) | response_raw %in% c("", "null", "na")),
              n_skipped          = sum(response_raw %in% NON_RESPONSE_LABELS, na.rm = TRUE),
              .groups = "drop")
}

# Standardised mean difference (pooled SD); for a proportion p, SD = sqrt(p(1-p)).
smd <- function(x, y) {
  s <- sqrt((var(x, na.rm = TRUE) + var(y, na.rm = TRUE)) / 2)
  if (is.na(s) || s == 0) return(0)
  (mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)) / s
}

compute_vif_table <- function(data, predictors) {
  x <- data %>% select(all_of(predictors)) %>% mutate(across(everything(), as.numeric))
  map_dfr(predictors, function(pred) {
    fit <- lm(as.formula(paste(pred, "~", paste(setdiff(predictors, pred), collapse = " + "))), data = x)
    r2  <- summary(fit)$r.squared
    tibble(term = pred, vif = if (is.na(r2) || r2 >= 1) Inf else 1 / (1 - r2))
  })
}

# ---- 1. Item-level responses ------------------------------------------------------
controls <- load_controls()
patients <- load_patients()
check_columns(patients, c("Matrici_Attentive", "SPAN_diretto", "SPAN_inverso"), PSY_FILE)

control_items <- extract_item_responses(controls, "external_code", "Diagnosi",
                                        CONTROL_ITEM_STUB, CONTROL_ITEMS, "adult") %>%
  mutate(subject_key = paste(source_dataset, id, sep = "::"))
patient_items <- extract_item_responses(patients, "external_code", "Diagnosi",
                                        PATIENT_ITEM_STUB, PATIENT_ITEMS, "psy") %>%
  mutate(subject_key = paste(source_dataset, id, sep = "::"))

item_level_all <- bind_rows(control_items, patient_items) %>%
  filter(diagnosis %in% c("Control", "Bipolar", "Schizophrenia"))

cat("\nItems with 'r.ic' responses (MatriKS88 / MatriKS 12+ item numbers):\n")
print(item_level_all %>% filter(str_to_lower(response_raw) == "r.ic") %>%
        count(item, name = "n_r.ic_responses") %>%
        mutate(matriks88_item = PATIENT_ITEMS[item], matriks12_item = CONTROL_ITEMS[item]))

participant_summary_all <- summarise_participant_items(item_level_all) %>%
  mutate(subject_key = paste(source_dataset, id, sep = "::"))
print(count(participant_summary_all, diagnosis))

# Demographics and item administration, used for the age tie-break of the
# matching and for the balance table.
participant_info <- bind_rows(
  transmute(controls, subject_key = paste("adult", external_code, sep = "::"),
            age, gender, education = anni_scolarita),
  transmute(patients, subject_key = paste("psy", external_code, sep = "::"),
            age, gender, education = anni_scolarita)
) %>%
  left_join(bind_rows(
    count_item_administration(controls, CONTROL_ITEM_STUB, CONTROL_ITEMS, "adult"),
    count_item_administration(patients, PATIENT_ITEM_STUB, PATIENT_ITEMS, "psy")
  ), by = "subject_key")

participant_summary_all <- participant_summary_all %>%
  left_join(participant_info, by = "subject_key") %>%
  mutate(accuracy_of_18 = total_correct / length(PATIENT_ITEMS))

# ---- 2. Performance matching (Table 4) ----------------------------------------------
matched_participants <- match_three_groups(participant_summary_all)
print(count(matched_participants, diagnosis))
cat("\nBipolar patients not matched:\n")
print(participant_summary_all %>%
        filter(diagnosis == "Bipolar", !subject_key %in% matched_participants$subject_key) %>%
        select(subject_key, total_correct, total_items, total_accuracy, age, gender, education))

# Matched sample, read by script 07 (one row per participant).
write.csv(matched_participants %>%
            transmute(subject_key, id, diagnosis = as.character(diagnosis), match_set,
                      total_correct, total_items, total_accuracy),
          file.path(OUTPUT_DIR, "matched_participants.csv"), row.names = FALSE)

matching_quality <- matched_participants %>%
  group_by(match_set) %>%
  summarise(range_accuracy = max(total_accuracy) - min(total_accuracy), .groups = "drop")
cat("\nAccuracy range within matched triplets:\n"); print(summary(matching_quality$range_accuracy))
cat("\nANOVA of accuracy by group in the matched sample:\n")
print(summary(aov(total_accuracy ~ diagnosis, data = matched_participants)))

matched_item_data <- item_level_all %>%
  semi_join(select(matched_participants, subject_key), by = "subject_key") %>%
  left_join(select(matched_participants, subject_key, match_set), by = "subject_key")

matching_overview <- bind_rows(
  mutate(participant_summary_all, sample = "full",    diagnosis = as.character(diagnosis)),
  mutate(matched_participants,    sample = "matched", diagnosis = as.character(diagnosis))
) %>%
  group_by(sample, diagnosis) %>%
  summarise(n_participants      = n(),
            mean_total_accuracy = mean(total_accuracy), sd_total_accuracy = sd(total_accuracy),
            mean_total_errors   = mean(total_errors),   sd_total_errors   = sd(total_errors),
            .groups = "drop")
print(matching_overview)

# Table 4 (balance): matched groups and the full control sample. Matched controls
# are also located in the accuracy distribution of all controls (percentile).
all_controls <- filter(participant_summary_all, diagnosis == "Control")
balance_rows <- function(d) {
  tibble(
    n                    = nrow(d),
    accuracy_mean        = mean(d$total_accuracy),  accuracy_sd = sd(d$total_accuracy),
    accuracy_min         = min(d$total_accuracy),   accuracy_max = max(d$total_accuracy),
    accuracy_of_18_mean  = mean(d$accuracy_of_18),  accuracy_of_18_sd = sd(d$accuracy_of_18),
    correct_le4_n        = sum(d$total_correct <= 4),   # not above chance (1/8) on 18 items
    items_answered_mean  = mean(d$total_items),
    skipped_total        = sum(d$n_skipped),
    not_administered_total = sum(d$n_not_administered),
    age_mean             = mean(d$age),  age_sd = sd(d$age),
    age_min              = min(d$age),   age_max = max(d$age),
    female_n             = sum(d$gender == "F"),  female_pct = 100 * mean(d$gender == "F"),
    education_mean       = mean(d$education, na.rm = TRUE),
    education_sd         = sd(d$education, na.rm = TRUE),
    education_n          = sum(!is.na(d$education)),
    accuracy_percentile_in_controls_median =
      median(vapply(d$total_accuracy, function(a) 100 * mean(all_controls$total_accuracy <= a), numeric(1)))
  )
}
balance_table <- bind_rows(
  matched_participants %>% mutate(group = paste(diagnosis, "(matched)")) %>%
    group_by(group) %>% group_modify(~ balance_rows(.x)) %>% ungroup(),
  balance_rows(all_controls) %>% mutate(group = "Control (all)")
)
balance_smd <- {
  g <- split(matched_participants, as.character(matched_participants$diagnosis))
  pairs <- list(c("Bipolar", "Schizophrenia"), c("Control", "Schizophrenia"), c("Control", "Bipolar"))
  map_dfr(pairs, function(pr) {
    a <- g[[pr[1]]]; b <- g[[pr[2]]]
    tibble(comparison = paste(pr[1], "vs", pr[2]),
           accuracy = smd(a$total_accuracy, b$total_accuracy),
           accuracy_of_18 = smd(a$accuracy_of_18, b$accuracy_of_18),
           age = smd(a$age, b$age),
           female = smd(a$gender == "F", b$gender == "F"),
           education = smd(a$education, b$education))
  })
}
cat("\nBalance of the matched sample:\n"); print(balance_table, width = Inf)
cat("\nStandardised mean differences in the matched sample:\n"); print(balance_smd, width = Inf)
write_xlsx(list(balance = balance_table, smd = balance_smd,
                matched_participants = mutate(matched_participants, diagnosis = as.character(diagnosis))),
           file.path(OUTPUT_DIR, "matched_sample_balance.xlsx"))

# ---- 3-4. Model data --------------------------------------------------------------
# Model 1: dummy-coded diagnosis, controls as reference.
model1_subject_data <- matched_participants %>%
  mutate(dx_bipolar       = as.numeric(diagnosis == "Bipolar"),
         dx_schizophrenic = as.numeric(diagnosis == "Schizophrenia")) %>%
  arrange(diagnosis, id)

model1_error_data <- matched_item_data %>%
  filter(response_type %in% OUTCOME_LEVELS) %>%
  left_join(select(model1_subject_data, subject_key, dx_bipolar, dx_schizophrenic), by = "subject_key")

# Model 2: patients only. Cognitive scores are z-standardised across all
# patients with a valid score on that measure (before listwise deletion);
# diagnosis coded 0 = Bipolar, 1 = Schizophrenia; interactions = product terms.
clinical_covariates <- patients %>%
  transmute(
    id                = external_code,
    diagnosis         = harmonize_diagnosis(Diagnosi),
    subject_key       = paste("psy", id, sep = "::"),
    matrici_attentive = to_num(Matrici_Attentive),
    span_forward      = to_num(SPAN_diretto),
    span_backward     = to_num(SPAN_inverso)
  ) %>%
  filter(diagnosis %in% c("Bipolar", "Schizophrenia")) %>%
  mutate(dx_schizophrenic    = as.numeric(diagnosis == "Schizophrenia"),
         matrici_attentive_z = safe_z(matrici_attentive),
         span_forward_z      = safe_z(span_forward),
         span_backward_z     = safe_z(span_backward),
         dx_x_attentive      = dx_schizophrenic * matrici_attentive_z,
         dx_x_span_forward   = dx_schizophrenic * span_forward_z,
         dx_x_span_backward  = dx_schizophrenic * span_backward_z)

cat("\nMissing cognitive scores (these patients are excluded from Model 2):\n")
print(colSums(is.na(select(clinical_covariates, matrici_attentive, span_forward, span_backward))))

MODEL1_PREDICTORS <- c("dx_bipolar", "dx_schizophrenic")
MODEL2_PREDICTORS <- c("dx_schizophrenic", "matrici_attentive_z", "span_forward_z", "span_backward_z",
                       "dx_x_attentive", "dx_x_span_forward", "dx_x_span_backward")

model2_subject_data <- clinical_covariates %>%
  inner_join(select(participant_summary_all, subject_key, total_items, total_correct,
                    total_errors, total_accuracy), by = "subject_key") %>%
  arrange(diagnosis, id)

model2_error_data <- patient_items %>%
  filter(diagnosis %in% c("Bipolar", "Schizophrenia"), response_type %in% OUTCOME_LEVELS) %>%
  inner_join(select(model2_subject_data, subject_key, all_of(MODEL2_PREDICTORS)), by = "subject_key")

# Patients included in Model 2 (all three cognitive scores available; the others
# are dropped by prepare_model_dataset() below) vs patients excluded for missing
# scores: missingness pattern, and comparison on demographic, clinical and MatriKS
# variables (Supplementary Table; SMD = excluded - included).
model2_inclusion <- clinical_covariates %>%
  transmute(subject_key, diagnosis,
            missing_attention = is.na(matrici_attentive),
            missing_span      = is.na(span_forward) | is.na(span_backward),
            included_model2   = !missing_attention & !missing_span) %>%
  left_join(participant_summary_all %>%
              select(subject_key, age, gender, education, total_correct, total_errors,
                     total_accuracy, all_of(OUTCOME_LEVELS)), by = "subject_key") %>%
  left_join(transmute(patients, subject_key = paste("psy", external_code, sep = "::"),
                      BPRS_PRE, HONOS_PRE), by = "subject_key") %>%
  mutate(female = as.numeric(gender == "F"), correct_le4 = as.numeric(total_correct <= 4))

inclusion_vars <- c("age", "female", "education", "BPRS_PRE", "HONOS_PRE",
                    "total_accuracy", "total_errors", "correct_le4")
inclusion_rows <- function(d) {
  inc_d <- filter(d, included_model2); exc_d <- filter(d, !included_model2)
  bind_rows(
    map_dfr(inclusion_vars, function(v) {
      inc <- inc_d[[v]]; exc <- exc_d[[v]]
      tibble(variable = v,
             included_n = sum(!is.na(inc)), included_mean = mean(inc, na.rm = TRUE),
             included_sd = sd(inc, na.rm = TRUE),
             excluded_n = sum(!is.na(exc)), excluded_mean = mean(exc, na.rm = TRUE),
             excluded_sd = sd(exc, na.rm = TRUE),
             smd_excluded_minus_included = smd(exc, inc))
    }),
    # Error-type composition: share of each type among all errors of the subgroup
    # (n = number of errors)
    map_dfr(OUTCOME_LEVELS, function(o) {
      share <- function(x) sum(x[[o]]) / sum(x[, OUTCOME_LEVELS])
      tibble(variable = paste0("share_of_errors_", o),
             included_n = sum(inc_d[, OUTCOME_LEVELS]), included_mean = share(inc_d),
             excluded_n = sum(exc_d[, OUTCOME_LEVELS]), excluded_mean = share(exc_d))
    })
  )
}
model2_inclusion_table <- bind_rows(
  inclusion_rows(model2_inclusion) %>% mutate(group = "All patients"),
  model2_inclusion %>% group_by(group = diagnosis) %>% group_modify(~ inclusion_rows(.x)) %>% ungroup()
) %>% select(group, everything())
model2_missing_pattern <- count(model2_inclusion, diagnosis, included_model2, missing_attention, missing_span)
model2_inclusion_test  <- fisher.test(table(model2_inclusion$diagnosis, model2_inclusion$included_model2))
cat("\nModel 2: patients included / excluded for missing cognitive scores:\n")
print(model2_missing_pattern)
cat(sprintf("Fisher's exact test, diagnosis x inclusion: OR = %.2f, p = %.3f\n",
            model2_inclusion_test$estimate, model2_inclusion_test$p.value))
print(model2_inclusion_table, n = Inf, width = Inf)
write_xlsx(list(inclusion_comparison = model2_inclusion_table,
                missing_pattern      = model2_missing_pattern,
                fisher_test          = tibble(odds_ratio = unname(model2_inclusion_test$estimate),
                                              p_value    = model2_inclusion_test$p.value),
                patients             = select(model2_inclusion, -gender)),
           file.path(OUTPUT_DIR, "model2_missing_covariates.xlsx"))

# Optional sensitivity analysis (MATRIKS_EXCLUDE_ITEMS, see 00_setup.R): the
# excluded item(s) are removed from both models, for all groups. Participants and
# matching are unchanged. The number of items passed to Stan stays 18: an item
# without responses only adds an unconstrained random intercept, which does not
# change the posterior of the other parameters.
if (length(MODEL_EXCLUDED_ITEMS) > 0) {
  excluded_idx <- match(MODEL_EXCLUDED_ITEMS, PATIENT_ITEMS)
  message("SENSITIVITY ANALYSIS: excluding MatriKS88 item(s) ", paste(MODEL_EXCLUDED_ITEMS, collapse = ", "),
          " (MatriKS 12+ item(s) ", paste(CONTROL_ITEMS[excluded_idx], collapse = ", "), ") from Models 1-2.")
  model1_error_data <- filter(model1_error_data, !item %in% excluded_idx)
  model2_error_data <- filter(model2_error_data, !item %in% excluded_idx)
}

vif_table_model2 <- compute_vif_table(model2_subject_data, MODEL2_PREDICTORS)
cat("\nModel 2 VIFs (participant level):\n"); print(vif_table_model2)

write_xlsx(
  list(matching_overview    = matching_overview,
       matching_quality     = matching_quality,
       matched_participants = mutate(matched_participants, diagnosis = as.character(diagnosis)),
       model1_subject_data  = mutate(model1_subject_data, diagnosis = as.character(diagnosis)),
       model2_subject_data  = model2_subject_data,
       model2_vif           = vif_table_model2),
  file.path(OUTPUT_DIR, "data_and_matching_outputs.xlsx")
)

# ---- 5. Stan data, fitting, summaries -------------------------------------------------
prepare_model_dataset <- function(data, predictors, model_label) {
  prepared <- data %>%
    filter(!is.na(response_type), !is.na(subject_key), !is.na(item)) %>%
    filter(if_all(all_of(predictors), ~ !is.na(.x))) %>%
    mutate(response_type = factor(response_type, levels = OUTCOME_LEVELS))
  if (nrow(prepared) == 0) stop(model_label, ": no complete observations.", call. = FALSE)
  cat(sprintf("\n%s - errors by group and type:\n", model_label))
  print(table(prepared$diagnosis, prepared$response_type))
  prepared
}

make_stan_data <- function(error_data, predictors) {
  X <- as.matrix(mutate(select(error_data, all_of(predictors)), across(everything(), as.numeric)))
  list(N = nrow(error_data), K = length(OUTCOME_LEVELS), P = ncol(X),
       J = n_distinct(error_data$subject_key), I = length(PATIENT_ITEMS),
       subject = as.integer(factor(error_data$subject_key)),
       item    = as.integer(error_data$item),
       y       = as.integer(factor(error_data$response_type, levels = OUTCOME_LEVELS)),
       X       = X)
}

check_stan_csv_dir <- function(path) {
  if (grepl("[^ -~]", path)) {
    stop("The CmdStan output folder contains non-ASCII characters: ", path,
         "\nSet an ASCII-only folder, e.g. Sys.setenv(MATRIKS_STAN_CSV_DIR = 'C:/stan_csv').",
         call. = FALSE)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

fit_model <- function(model, stan_data, label) {
  fit <- model$sample(
    data = stan_data, seed = SEED,
    chains = MCMC$chains, parallel_chains = MCMC$parallel_chains,
    iter_warmup = MCMC$iter_warmup, iter_sampling = MCMC$iter_sampling,
    adapt_delta = MCMC$adapt_delta, max_treedepth = MCMC$max_treedepth,
    refresh = 500, output_dir = check_stan_csv_dir(STAN_CSV_DIR), output_basename = label
  )
  fit$save_object(file.path(OUTPUT_DIR, paste0(label, "_fit.rds")))   # saved before any post-processing
  fit
}

summarise_draw_set <- function(fit, variable_name) {
  draws <- as_draws_df(fit$draws(variables = variable_name)) %>% as_tibble()
  draws %>%
    select(-any_of(c(".chain", ".iteration", ".draw"))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    group_by(variable) %>%
    summarise(mean = mean(value), median = median(value), sd = sd(value),
              q2.5 = quantile(value, 0.025), q97.5 = quantile(value, 0.975),
              prob_gt_zero = mean(value > 0), prob_lt_zero = mean(value < 0),
              .groups = "drop") %>%
    left_join(as_tibble(fit$summary(variables = variable_name)) %>%
                select(variable, rhat, ess_bulk, ess_tail), by = "variable")
}

index_of <- function(variable, pattern, pos) as.integer(str_match(variable, pattern)[, pos])

extract_model_results <- function(fit, predictor_names) {
  beta_pat <- "beta\\[(\\d+),(\\d+)\\]"
  list(
    beta = summarise_draw_set(fit, "beta") %>%
      mutate(term       = predictor_names[index_of(variable, beta_pat, 2)],
             comparison = NONREFERENCE_LEVELS[index_of(variable, beta_pat, 3)],
             odds_ratio_median = exp(median), odds_ratio_q2.5 = exp(q2.5),
             odds_ratio_q97.5  = exp(q97.5),  cri95_excludes_0 = q2.5 > 0 | q97.5 < 0) %>%
      select(term, comparison, everything()),
    alpha = summarise_draw_set(fit, "alpha") %>%
      mutate(comparison = NONREFERENCE_LEVELS[index_of(variable, "alpha\\[(\\d+)\\]", 2)]) %>%
      select(comparison, everything()),
    mu_beta = summarise_draw_set(fit, "mu_beta") %>%
      mutate(term = predictor_names[index_of(variable, "mu_beta\\[(\\d+)\\]", 2)]) %>%
      select(term, everything()),
    tau_beta = summarise_draw_set(fit, "tau_beta") %>%
      mutate(term = predictor_names[index_of(variable, "tau_beta\\[(\\d+)\\]", 2)]) %>%
      select(term, everything()),
    sigma_u = summarise_draw_set(fit, "sigma_u") %>%        # participant random-intercept SDs
      mutate(comparison = NONREFERENCE_LEVELS[index_of(variable, "sigma_u\\[(\\d+)\\]", 2)]) %>%
      select(comparison, everything()),
    sigma_v = summarise_draw_set(fit, "sigma_v") %>%        # item random-intercept SDs
      mutate(comparison = NONREFERENCE_LEVELS[index_of(variable, "sigma_v\\[(\\d+)\\]", 2)]) %>%
      select(comparison, everything())
  )
}

compute_model_diagnostics <- function(fit, model_label, stan_data) {
  summ <- fit$summary(variables = c("alpha", "mu_beta", "tau_beta", "beta", "sigma_u", "sigma_v"))
  dsum <- fit$diagnostic_summary(quiet = TRUE)
  tibble(model = model_label, observations = stan_data$N, participants = stan_data$J,
         items = stan_data$I, predictors = stan_data$P,
         max_rhat = max(summ$rhat, na.rm = TRUE),
         min_bulk_ess = min(summ$ess_bulk, na.rm = TRUE),
         min_tail_ess = min(summ$ess_tail, na.rm = TRUE),
         divergences = sum(dsum$num_divergent),
         treedepth_hits = sum(dsum$num_max_treedepth),
         min_ebfmi = min(dsum$ebfmi, na.rm = TRUE))
}

save_diagnostic_plots <- function(fit, predictor_names, model_label) {
  P <- length(predictor_names); K1 <- length(NONREFERENCE_LEVELS)
  trace_vars <- c(paste0("alpha[", seq_len(K1), "]"), paste0("beta[1,", seq_len(K1), "]"),
                  paste0("tau_beta[", seq_len(P), "]"), paste0("sigma_u[", seq_len(K1), "]"),
                  paste0("sigma_v[", seq_len(K1), "]"))
  trace_plot <- mcmc_trace(fit$draws(variables = trace_vars), facet_args = list(ncol = 4)) +
    ggtitle(sprintf("%s trace plots", model_label))

  beta_draws <- fit$draws(variables = "beta")
  idx <- str_match(dimnames(beta_draws)$variable, "beta\\[(\\d+),(\\d+)\\]")
  dimnames(beta_draws)$variable <- paste0(predictor_names[as.integer(idx[, 2])], " | ",
                                          NONREFERENCE_LEVELS[as.integer(idx[, 3])])
  interval_plot <- mcmc_intervals(beta_draws, prob = 0.8, prob_outer = 0.95) +
    ggtitle(sprintf("%s coefficients (log-odds vs R)", model_label))

  ggsave(file.path(DIAG_DIR, paste0(model_label, "_trace.png")), trace_plot,
         width = 12, height = 9, dpi = 300)
  ggsave(file.path(DIAG_DIR, paste0(model_label, "_intervals.png")), interval_plot,
         width = 11, height = max(5, 0.35 * P * K1), dpi = 300)
}

multinomial_model <- cmdstan_model(STAN_FILE)

model1_stan_data <- make_stan_data(prepare_model_dataset(model1_error_data, MODEL1_PREDICTORS, "Model 1"),
                                   MODEL1_PREDICTORS)
model2_stan_data <- make_stan_data(prepare_model_dataset(model2_error_data, MODEL2_PREDICTORS, "Model 2"),
                                   MODEL2_PREDICTORS)
saveRDS(list(model1 = model1_stan_data, model2 = model2_stan_data),     # prior sensitivity (script 08)
        file.path(OUTPUT_DIR, "model_stan_data.rds"))

model1_fit <- fit_model(multinomial_model, model1_stan_data, "model1")
model2_fit <- fit_model(multinomial_model, model2_stan_data, "model2")

model1_results <- extract_model_results(model1_fit, MODEL1_PREDICTORS)
model2_results <- extract_model_results(model2_fit, MODEL2_PREDICTORS)

diagnostics <- bind_rows(compute_model_diagnostics(model1_fit, "model1", model1_stan_data),
                         compute_model_diagnostics(model2_fit, "model2", model2_stan_data))
print(diagnostics, width = Inf)

save_diagnostic_plots(model1_fit, MODEL1_PREDICTORS, "model1")
save_diagnostic_plots(model2_fit, MODEL2_PREDICTORS, "model2")

write_xlsx(
  list(model1_betas      = model1_results$beta,     model1_intercepts = model1_results$alpha,
       model1_mu_beta    = model1_results$mu_beta,  model1_tau_beta   = model1_results$tau_beta,
       model1_random_sd  = model1_results$sigma_u,  model1_item_sd    = model1_results$sigma_v,
       model1_diagnostics = filter(diagnostics, model == "model1"),
       model2_betas      = model2_results$beta,     model2_intercepts = model2_results$alpha,
       model2_mu_beta    = model2_results$mu_beta,  model2_tau_beta   = model2_results$tau_beta,
       model2_random_sd  = model2_results$sigma_u,  model2_item_sd    = model2_results$sigma_v,
       model2_diagnostics = filter(diagnostics, model == "model2")),
  file.path(OUTPUT_DIR, "multinomial_model_results.xlsx")
)
