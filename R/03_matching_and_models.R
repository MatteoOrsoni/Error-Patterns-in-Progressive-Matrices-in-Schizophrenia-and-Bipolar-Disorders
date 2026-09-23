###############################################################################
# 03_matching_and_models.R
# Bayesian multinomial mixed-effects models of error type (Stan via cmdstanr).
#   1. Item-level responses on the 18 common items
#   2. Performance matching: triplets (one Control, one Bipolar, one
#      Schizophrenia participant) on total accuracy       -> Figure 1
#   3. Model 1: diagnosis (matched three-group sample)
#   4. Model 2: attention and verbal memory x diagnosis (patients only)
#   5. Posterior summaries, sampling diagnostics, diagnostic plots
# Outputs (in outputs/): Figure1_matched_accuracy_density.png,
#   data_and_matching_outputs.xlsx, multinomial_model_results.xlsx,
#   model1_fit.rds, model2_fit.rds, diagnostics/*.png
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

# Greedy nearest-neighbour matching in triplets, without replacement and
# without caliper. The smallest group is the anchor; its participants are
# processed in ascending order of accuracy (ties broken by ID) and each is
# paired with the closest still-available participant of each other group.
# The procedure is deterministic.
match_three_groups <- function(data, score_col = "total_accuracy") {
  groups <- c("Control", "Bipolar", "Schizophrenia")
  dat <- data %>%
    filter(diagnosis %in% groups) %>%
    mutate(diagnosis = factor(diagnosis, levels = groups))
  stopifnot(all(groups %in% as.character(dat$diagnosis)))

  anchor_group <- dat %>% count(diagnosis) %>% arrange(n, diagnosis) %>%
    slice(1) %>% pull(diagnosis) %>% as.character()
  available    <- split(dat, dat$diagnosis)
  anchor_tbl   <- available[[anchor_group]] %>% arrange(.data[[score_col]], subject_key)
  other_groups <- setdiff(groups, anchor_group)
  used         <- setNames(rep(list(character(0)), length(other_groups)), other_groups)

  matched_sets <- list()
  for (i in seq_len(nrow(anchor_tbl))) {
    anchor_row <- anchor_tbl[i, ]
    picks <- list(anchor_row)
    for (g in other_groups) {
      pool <- available[[g]] %>%
        filter(!(subject_key %in% used[[g]])) %>%
        mutate(distance = abs(.data[[score_col]] - anchor_row[[score_col]][1])) %>%
        arrange(distance, .data[[score_col]], subject_key)
      if (nrow(pool) == 0) { picks <- NULL; break }
      chosen <- slice(pool, 1)
      used[[g]] <- c(used[[g]], chosen$subject_key)
      picks[[length(picks) + 1]] <- select(chosen, -distance)
    }
    if (!is.null(picks)) {
      matched_sets[[length(matched_sets) + 1]] <-
        bind_rows(picks) %>% mutate(match_set = length(matched_sets) + 1)
    }
  }
  out <- bind_rows(matched_sets)
  if (nrow(out) == 0) stop("Matching produced no complete triplet.", call. = FALSE)
  out
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

# ---- 2. Performance matching (Figure 1) ------------------------------------------
matched_participants <- match_three_groups(participant_summary_all)

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

figure1 <- ggplot(
  matched_participants %>%
    mutate(diagnosis = recode(as.character(diagnosis), "Schizophrenia" = "Schizophrenic")),
  aes(x = total_accuracy, color = diagnosis, fill = diagnosis)
) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = DX_COLORS) +
  scale_fill_manual(values  = DX_COLORS) +
  labs(title = "Distribution of Total MatriKS Accuracy (Matched Sample)",
       x = "Total Accuracy", y = "Density", color = "Conditions", fill = "Conditions") +
  DX_THEME +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggsave(file.path(OUTPUT_DIR, "Figure1_matched_accuracy_density.png"), figure1,
       width = 10, height = 6, dpi = 300)

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
