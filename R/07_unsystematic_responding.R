###############################################################################
# 07_unsystematic_responding.R
# Unsystematic responding (Sections 2.3.3.1-2.3.3.2; Appendices C-D):
#   1. Item-level responses and response times (RT) on the 18 common items
#   2. Matched sample: the participants matched by script 03 (optimal matching
#      in triplets with caliper, outputs/matched_participants.csv)
#   3. RT covariates: item-centred RT, split into a between-participant part
#      (participant's mean speed) and a within-participant part (deviation of
#      each response from the participant's mean)
#   4. Descriptive RT analyses                         -> Appendix C (Table C1, Figure C)
#   5. Model 1 with RT covariates (matched sample, errors only):
#        model1rt_ref         : diagnosis only, same observations as the RT
#                               models (reference for the RT models)
#        model1rt_additive    : diagnosis + RT between + RT within
#        model1rt_interaction : + diagnosis x RT between + diagnosis x RT within
#      Sensitivity: the RT models refitted without responses whose RT is above
#      the 99th percentile of the matched sample (covariates recomputed).
#   6. Posterior contrasts, appendix table (Table C2) and exports
#   7. Model 1 with a random-choice (mixture) component (chance-responding
#      benchmark): models A, B (errors only), C, C0 (all answered items), and
#      C and C0 with random choice allowed to favour the D distractor, LOO
#      (Appendix D, Tables D1-D2)
# Set MATRIKS_REUSE_FITS=1 to reuse the model1rt_*_fit.rds files of a previous run
# instead of refitting (post-processing only).
# Outputs (in outputs/): unsystematic_responding_descriptives.xlsx,
#   FigureC_response_times.png, unsystematic_responding_models.xlsx,
#   model1rt_*_fit.rds, diagnostics/model1rt_*.png, random_choice_mixture_models.xlsx,
#   model1mix_ref_fit.rds, mixture_*_fit.rds
###############################################################################

if (!exists("MATRIKS_SETUP_DONE")) source(file.path("R", "00_setup.R"))
require_packages(c(BAYES_PACKAGES, "splines", "loo"))
suppressPackageStartupMessages({ library(cmdstanr); library(posterior); library(bayesplot) })

OUTCOME_LEVELS      <- c("R", "D", "IC", "WP")          # R = reference category
NONREFERENCE_LEVELS <- c("D vs R", "IC vs R", "WP vs R")
DIAG_DIR <- file.path(OUTPUT_DIR, "diagnostics")
dir.create(DIAG_DIR, showWarnings = FALSE, recursive = TRUE)

CONTROL_TIME_STUB <- "time_mat50_item."     # RT in seconds (integer), per item
PATIENT_TIME_STUB <- "time_mat88_item."
RT_UPPER_QUANTILE <- 0.99                   # sensitivity analysis: RT above this quantile removed

# ---- Helpers (classification identical to script 03) --------------------------------
harmonize_diagnosis <- function(x) {
  x_chr <- str_to_lower(trimws(as.character(x)))
  case_when(
    x_chr == "control"                                ~ "Control",
    x_chr == "bipolar"                                ~ "Bipolar",
    x_chr %in% c("schizophrenia", "schizophrenic")    ~ "Schizophrenia",
    TRUE                                              ~ as.character(x)
  )
}

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

# Long format: one row per participant x item, with the response and its RT.
# `item` is the harmonised index 1..18, identical for the two versions.
extract_item_responses_rt <- function(data, id_col, dx_col, item_stub, time_stub, item_ids,
                                      dataset_label) {
  item_cols <- paste0(item_stub, item_ids)
  time_cols <- paste0(time_stub, item_ids)
  check_columns(data, c(item_cols, time_cols), dataset_label)
  out <- data %>%
    transmute(id        = as.character(.data[[id_col]]),
              diagnosis = harmonize_diagnosis(.data[[dx_col]]),
              across(all_of(item_cols), as.character),
              across(all_of(time_cols), to_num)) %>%
    pivot_longer(-c(id, diagnosis), names_to = c(".value", "item_raw"),
                 names_pattern = "^(response|time)_mat\\d+_item\\.(\\d+)$") %>%
    rename(response_raw = response, rt = time) %>%
    mutate(item_raw       = as.integer(item_raw),
           item           = match(item_raw, item_ids),
           source_dataset = dataset_label,
           subject_key    = paste(source_dataset, id, sep = "::"),
           response_type  = classify_response(response_raw, RIC_POLICY))
  filter(out, !is.na(response_type), response_type != "UNKNOWN")
}

summarise_participant_items <- function(item_level_data) {
  item_level_data %>%
    group_by(id, diagnosis, source_dataset) %>%
    summarise(total_items   = n(),
              total_correct = sum(response_type == "correct"),
              total_errors  = sum(response_type != "correct"),
              .groups = "drop") %>%
    mutate(total_accuracy = total_correct / total_items)
}

# ---- 1. Item-level responses and RT -------------------------------------------------
controls <- load_controls()
patients <- load_patients()

item_level_all <- bind_rows(
  extract_item_responses_rt(controls, "external_code", "Diagnosi", CONTROL_ITEM_STUB,
                            CONTROL_TIME_STUB, CONTROL_ITEMS, "adult"),
  extract_item_responses_rt(patients, "external_code", "Diagnosi", PATIENT_ITEM_STUB,
                            PATIENT_TIME_STUB, PATIENT_ITEMS, "psy")
) %>%
  filter(diagnosis %in% c("Control", "Bipolar", "Schizophrenia")) %>%
  mutate(is_error = response_type != "correct")

# ---- 2. Matched sample (from script 03) -----------------------------------------------
# The matched participants are those selected by script 03 (match_three_groups()
# in 00_setup.R: optimal matching in triplets on accuracy, caliper 0.10), read
# from the file written by script 03 so that both scripts use the same sample.
participant_summary_all <- summarise_participant_items(item_level_all) %>%
  mutate(subject_key = paste(source_dataset, id, sep = "::"))
matched_file <- file.path(OUTPUT_DIR, "matched_participants.csv")
if (!file.exists(matched_file)) stop("Run script 03 first: ", matched_file, " not found.", call. = FALSE)
matched_file_data <- read.csv(matched_file, stringsAsFactors = FALSE)
matched_participants <- matched_file_data %>%
  select(subject_key, match_set, accuracy_03 = total_accuracy) %>%
  inner_join(participant_summary_all, by = "subject_key")
stopifnot(nrow(matched_participants) == nrow(matched_file_data),
          isTRUE(all.equal(matched_participants$accuracy_03, matched_participants$total_accuracy)))
matched_item_data <- item_level_all %>%
  semi_join(select(matched_participants, subject_key), by = "subject_key")

cat("\nResponses without RT (excluded from the RT models):\n")
print(item_level_all %>% filter(is.na(rt)) %>%
        count(diagnosis, id, matched = subject_key %in% matched_participants$subject_key,
              is_error))

# ---- 3. RT covariates -----------------------------------------------------------------
# rt_c       = (RT - mean RT of that item in the matched controls) / SD of the
#              item-centred RT of the matched controls; 0 = the typical time of a
#              matched control on that item (the reference group of Model 1), one
#              unit = one SD of the matched controls (in seconds: `rt_scale`).
# rt_between = participant's mean rt_c over all answered items (correct and
#              errors): the participant's general speed.
# rt_within  = rt_c - rt_between: how much faster (< 0) or slower (> 0) than
#              usual for that participant a given response was.
# Responses without RT, or above `max_rt` (sensitivity), are removed first.
add_rt_covariates <- function(item_data, reference_keys, max_rt = Inf) {
  d <- filter(item_data, !is.na(rt), rt <= max_rt)
  is_ref <- d$subject_key %in% reference_keys
  item_means <- d[is_ref, ] %>%
    group_by(item) %>% summarise(item_mean_rt = mean(rt), .groups = "drop")
  d <- d %>% left_join(item_means, by = "item") %>% mutate(rt_item_c = rt - item_mean_rt)
  rt_scale <- sd(d$rt_item_c[is_ref])
  d <- d %>%
    mutate(rt_c = rt_item_c / rt_scale) %>%
    group_by(subject_key) %>%
    mutate(rt_between = mean(rt_c), rt_within = rt_c - rt_between) %>%
    ungroup()
  list(data = d, rt_scale = rt_scale, item_means = item_means, max_rt = max_rt)
}

matched_control_keys <- matched_participants$subject_key[matched_participants$diagnosis == "Control"]
rt_main <- add_rt_covariates(item_level_all, matched_control_keys)
rt_cut  <- unname(quantile(matched_item_data$rt, RT_UPPER_QUANTILE, na.rm = TRUE))
rt_sens <- add_rt_covariates(item_level_all, matched_control_keys, max_rt = rt_cut)
cat(sprintf("\nRT scale (SD of item-centred RT in matched controls): %.1f s; %s percentile cut-off for the sensitivity analysis: %.0f s\n",
            rt_main$rt_scale, RT_UPPER_QUANTILE * 100, rt_cut))

# ---- 4. Descriptive RT analyses (Appendix C) ----------------------------------------------
rt_stats <- function(x) {
  x <- x[!is.na(x)]
  tibble(n = length(x), mean = mean(x), sd = sd(x), min = min(x),
         q10 = quantile(x, .10), q25 = quantile(x, .25), median = median(x),
         q75 = quantile(x, .75), q90 = quantile(x, .90), q99 = quantile(x, .99), max = max(x))
}

rt_samples <- bind_rows(mutate(item_level_all,    sample = "full"),
                        mutate(matched_item_data, sample = "matched"))

desc_by_group <- rt_samples %>%
  group_by(sample, diagnosis) %>%
  group_modify(~ bind_cols(rt_stats(.x$rt), tibble(n_missing_rt = sum(is.na(.x$rt))))) %>%
  ungroup()
desc_by_response <- rt_samples %>%
  mutate(response = if_else(is_error, paste("error", response_type), "correct")) %>%
  group_by(sample, diagnosis, response) %>%
  group_modify(~ rt_stats(.x$rt)) %>%
  ungroup()
desc_participant_speed <- rt_main$data %>%
  semi_join(select(matched_participants, subject_key), by = "subject_key") %>%
  distinct(subject_key, diagnosis, rt_between) %>%
  group_by(diagnosis) %>%
  summarise(n = n(), mean_rt_between = mean(rt_between), sd_rt_between = sd(rt_between),
            median_rt_between = median(rt_between), min = min(rt_between), max = max(rt_between),
            .groups = "drop")
print(desc_by_group, width = Inf)
print(desc_participant_speed, width = Inf)

write_xlsx(list(rt_by_group           = desc_by_group,
                rt_by_response_type   = desc_by_response,
                participant_speed     = desc_participant_speed,
                item_mean_rt_matched_controls = rt_main$item_means,
                rt_scale              = tibble(rt_scale_seconds = rt_main$rt_scale,
                                               sensitivity_cutoff_seconds = rt_cut)),
           file.path(OUTPUT_DIR, "unsystematic_responding_descriptives.xlsx"))

# Figure C: (A) cumulative distribution of RT by group; (B) proportion of
# D errors among errors as a smooth function of RT (logistic regression on a
# natural spline of RT, matched sample). Axes truncated for readability only.
plot_data <- matched_item_data %>%
  filter(!is.na(rt)) %>%
  mutate(diagnosis = recode(diagnosis, "Schizophrenia" = "Schizophrenic"))

fig_a <- ggplot(plot_data, aes(x = rt, color = diagnosis)) +
  stat_ecdf(linewidth = 1) +
  coord_cartesian(xlim = c(0, 120)) +
  scale_color_manual(values = DX_COLORS) +
  labs(title = "A. Response times (matched sample)", x = "Response time (s)",
       y = "Cumulative proportion of responses", color = "Conditions") +
  DX_THEME

fig_b <- ggplot(filter(plot_data, is_error),
                aes(x = rt, y = as.numeric(response_type == "D"), color = diagnosis, fill = diagnosis)) +
  geom_smooth(method = "glm", formula = y ~ splines::ns(x, df = 3),
              method.args = list(family = binomial), alpha = 0.15) +
  geom_rug(sides = "b", alpha = 0.3, position = position_jitter(width = 0.3, height = 0)) +
  coord_cartesian(xlim = c(0, 60), ylim = c(0, 0.4)) +
  scale_color_manual(values = DX_COLORS) + scale_fill_manual(values = DX_COLORS) +
  guides(color = "none", fill = "none") +   # legend shown once, under panel A
  labs(title = "B. D errors by response time (errors, matched sample)",
       x = "Response time (s)", y = "Proportion of D among errors",
       color = "Conditions", fill = "Conditions") +
  DX_THEME

ggsave(file.path(OUTPUT_DIR, "FigureC_response_times.png"),
       (fig_a / fig_b) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
       width = 11, height = 12, dpi = 300)

# ---- 5. Model 1 with RT covariates ------------------------------------------------------
PREDICTORS <- list(
  ref         = c("dx_bipolar", "dx_schizophrenic"),
  additive    = c("dx_bipolar", "dx_schizophrenic", "rt_between", "rt_within"),
  interaction = c("dx_bipolar", "dx_schizophrenic", "rt_between", "rt_within",
                  "bd_x_rt_between", "sz_x_rt_between", "bd_x_rt_within", "sz_x_rt_within")
)

build_model_data <- function(rt_obj) {
  rt_obj$data %>%
    semi_join(select(matched_participants, subject_key), by = "subject_key") %>%
    filter(response_type %in% OUTCOME_LEVELS) %>%
    mutate(dx_bipolar       = as.numeric(diagnosis == "Bipolar"),
           dx_schizophrenic = as.numeric(diagnosis == "Schizophrenia"),
           bd_x_rt_between  = dx_bipolar * rt_between,
           sz_x_rt_between  = dx_schizophrenic * rt_between,
           bd_x_rt_within   = dx_bipolar * rt_within,
           sz_x_rt_within   = dx_schizophrenic * rt_within)
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
    refresh = 1000, output_dir = check_stan_csv_dir(STAN_CSV_DIR), output_basename = label
  )
  fit$save_object(file.path(OUTPUT_DIR, paste0(label, "_fit.rds")))
  fit
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

model_data_main <- build_model_data(rt_main)
model_data_sens <- build_model_data(rt_sens)
cat("\nRT models - errors by group and type (main):\n")
print(table(model_data_main$diagnosis, model_data_main$response_type))
cat("\nRT models - errors by group and type (sensitivity, RT <= cut-off):\n")
print(table(model_data_sens$diagnosis, model_data_sens$response_type))

model_specs <- tribble(
  ~label,                       ~predictors,   ~data,
  "model1rt_ref",               "ref",         "main",
  "model1rt_additive",          "additive",    "main",
  "model1rt_interaction",       "interaction", "main",
  "model1rt_additive_sens",     "additive",    "sens",
  "model1rt_interaction_sens",  "interaction", "sens"
)

REUSE_FITS <- identical(Sys.getenv("MATRIKS_REUSE_FITS"), "1")
multinomial_model <- cmdstan_model(STAN_FILE)
fits <- list(); stan_data <- list()
for (i in seq_len(nrow(model_specs))) {
  s <- model_specs[i, ]
  d <- if (s$data == "main") model_data_main else model_data_sens
  stan_data[[s$label]] <- make_stan_data(d, PREDICTORS[[s$predictors]])
  rds <- file.path(OUTPUT_DIR, paste0(s$label, "_fit.rds"))
  fits[[s$label]] <- if (REUSE_FITS && file.exists(rds)) readRDS(rds) else
    fit_model(multinomial_model, stan_data[[s$label]], s$label)
}

diagnostics <- map_dfr(model_specs$label, ~ compute_model_diagnostics(fits[[.x]], .x, stan_data[[.x]]))
print(diagnostics, width = Inf)

# ---- 6. Posterior contrasts -------------------------------------------------------------
draws_mat <- function(fit, variable) {
  m <- as_draws_matrix(fit$draws(variables = variable))
  matrix(as.numeric(m), nrow = nrow(m), dimnames = list(NULL, colnames(m)))
}

get_beta_array <- function(fit, predictor_names) {
  d  <- draws_mat(fit, "beta")
  P  <- length(predictor_names); K1 <- length(NONREFERENCE_LEVELS)
  arr <- array(NA_real_, dim = c(nrow(d), P, K1),
               dimnames = list(NULL, predictor_names, NONREFERENCE_LEVELS))
  for (p in seq_len(P)) for (k in seq_len(K1)) arr[, p, k] <- d[, sprintf("beta[%d,%d]", p, k)]
  arr
}

summ_draws <- function(x) {
  qs <- unname(quantile(x, c(0.025, 0.975))); med <- median(x)
  tibble(median = med, q2.5 = qs[1], q97.5 = qs[2],
         OR_median = exp(med), OR_q2.5 = exp(qs[1]), OR_q97.5 = exp(qs[2]),
         P_gt_0 = mean(x > 0), P_lt_0 = mean(x < 0),
         cri95_excludes_0 = qs[1] > 0 | qs[2] < 0)
}

# A contrast is a named weight vector over predictors; its value for each error
# type is sum_p weight_p * beta[p, type], computed draw by draw.
contrast_summary <- function(beta, weights, label) {
  m <- matrix(0, dim(beta)[1], dim(beta)[3], dimnames = list(NULL, NONREFERENCE_LEVELS))
  for (p in names(weights)) m <- m + weights[[p]] * beta[, p, ]
  bind_rows(
    map_dfr(colnames(m), ~ bind_cols(tibble(contrast = label, error_type = .x), summ_draws(m[, .x]))),
    bind_cols(tibble(contrast = label, error_type = "Average of the 3 (vs R)"), summ_draws(rowMeans(m)))
  )
}

# Differences between error types within a contrast (contrast x error-type
# interaction): D - IC, D - WP, IC - WP.
type_diff_summary <- function(beta, weights, label) {
  m <- matrix(0, dim(beta)[1], dim(beta)[3])
  for (p in names(weights)) m <- m + weights[[p]] * beta[, p, ]
  short <- sub(" vs R$", "", NONREFERENCE_LEVELS)
  map_dfr(combn(seq_along(short), 2, simplify = FALSE), function(ix) {
    bind_cols(tibble(contrast = label, type_difference = paste(short[ix[1]], "-", short[ix[2]])),
              summ_draws(m[, ix[1]] - m[, ix[2]]))
  })
}

# Typical speed of the matched schizophrenia group (median rt_between).
sz_speed <- desc_participant_speed$median_rt_between[desc_participant_speed$diagnosis == "Schizophrenia"]

contrast_sets <- list(
  ref = list(
    "BD vs Control" = c(dx_bipolar = 1),
    "SZ vs Control" = c(dx_schizophrenic = 1),
    "SZ vs BD"      = c(dx_schizophrenic = 1, dx_bipolar = -1)),
  additive = list(
    "BD vs Control (equal RT)"          = c(dx_bipolar = 1),
    "SZ vs Control (equal RT)"          = c(dx_schizophrenic = 1),
    "SZ vs BD (equal RT)"               = c(dx_schizophrenic = 1, dx_bipolar = -1),
    "RT between (+1 SD slower person)"  = c(rt_between = 1),
    "RT within (+1 SD slower response)" = c(rt_within = 1)),
  interaction = list(
    "BD vs Control (RT at matched-control item mean)" = c(dx_bipolar = 1),
    "SZ vs Control (RT at matched-control item mean)" = c(dx_schizophrenic = 1),
    "SZ vs BD (RT at matched-control item mean)"      = c(dx_schizophrenic = 1, dx_bipolar = -1),
    "SZ vs Control (person speed = SZ median)" = c(dx_schizophrenic = 1, sz_x_rt_between = sz_speed),
    "RT between slope: Control"  = c(rt_between = 1),
    "RT between slope: BD"       = c(rt_between = 1, bd_x_rt_between = 1),
    "RT between slope: SZ"       = c(rt_between = 1, sz_x_rt_between = 1),
    "RT between: SZ - Control (interaction)" = c(sz_x_rt_between = 1),
    "RT between: BD - Control (interaction)" = c(bd_x_rt_between = 1),
    "RT within slope: Control"   = c(rt_within = 1),
    "RT within slope: BD"        = c(rt_within = 1, bd_x_rt_within = 1),
    "RT within slope: SZ"        = c(rt_within = 1, sz_x_rt_within = 1),
    "RT within: SZ - Control (interaction)" = c(sz_x_rt_within = 1),
    "RT within: BD - Control (interaction)" = c(bd_x_rt_within = 1))
)
GROUP_CONTRASTS <- c("BD vs Control", "SZ vs Control", "SZ vs BD")   # first three of each set

contrasts_all <- map_dfr(seq_len(nrow(model_specs)), function(i) {
  s <- model_specs[i, ]
  beta <- get_beta_array(fits[[s$label]], PREDICTORS[[s$predictors]])
  imap_dfr(contrast_sets[[s$predictors]], ~ contrast_summary(beta, .x, .y)) %>%
    mutate(model = s$label, .before = 1)
})
print(select(contrasts_all, model, contrast, error_type, OR_median, OR_q2.5, OR_q97.5, P_gt_0),
      n = Inf, width = Inf)

type_diffs_all <- map_dfr(seq_len(nrow(model_specs)), function(i) {
  s <- model_specs[i, ]
  beta <- get_beta_array(fits[[s$label]], PREDICTORS[[s$predictors]])
  imap_dfr(contrast_sets[[s$predictors]][1:3], ~ type_diff_summary(beta, .x, .y)) %>%
    mutate(model = s$label, .before = 1)
})

pooling <- map_dfr(seq_len(nrow(model_specs)), function(i) {
  s <- model_specs[i, ]
  map_dfr(c("mu_beta", "tau_beta"), function(v) {
    m <- draws_mat(fits[[s$label]], v)
    map_dfr(seq_len(ncol(m)), ~ bind_cols(tibble(model = s$label, parameter = v,
                                                 term = PREDICTORS[[s$predictors]][.x]),
                                          summ_draws(m[, .x])))
  })
})

random_sd <- map_dfr(model_specs$label, function(lab) {
  map_dfr(c("sigma_u", "sigma_v", "alpha"), function(v) {
    m <- draws_mat(fits[[lab]], v)
    map_dfr(seq_len(ncol(m)), ~ bind_cols(tibble(model = lab, parameter = v,
                                                 comparison = NONREFERENCE_LEVELS[.x]),
                                          summ_draws(m[, .x])))
  })
})

for (lab in model_specs$label) {
  p_names <- PREDICTORS[[model_specs$predictors[model_specs$label == lab]]]
  beta_draws <- fits[[lab]]$draws(variables = "beta")
  idx <- str_match(dimnames(beta_draws)$variable, "beta\\[(\\d+),(\\d+)\\]")
  dimnames(beta_draws)$variable <- paste0(p_names[as.integer(idx[, 2])], " | ",
                                          NONREFERENCE_LEVELS[as.integer(idx[, 3])])
  ggsave(file.path(DIAG_DIR, paste0(lab, "_intervals.png")),
         mcmc_intervals(beta_draws, prob = 0.8, prob_outer = 0.95) +
           ggtitle(sprintf("%s coefficients (log-odds vs R)", lab)),
         width = 11, height = max(5, 0.35 * length(p_names) * 3), dpi = 300)
}

# ---- Appendix C, Table C2 (same layout as Table A1) ----------------------------------
fmt_num <- function(x) sub("^-", "\u2212", formatC(round(x, 2) + 0, format = "f", digits = 2))
fmt_ci  <- function(lo, hi) sprintf("[%s, %s]", fmt_num(lo), fmt_num(hi))
fmt_pp  <- function(x) ifelse(x > .995, ">.99", ifelse(x < .005, "<.01",
                                                    sub("^0", "", formatC(x, format = "f", digits = 2))))
table_row <- function(d, label, or = TRUE, prob = TRUE) {
  tibble(Parameter = label, Estimate = fmt_num(d$median), `95% CrI` = fmt_ci(d$q2.5, d$q97.5),
         OR = if (or) fmt_num(d$OR_median) else "\u2013",
         `OR 95% CrI` = if (or) fmt_ci(d$OR_q2.5, d$OR_q97.5) else "\u2013",
         P_beta_gt_0 = if (prob) fmt_pp(d$P_gt_0) else "\u2013")
}
heading <- function(label) tibble(Parameter = label, Estimate = "", `95% CrI` = "", OR = "",
                                  `OR 95% CrI` = "", P_beta_gt_0 = "")
type_rows <- function(d, col) map_dfr(seq_len(nrow(d)), ~ table_row(d[.x, ], paste0("    ", sub(" - ", " \u2212 ", d[[col]][.x]))))
contrast_block <- function(model, contrast, label, average = FALSE) {
  d <- filter(contrasts_all, model == !!model, contrast == !!contrast)
  if (!average) d <- filter(d, error_type != "Average of the 3 (vs R)")
  d$error_type <- sub("Average of the 3 \\(vs R\\)", "Average of the three", d$error_type)
  bind_rows(heading(paste0("  ", label)), type_rows(d, "error_type"))
}
param_rows <- function(model, parameter, labels, or, prob) {
  d <- filter(random_sd, model == !!model, parameter == !!parameter)
  map_dfr(seq_len(nrow(d)), ~ table_row(d[.x, ], paste0("    ", labels[.x]), or = or, prob = prob))
}
pool_rows <- function(model, parameter, labels, prob) {
  d <- filter(pooling, model == !!model, parameter == !!parameter)
  map_dfr(seq_len(nrow(d)), ~ table_row(d[.x, ], paste0("    ", labels[d$term[.x]]), or = FALSE, prob = prob))
}
term_labels <- c(dx_bipolar = "Bipolar", dx_schizophrenic = "Schizophrenia",
                 rt_between = "RT between", rt_within = "RT within")
ia_sets <- contrast_sets$interaction

appendix_table <- bind_rows(
  heading("A. Group contrasts at equal response time (additive model): log-odds of each error type relative to R"),
  contrast_block("model1rt_additive", "BD vs Control (equal RT)", "Bipolar vs controls"),
  contrast_block("model1rt_additive", "SZ vs Control (equal RT)", "Schizophrenia vs controls"),
  contrast_block("model1rt_additive", "SZ vs BD (equal RT)", "Schizophrenia vs bipolar", average = TRUE),
  heading("B. Differences between error types within each contrast (additive model)"),
  map_dfr(c("BD vs Control (equal RT)" = "Bipolar vs controls",
            "SZ vs Control (equal RT)" = "Schizophrenia vs controls",
            "SZ vs BD (equal RT)" = "Schizophrenia vs bipolar") %>% imap(~ c(.y, .x)),
          ~ bind_rows(heading(paste0("  ", .x[2])),
                      type_rows(filter(type_diffs_all, model == "model1rt_additive", contrast == .x[1]),
                                "type_difference"))),
  heading(sprintf("C. Response-time effects (additive model), per SD (%.0f s) slower", rt_main$rt_scale)),
  contrast_block("model1rt_additive", "RT between (+1 SD slower person)", "Between-participant RT (participant's mean speed)"),
  contrast_block("model1rt_additive", "RT within (+1 SD slower response)", "Within-participant RT (response vs participant's mean)"),
  heading("D. Diagnosis x response-time interaction model"),
  contrast_block("model1rt_interaction", "SZ vs Control (RT at matched-control item mean)", "Schizophrenia vs controls, at the typical RT of matched controls"),
  contrast_block("model1rt_interaction", "SZ vs Control (person speed = SZ median)", "Schizophrenia vs controls, at the typical speed of schizophrenia patients"),
  contrast_block("model1rt_interaction", "BD vs Control (RT at matched-control item mean)", "Bipolar vs controls, at the typical RT of matched controls"),
  contrast_block("model1rt_interaction", "RT within slope: Control", "Within-participant RT slope, controls"),
  contrast_block("model1rt_interaction", "RT within slope: BD", "Within-participant RT slope, bipolar"),
  contrast_block("model1rt_interaction", "RT within slope: SZ", "Within-participant RT slope, schizophrenia"),
  contrast_block("model1rt_interaction", "RT within: SZ - Control (interaction)", "Schizophrenia x within-participant RT (vs controls)"),
  contrast_block("model1rt_interaction", "RT between slope: Control", "Between-participant RT slope, controls"),
  contrast_block("model1rt_interaction", "RT between slope: BD", "Between-participant RT slope, bipolar"),
  contrast_block("model1rt_interaction", "RT between slope: SZ", "Between-participant RT slope, schizophrenia"),
  contrast_block("model1rt_interaction", "RT between: SZ - Control (interaction)", "Schizophrenia x between-participant RT (vs controls)"),
  heading("E. Baseline, pooling and variance parameters (additive model)"),
  heading("  Intercepts: controls at the typical RT of matched controls, log-odds relative to R"),
  param_rows("model1rt_additive", "alpha", NONREFERENCE_LEVELS, or = TRUE, prob = FALSE),
  heading("  Pooling across the three contrasts"),
  pool_rows("model1rt_additive", "mu_beta", setNames(paste("\u03bc,", term_labels), names(term_labels)), prob = TRUE),
  pool_rows("model1rt_additive", "tau_beta", setNames(paste("\u03c4,", term_labels), names(term_labels)), prob = FALSE),
  heading("  Random-intercept SDs (log-odds scale)"),
  param_rows("model1rt_additive", "sigma_u", paste("Participants,", NONREFERENCE_LEVELS), or = FALSE, prob = FALSE),
  param_rows("model1rt_additive", "sigma_v", paste("Items,", NONREFERENCE_LEVELS), or = FALSE, prob = FALSE)
)
names(appendix_table)[names(appendix_table) == "P_beta_gt_0"] <- "P(\u03b2 > 0)"

write_xlsx(
  list(table_C2_appendix = appendix_table,
       contrasts         = contrasts_all,
       type_differences  = type_diffs_all,
       pooling           = pooling,
       random_effects    = random_sd,
       diagnostics       = diagnostics,
       rt_settings       = tibble(rt_scale_seconds = rt_main$rt_scale,
                                  sensitivity_cutoff_seconds = rt_cut,
                                  sz_median_person_speed = sz_speed)),
  file.path(OUTPUT_DIR, "unsystematic_responding_models.xlsx")
)

# ---- 7. Model 1 with a random-choice (mixture) component ----------------------------------
# Chance-responding benchmark. Each response is modelled as
# a mixture of the systematic process of Model 1 (probability 1 - lambda) and of
# a random choice among the response options of the item (probability lambda),
# under which each error type is chosen in proportion to its number of
# distractors. lambda differs by group (Control, BD, SZ) and between participants
# (random effect on the logit scale). No RT covariates. Matched sample, all
# participants (same errors as Model 1). See stan/mixture_random_choice.stan.
#   model1mix_ref         : Model 1 (no random component), for the LOO comparison
#   mixture_errors_dx     : errors only, diagnosis in the systematic part   (model A)
#   mixture_errors_chance : errors only, no diagnosis in the systematic part,
#                           group differences only through lambda           (model B)
#   mixture_all_dx        : all answered items (correct responses included,
#                           random choice over the 8 options)               (model C)
#   mixture_all_chance    : as C, without diagnosis in the error-type part  (model C0)
#   mixture_all_dx_salience, mixture_all_chance_salience: as C and C0, but random
#                           choice may favour the D distractor, which has weight
#                           exp(delta) relative to the other options (perceptual
#                           salience of D; delta = 0 in C and C0)        (models Cs, C0s)
MIXTURE_STAN_FILE <- file.path("stan", "mixture_random_choice.stan")
MIX_GROUPS        <- c("Control", "Bipolar", "Schizophrenia")
N_DISTRACTORS     <- 7

# Distractor composition of each item (columns in OUTCOME_LEVELS order: R, D, IC,
# WP). Two R, one D, two IC and two WP distractors, except harmonised item 1
# (MatriKS88 34 / MatriKS 12+ 4): one R, one D, three IC, two WP. In the three
# items in which repetition errors are recorded with the "r.ic" label, both
# repetition distractors are counted as R, consistent with RIC_POLICY = "R".
DISTRACTOR_COUNTS <- matrix(rep(c(2, 1, 2, 2), each = length(PATIENT_ITEMS)),
                            ncol = length(OUTCOME_LEVELS), dimnames = list(NULL, OUTCOME_LEVELS))
DISTRACTOR_COUNTS[match(34, PATIENT_ITEMS), ] <- c(1, 1, 3, 2)
if (RIC_POLICY != "R") stop("Script 07 requires RIC_POLICY = \"R\" (the default).", call. = FALSE)
stopifnot(all(rowSums(DISTRACTOR_COUNTS) == N_DISTRACTORS))

mix_item_data <- matched_item_data %>%
  mutate(dx_bipolar       = as.numeric(diagnosis == "Bipolar"),
         dx_schizophrenic = as.numeric(diagnosis == "Schizophrenia"))
mix_error_data <- filter(mix_item_data, response_type %in% OUTCOME_LEVELS)

make_mixture_data <- function(d, predictors, include_correct, salient = 0L) {
  keys  <- sort(unique(d$subject_key))
  group <- d$diagnosis[match(keys, d$subject_key)]
  X  <- matrix(as.numeric(as.matrix(d[, predictors])), nrow = nrow(d), ncol = length(predictors))
  Xa <- if (include_correct) as.matrix(d[, MODEL1_PREDICTORS_MIX]) else matrix(0, nrow(d), 0)
  list(N = nrow(d), K = length(OUTCOME_LEVELS), P = length(predictors),
       J = length(keys), I = length(PATIENT_ITEMS), G = length(MIX_GROUPS),
       subject = match(d$subject_key, keys), item = as.integer(d$item),
       subject_group = match(group, MIX_GROUPS),
       include_correct = as.integer(include_correct),
       y = ifelse(d$response_type == "correct", 0L, match(d$response_type, OUTCOME_LEVELS)),
       X = X, Pa = ncol(Xa), Xa = Xa,
       q = DISTRACTOR_COUNTS / N_DISTRACTORS, D = N_DISTRACTORS,
       salient = as.integer(salient))
}
MODEL1_PREDICTORS_MIX <- c("dx_bipolar", "dx_schizophrenic")

SALIENT_TYPE <- match("D", OUTCOME_LEVELS)
mixture_specs <- tribble(
  ~label,                        ~predictors,                 ~include_correct, ~salient,
  "mixture_errors_dx",           list(MODEL1_PREDICTORS_MIX), FALSE,            0L,
  "mixture_errors_chance",       list(character(0)),          FALSE,            0L,
  "mixture_all_dx",              list(MODEL1_PREDICTORS_MIX), TRUE,             0L,
  "mixture_all_chance",          list(character(0)),          TRUE,             0L,
  "mixture_all_dx_salience",     list(MODEL1_PREDICTORS_MIX), TRUE,             SALIENT_TYPE,
  "mixture_all_chance_salience", list(character(0)),          TRUE,             SALIENT_TYPE
)

mixture_model <- cmdstan_model(MIXTURE_STAN_FILE)
mix_fits <- list(); mix_data <- list()
get_fit <- function(label, fit_fun) {
  rds <- file.path(OUTPUT_DIR, paste0(label, "_fit.rds"))
  if (REUSE_FITS && file.exists(rds)) readRDS(rds) else fit_fun()
}

# Model 1 refitted here (same data as script 03) to obtain its log-likelihood.
mix_data$model1mix_ref <- make_stan_data(mix_error_data, MODEL1_PREDICTORS_MIX)
mix_fits$model1mix_ref <- get_fit("model1mix_ref", function()
  fit_model(multinomial_model, mix_data$model1mix_ref, "model1mix_ref"))
for (i in seq_len(nrow(mixture_specs))) {
  s <- mixture_specs[i, ]
  d <- if (s$include_correct) mix_item_data else mix_error_data
  mix_data[[s$label]] <- make_mixture_data(d, unlist(s$predictors), s$include_correct, s$salient)
  mix_fits[[s$label]] <- get_fit(s$label, function()
    fit_model(mixture_model, mix_data[[s$label]], s$label))
}

mixture_diagnostics <- map_dfr(names(mix_fits), function(lab) {
  fit <- mix_fits[[lab]]
  vars <- intersect(c("alpha", "mu_beta", "tau_beta", "beta", "sigma_u", "sigma_v", "gamma",
                      "sigma_w", "a0", "b_acc", "sigma_ua", "sigma_va", "delta"),
                    fit$metadata()$stan_variables)
  vars <- vars[map_lgl(vars, ~ prod(fit$metadata()$stan_variable_sizes[[.x]]) > 0)]
  summ <- fit$summary(variables = vars)
  dsum <- fit$diagnostic_summary(quiet = TRUE)
  tibble(model = lab, observations = mix_data[[lab]]$N, participants = mix_data[[lab]]$J,
         max_rhat = max(summ$rhat, na.rm = TRUE), min_bulk_ess = min(summ$ess_bulk, na.rm = TRUE),
         min_tail_ess = min(summ$ess_tail, na.rm = TRUE), divergences = sum(dsum$num_divergent),
         treedepth_hits = sum(dsum$num_max_treedepth), min_ebfmi = min(dsum$ebfmi, na.rm = TRUE))
})
print(mixture_diagnostics, width = Inf)

# Group contrasts in the systematic component (models with diagnosis).
mixture_contrasts <- map_dfr(c("model1mix_ref", "mixture_errors_dx", "mixture_all_dx",
                               "mixture_all_dx_salience"), function(lab) {
  beta <- get_beta_array(mix_fits[[lab]], MODEL1_PREDICTORS_MIX)
  bind_rows(
    imap_dfr(contrast_sets$ref, ~ contrast_summary(beta, .x, .y)),
    imap_dfr(contrast_sets$ref, ~ type_diff_summary(beta, .x, .y)) %>%
      rename(error_type = type_difference)
  ) %>% mutate(model = lab, .before = 1)
})
print(select(mixture_contrasts, model, contrast, error_type, OR_median, OR_q2.5, OR_q97.5, P_gt_0),
      n = Inf, width = Inf)

# Random-choice share: lambda of a typical participant and mean lambda of each
# group, their differences, and the expected share of errors attributed to the
# random component (posterior responsibility).
n_errors_by_group <- table(factor(mix_error_data$diagnosis, levels = MIX_GROUPS))
mixture_lambda <- map_dfr(mixture_specs$label, function(lab) {
  fit <- mix_fits[[lab]]
  lt <- draws_mat(fit, "lambda_typical"); lm <- draws_mat(fit, "lambda_mean")
  re <- sweep(draws_mat(fit, "random_errors"), 2, as.numeric(n_errors_by_group), "/")
  pieces <- list()
  for (g in seq_along(MIX_GROUPS)) {
    pieces[[length(pieces) + 1]] <- bind_cols(tibble(quantity = "lambda, typical participant", group = MIX_GROUPS[g]), summ_draws(lt[, g]))
    pieces[[length(pieces) + 1]] <- bind_cols(tibble(quantity = "lambda, group mean", group = MIX_GROUPS[g]), summ_draws(lm[, g]))
    pieces[[length(pieces) + 1]] <- bind_cols(tibble(quantity = "share of errors from random choice", group = MIX_GROUPS[g]), summ_draws(re[, g]))
  }
  for (pr in list(c(3, 1), c(2, 1), c(3, 2))) {
    lab_pr <- paste(MIX_GROUPS[pr[1]], "-", MIX_GROUPS[pr[2]])
    pieces[[length(pieces) + 1]] <- bind_cols(tibble(quantity = "lambda, group mean: difference", group = lab_pr), summ_draws(lm[, pr[1]] - lm[, pr[2]]))
    pieces[[length(pieces) + 1]] <- bind_cols(tibble(quantity = "share of errors from random choice: difference", group = lab_pr), summ_draws(re[, pr[1]] - re[, pr[2]]))
  }
  bind_rows(pieces) %>% select(-starts_with("OR_"), -cri95_excludes_0) %>% mutate(model = lab, .before = 1)
})
print(mixture_lambda, n = Inf, width = Inf)

# Posterior expected counts vs observed counts by group and response type.
observed_counts <- mix_item_data %>%
  count(diagnosis, response_type) %>%
  mutate(response_type = recode(response_type, correct = "correct"))
mixture_counts <- map_dfr(mixture_specs$label, function(lab) {
  ec <- draws_mat(mix_fits[[lab]], "expected_counts")
  types <- c("correct", OUTCOME_LEVELS)
  map_dfr(seq_along(MIX_GROUPS), function(g) map_dfr(seq_along(types), function(k) {
    x <- ec[, sprintf("expected_counts[%d,%d]", g, k)]
    if (!mixture_specs$include_correct[mixture_specs$label == lab] && k == 1) return(NULL)
    tibble(model = lab, group = MIX_GROUPS[g], response_type = types[k],
           expected_median = median(x), expected_q2.5 = quantile(x, .025), expected_q97.5 = quantile(x, .975))
  }))
}) %>%
  left_join(rename(observed_counts, group = diagnosis, observed = n), by = c("group", "response_type"))
print(mixture_counts, n = Inf, width = Inf)

# Accuracy part of models C and C0 (log-odds of a correct systematic response).
mixture_accuracy <- map_dfr(grep("^mixture_all", mixture_specs$label, value = TRUE), function(lab) {
  m <- draws_mat(mix_fits[[lab]], "b_acc")
  map_dfr(seq_len(ncol(m)), ~ bind_cols(tibble(model = lab, term = MODEL1_PREDICTORS_MIX[.x]), summ_draws(m[, .x])))
})

# Salience of D in random choice: delta (log weight of the D distractor relative
# to each other option; exp(delta) = 1 for uniform random choice).
mixture_salience <- map_dfr(mixture_specs$label[mixture_specs$salient > 0], function(lab) {
  bind_cols(tibble(model = lab, term = "delta (D weight in random choice)"),
            summ_draws(draws_mat(mix_fits[[lab]], "delta")[, 1]))
})
print(mixture_salience, width = Inf)

# LOO comparison: the errors-only models share the same data (Model 1, A, B);
# models C and C0 are compared with each other.
loo_list <- map(names(mix_fits), ~ mix_fits[[.x]]$loo(variables = "log_lik", cores = MCMC$parallel_chains))
names(loo_list) <- names(mix_fits)
loo_compare_table <- function(labels) {
  cmp <- loo::loo_compare(loo_list[labels])
  cols <- c("elpd_diff", "se_diff", "elpd_loo", "se_elpd_loo", "p_loo", "looic")
  cmp_df <- as.data.frame(unclass(cmp))                 # matrix (loo < 2.9) or data frame
  model_names <- if ("model" %in% names(cmp_df)) cmp_df$model else rownames(cmp)
  tibble(model = model_names, !!!map(set_names(cols), ~ as.numeric(cmp_df[[.x]]))) %>%
    mutate(pareto_k_above_0.7 = map_int(model, ~ sum(loo_list[[.x]]$diagnostics$pareto_k > 0.7)))
}
mixture_loo <- bind_rows(
  mutate(loo_compare_table(c("model1mix_ref", "mixture_errors_dx", "mixture_errors_chance")), set = "errors only", .before = 1),
  mutate(loo_compare_table(grep("^mixture_all", mixture_specs$label, value = TRUE)),
         set = "all responses", .before = 1)
)
print(mixture_loo, width = Inf)

# ---- Appendix D, Table D1 (same layout as Tables A1 and C2) --------------------------
mix_names <- c(model1mix_ref = "Model 1 (no random-choice component)",
               mixture_errors_dx = "Model A: errors only, diagnosis in the systematic component",
               mixture_errors_chance = "Model B: errors only, chance-only account",
               mixture_all_dx = "Model C: all responses, diagnosis in the systematic component",
               mixture_all_chance = "Model C0: all responses, chance-only account",
               mixture_all_dx_salience = "Model Cs: model C with salience of D in random choice",
               mixture_all_chance_salience = "Model C0s: model C0 with salience of D in random choice")
mix_block <- function(model, contrast, label, types) {
  d <- filter(mixture_contrasts, model == !!model, contrast == !!contrast, error_type %in% types)
  d$error_type <- sub("Average of the 3 \\(vs R\\)", "Average of the three", d$error_type)
  bind_rows(heading(paste0("  ", label)), type_rows(d, "error_type"))
}
contrast_labels <- c("BD vs Control" = "Bipolar vs controls", "SZ vs Control" = "Schizophrenia vs controls",
                     "SZ vs BD" = "Schizophrenia vs bipolar")
by_type_rows <- function(model) imap_dfr(contrast_labels, ~ mix_block(model, .y, .x,
  c(NONREFERENCE_LEVELS, if (.y == "SZ vs BD") "Average of the 3 (vs R)")))
diff_rows <- function(model) imap_dfr(contrast_labels, ~ mix_block(model, .y, .x, c("D - IC", "D - WP", "IC - WP")))
group_labels <- c(Control = "Controls", Bipolar = "Bipolar", Schizophrenia = "Schizophrenia",
                  "Bipolar - Control" = "Bipolar \u2212 controls",
                  "Schizophrenia - Control" = "Schizophrenia \u2212 controls",
                  "Schizophrenia - Bipolar" = "Schizophrenia \u2212 bipolar")
lambda_rows <- function(model, quantity) {
  d <- filter(mixture_lambda, model == !!model, quantity %in% c(!!quantity, paste0(!!quantity, ": difference")))
  map_dfr(seq_len(nrow(d)), ~ table_row(d[.x, ], paste0("    ", group_labels[[d$group[.x]]]), or = FALSE,
                                        prob = grepl("difference", d$quantity[.x])))
}
mix_param_rows <- function(model, variable, labels, or = FALSE, prob = FALSE) {
  m <- draws_mat(mix_fits[[model]], variable)
  map_dfr(seq_len(ncol(m)), ~ table_row(summ_draws(m[, .x]), paste0("    ", labels[.x]), or = or, prob = prob))
}

table_D1 <- bind_rows(
  heading("A. Group contrasts in the systematic component, model A (errors only): log-odds of each error type relative to R"),
  by_type_rows("mixture_errors_dx"),
  heading("B. Differences between error types within each contrast, model A"),
  diff_rows("mixture_errors_dx"),
  heading("C. Group contrasts in the systematic component, model C (all responses)"),
  by_type_rows("mixture_all_dx"),
  heading("D. Random-choice component, model A (errors only)"),
  heading("  Share of errors attributed to random choice"),
  lambda_rows("mixture_errors_dx", "share of errors from random choice"),
  heading("  Mixing weight \u03bb (mean of the participants of each group)"),
  lambda_rows("mixture_errors_dx", "lambda, group mean"),
  heading("E. Random-choice component, model C (all responses)"),
  heading("  Share of errors attributed to random choice"),
  lambda_rows("mixture_all_dx", "share of errors from random choice"),
  heading("  Mixing weight \u03bb (share of all responses; mean of the participants of each group)"),
  lambda_rows("mixture_all_dx", "lambda, group mean"),
  heading("F. Baseline, pooling and variance parameters, model A"),
  heading("  Intercepts of the systematic component: controls, log-odds relative to R"),
  mix_param_rows("mixture_errors_dx", "alpha", NONREFERENCE_LEVELS, or = TRUE),
  heading("  Pooling across the three contrasts"),
  mix_param_rows("mixture_errors_dx", "mu_beta", paste("\u03bc,", c("Bipolar", "Schizophrenia")), prob = TRUE),
  mix_param_rows("mixture_errors_dx", "tau_beta", paste("\u03c4,", c("Bipolar", "Schizophrenia"))),
  heading("  Random-intercept SDs (log-odds scale)"),
  mix_param_rows("mixture_errors_dx", "sigma_u", paste("Participants,", NONREFERENCE_LEVELS)),
  mix_param_rows("mixture_errors_dx", "sigma_v", paste("Items,", NONREFERENCE_LEVELS)),
  mix_param_rows("mixture_errors_dx", "sigma_w", "Participants, logit(\u03bb)"),
  heading("G. Salience of D in random choice, model Cs (all responses)"),
  heading("  Group contrasts in the systematic component"),
  by_type_rows("mixture_all_dx_salience"),
  heading("  Weight of the D distractor in random choice (relative to each other option)"),
  mix_param_rows("mixture_all_dx_salience", "delta", "\u03b4", or = TRUE, prob = TRUE),
  heading("  Share of errors attributed to random choice"),
  lambda_rows("mixture_all_dx_salience", "share of errors from random choice")
)
names(table_D1)[names(table_D1) == "P_beta_gt_0"] <- "P(\u03b2 > 0)"

# Table D2: model comparison (PSIS-LOO) and observed vs expected counts.
table_D2_loo <- mixture_loo %>%
  transmute(Set = set, Model = mix_names[model],
            `elpd (SE)` = sprintf("%s (%.1f)", fmt_num(elpd_loo), se_elpd_loo),
            `elpd difference (SE)` = sprintf("%s (%.1f)", fmt_num(elpd_diff), se_diff),
            p_loo = sprintf("%.1f", p_loo), `Pareto k > 0.7` = pareto_k_above_0.7)
table_D2_counts <- mixture_counts %>%
  filter(model %in% c("mixture_errors_dx", "mixture_errors_chance"), response_type != "correct") %>%
  mutate(value = sprintf("%.1f [%.1f, %.1f]", expected_median, expected_q2.5, expected_q97.5),
         model = recode(model, mixture_errors_dx = "Expected, model A",
                        mixture_errors_chance = "Expected, model B (chance only)")) %>%
  select(group, response_type, observed, model, value) %>%
  pivot_wider(names_from = model, values_from = value) %>%
  mutate(group = factor(group, MIX_GROUPS), response_type = factor(response_type, OUTCOME_LEVELS)) %>%
  arrange(group, response_type) %>%
  mutate(group = recode(as.character(group), Control = "Controls")) %>%
  rename(Group = group, `Error type` = response_type, Observed = observed)

write_xlsx(
  list(table_D1_appendix = table_D1,
       table_D2_loo     = table_D2_loo,
       table_D2_counts  = table_D2_counts,
       group_contrasts = mixture_contrasts,
       random_share    = mixture_lambda,
       expected_counts = mixture_counts,
       accuracy_part   = mixture_accuracy,
       salience        = mixture_salience,
       loo             = mixture_loo,
       diagnostics     = mixture_diagnostics,
       distractors     = as_tibble(DISTRACTOR_COUNTS) %>%
         mutate(item = seq_len(n()), matriks88_item = PATIENT_ITEMS, matriks12_item = CONTROL_ITEMS, .before = 1)),
  file.path(OUTPUT_DIR, "random_choice_mixture_models.xlsx")
)
