###############################################################################
# 09_distractor_distance.R
# Degree of incorrectness (Section 2.3.3.3; Appendix H):
#   1. Distractor coding: number of departures of each of the 126 distractors
#      (18 items x 7) from the correct response, consensus coding of two raters
#      (data/distractor_distance.csv; coding rules in Table G)
#   2. Descriptives: distance by error category and distractor type; distance
#      of the distractors chosen by each group (matched sample)
#   3. Model 1 at the level of the chosen distractor (matched sample, errors only;
#      stan/multinomial_distance.stan):
#        model1dist_ref         : diagnosis only. Same model as Model 1, with the
#                                 number of distractors of each category on each
#                                 item as an offset
#        model1dist_additive    : + distance of the distractor
#        model1dist_interaction : + diagnosis x distance
#      Sensitivity: additive and interaction models with the number of features
#      changed (0-7) instead of the distance; reference and interaction models
#      without the three items in which repetition errors are recorded with the
#      "r.ic" label.
#      Model 1 (stan/multinomial_exchangeable.stan, as in script 03) is refitted to
#      the same errors, in the same order, for the comparison with the original model.
#   4. Posterior contrasts, PSIS-LOO comparison (chosen distractor: distractor-level
#      models; error category: Model 1 and the distractor-level models), appendix
#      tables (Tables H1-H2)
# Requires outputs/matched_participants.csv (script 03).
# Set MATRIKS_REUSE_FITS=1 to reuse the model1dist_*_fit.rds files of a previous run.
# Outputs (in outputs/): distractor_distance_descriptives.xlsx,
#   distractor_distance_models.xlsx, model1dist_*_fit.rds
###############################################################################

if (!exists("MATRIKS_SETUP_DONE")) source(file.path("R", "00_setup.R"))
require_packages(c(BAYES_PACKAGES, "loo"))
suppressPackageStartupMessages({ library(cmdstanr); library(posterior); library(bayesplot) })

OUTCOME_LEVELS      <- c("R", "D", "IC", "WP")          # R = reference category
NONREFERENCE_LEVELS <- c("D vs R", "IC vs R", "WP vs R")
DISTANCE_FILE       <- file.path("data", "distractor_distance.csv")   # part of the repository
DISTANCE_STAN_FILE  <- file.path("stan", "multinomial_distance.stan")
MODEL1_STAN_FILE    <- file.path("stan", "multinomial_exchangeable.stan")
# The distractor-level models identify an "r.ic" response as the choice of a
# repetition distractor, so they are defined only for RIC_POLICY = "R".
if (RIC_POLICY != "R") stop("Script 09 requires RIC_POLICY = \"R\" (the default).", call. = FALSE)

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

# Long format: one row per participant x item; `item` is the harmonised index
# 1..18 and `label` the response label recorded by MatriKS (lower case).
extract_item_labels <- function(data, item_stub, item_ids, dataset_label) {
  item_cols <- paste0(item_stub, item_ids)
  check_columns(data, item_cols, dataset_label)
  data %>%
    transmute(id        = as.character(external_code),
              diagnosis = harmonize_diagnosis(Diagnosi),
              across(all_of(item_cols), as.character)) %>%
    pivot_longer(all_of(item_cols), names_to = "item_col", values_to = "response_raw") %>%
    mutate(item          = match(as.integer(str_extract(item_col, "\\d+$")), item_ids),
           subject_key   = paste(dataset_label, id, sep = "::"),
           label         = str_to_lower(trimws(response_raw)),
           response_type = classify_response(response_raw, RIC_POLICY)) %>%
    filter(!is.na(response_type), response_type != "UNKNOWN")
}

# ---- 1. Distractor coding -------------------------------------------------------------
check_file(DISTANCE_FILE)
distractors <- read.csv(DISTANCE_FILE, stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  mutate(item = match(matriks88_item, PATIENT_ITEMS),
         distractor = str_to_lower(trimws(distractor)))
stopifnot(!anyNA(distractors$item),
          all(distractors$category %in% OUTCOME_LEVELS),
          all(distractors$matriks12_item == CONTROL_ITEMS[distractors$item]),
          all(distractors$distance == rowSums(distractors[, c("missing", "added", "shape", "size",
                                                              "fill", "line", "orientation")])),
          n_distinct(distractors$item) == length(PATIENT_ITEMS))
distractors <- arrange(distractors, item, distractor)
cat("\nDistractors per item and category:\n")
print(table(CONTROL_ITEMS[distractors$item], distractors$category))

# ---- 2. Matched sample and chosen distractors ---------------------------------------------
controls <- load_controls()
patients <- load_patients()
item_level_all <- bind_rows(
  extract_item_labels(controls, CONTROL_ITEM_STUB, CONTROL_ITEMS, "adult"),
  extract_item_labels(patients, PATIENT_ITEM_STUB, PATIENT_ITEMS, "psy")
) %>%
  filter(diagnosis %in% c("Control", "Bipolar", "Schizophrenia"))

matched_file <- file.path(OUTPUT_DIR, "matched_participants.csv")
if (!file.exists(matched_file)) stop("Run script 03 first: ", matched_file, " not found.", call. = FALSE)
matched_participants <- read.csv(matched_file, stringsAsFactors = FALSE)

error_data <- item_level_all %>%
  semi_join(select(matched_participants, subject_key), by = "subject_key") %>%
  filter(response_type %in% OUTCOME_LEVELS) %>%
  mutate(dx_bipolar       = as.numeric(diagnosis == "Bipolar"),
         dx_schizophrenic = as.numeric(diagnosis == "Schizophrenia"))
cat("\nErrors by group and type (must match Model 1 of script 03):\n")
print(table(error_data$diagnosis, error_data$response_type))

# The "r.ic" label does not say which repetition distractor was chosen: in the
# three items where it occurs, no other repetition label is recorded. A response
# labelled "r.ic" is therefore treated as the choice of one of the item's
# repetition distractors (the likelihood sums over them). Every other label
# identifies exactly one distractor.
RIC_ITEMS <- sort(unique(error_data$item[error_data$label == "r.ic"]))
cat("\nItems with 'r.ic' responses (MatriKS 12+ numbering):", CONTROL_ITEMS[RIC_ITEMS], "\n")

chosen_index <- map(seq_len(nrow(error_data)), function(n) {
  opts <- filter(distractors, item == error_data$item[n])
  idx <- if (error_data$label[n] == "r.ic") which(opts$category == "R") else
    which(opts$distractor == error_data$label[n])
  if (length(idx) == 0) stop("Response label '", error_data$label[n], "' not found among the distractors of item ",
                             CONTROL_ITEMS[error_data$item[n]], " (MatriKS 12+).", call. = FALSE)
  idx
})
error_data$chosen_distance <- map_dbl(seq_len(nrow(error_data)), function(n) {
  mean(filter(distractors, item == error_data$item[n])$distance[chosen_index[[n]]])
})

# ---- 3. Descriptives (how distance and category relate) -------------------------------------
dist_stats <- function(x) tibble(n = length(x), mean = mean(x), sd = sd(x), median = median(x),
                                 min = min(x), max = max(x))
distance_by_category <- distractors %>% group_by(category) %>%
  group_modify(~ dist_stats(.x$distance)) %>% ungroup() %>%
  mutate(category = factor(category, levels = c("IC", "R", "D", "WP"))) %>% arrange(category)
distance_by_type <- distractors %>% group_by(category, distractor) %>%
  group_modify(~ dist_stats(.x$distance)) %>% ungroup()
features_by_category <- distractors %>% group_by(category) %>%
  group_modify(~ dist_stats(.x$features_changed)) %>% ungroup()
# Position of the D distractor within its item: rank among the 7 distractors
# (1 = closest; ties share the mean rank) and whether it is the furthest.
d_position <- distractors %>% group_by(item) %>%
  mutate(rank_in_item = rank(distance), furthest = distance == max(distance)) %>%
  ungroup() %>% filter(category == "D") %>%
  transmute(matriks12_item, distance, rank_in_item, furthest)
chosen_by_group <- error_data %>% group_by(diagnosis) %>%
  group_modify(~ dist_stats(.x$chosen_distance)) %>% ungroup()
chosen_by_group_type <- error_data %>% group_by(diagnosis, response_type) %>%
  summarise(n = n(), mean_distance = mean(chosen_distance), .groups = "drop")
print(distance_by_category); print(chosen_by_group)

write_xlsx(list(distractors           = select(distractors, -item),
                distance_by_category  = mutate(distance_by_category, category = as.character(category)),
                distance_by_type      = distance_by_type,
                features_by_category  = features_by_category,
                d_distractor_position = d_position,
                chosen_by_group       = chosen_by_group,
                chosen_by_group_type  = chosen_by_group_type),
           file.path(OUTPUT_DIR, "distractor_distance_descriptives.xlsx"))

# ---- 4. Models -------------------------------------------------------------------------
PREDICTORS <- c("dx_bipolar", "dx_schizophrenic")

make_distance_stan_data <- function(error_data, chosen_index, distance_var, use_distance, use_interaction) {
  X <- as.matrix(mutate(select(error_data, all_of(PREDICTORS)), across(everything(), as.numeric)))
  I <- length(PATIENT_ITEMS)
  n_opt <- as.integer(table(factor(distractors$item, levels = seq_len(I))))
  D_max <- max(n_opt)
  opt_cat  <- matrix(0L, I, D_max)
  opt_dist <- matrix(0,  I, D_max)
  for (i in seq_len(I)) {
    o <- filter(distractors, item == i)
    opt_cat[i, seq_len(nrow(o))]  <- match(o$category, OUTCOME_LEVELS)
    opt_dist[i, seq_len(nrow(o))] <- o[[distance_var]]
  }
  chosen <- matrix(0L, nrow(error_data), D_max)
  for (n in seq_len(nrow(error_data))) chosen[n, chosen_index[[n]]] <- 1L
  list(N = nrow(error_data), K = length(OUTCOME_LEVELS), P = ncol(X),
       J = n_distinct(error_data$subject_key), I = I, D_max = D_max,
       subject = as.integer(factor(error_data$subject_key)),
       item    = as.integer(error_data$item),
       X = X, n_opt = n_opt, opt_cat = opt_cat, opt_dist = opt_dist, chosen = chosen,
       use_distance = as.integer(use_distance), use_interaction = as.integer(use_interaction))
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

model_specs <- tribble(
  ~label,                               ~distance_var,      ~use_distance, ~use_interaction, ~data,
  "model1dist_ref",                     "distance",         FALSE,         FALSE,            "all",
  "model1dist_additive",                "distance",         TRUE,          FALSE,            "all",
  "model1dist_interaction",             "distance",         TRUE,          TRUE,             "all",
  "model1dist_additive_features",       "features_changed", TRUE,          FALSE,            "all",
  "model1dist_interaction_features",    "features_changed", TRUE,          TRUE,             "all",
  "model1dist_ref_noric",               "distance",         FALSE,         FALSE,            "no_ric",
  "model1dist_interaction_noric",       "distance",         TRUE,          TRUE,             "no_ric"
)

keep_rows <- list(all = rep(TRUE, nrow(error_data)), no_ric = !error_data$item %in% RIC_ITEMS)

REUSE_FITS <- identical(Sys.getenv("MATRIKS_REUSE_FITS"), "1")
distance_model <- cmdstan_model(DISTANCE_STAN_FILE)
fits <- list(); stan_data <- list()
for (i in seq_len(nrow(model_specs))) {
  s <- model_specs[i, ]
  rows <- keep_rows[[s$data]]
  stan_data[[s$label]] <- make_distance_stan_data(error_data[rows, ], chosen_index[rows], s$distance_var,
                                                  s$use_distance, s$use_interaction)
  rds <- file.path(OUTPUT_DIR, paste0(s$label, "_fit.rds"))
  fits[[s$label]] <- if (REUSE_FITS && file.exists(rds)) readRDS(rds) else
    fit_model(distance_model, stan_data[[s$label]], s$label)
}

# Model 1 on the same errors, in the same order: outcome = error category.
stan_data$model1dist_model1 <- with(stan_data$model1dist_ref, list(
  N = N, K = K, P = P, J = J, I = I, subject = subject, item = item, X = X,
  y = match(error_data$response_type, OUTCOME_LEVELS)))
rds <- file.path(OUTPUT_DIR, "model1dist_model1_fit.rds")
fits$model1dist_model1 <- if (REUSE_FITS && file.exists(rds)) readRDS(rds) else
  fit_model(cmdstan_model(MODEL1_STAN_FILE), stan_data$model1dist_model1, "model1dist_model1")
ALL_MODELS <- c("model1dist_model1", model_specs$label)

compute_model_diagnostics <- function(fit, model_label, stan_data) {
  vars <- c("alpha", "mu_beta", "tau_beta", "beta", "sigma_u", "sigma_v")
  if (isTRUE(stan_data$use_distance == 1))    vars <- c(vars, "gamma")
  if (isTRUE(stan_data$use_interaction == 1)) vars <- c(vars, "delta")
  summ <- fit$summary(variables = vars)
  dsum <- fit$diagnostic_summary(quiet = TRUE)
  tibble(model = model_label, observations = stan_data$N, participants = stan_data$J,
         items = stan_data$I,
         max_rhat = max(summ$rhat, na.rm = TRUE),
         min_bulk_ess = min(summ$ess_bulk, na.rm = TRUE),
         min_tail_ess = min(summ$ess_tail, na.rm = TRUE),
         divergences = sum(dsum$num_divergent),
         treedepth_hits = sum(dsum$num_max_treedepth),
         min_ebfmi = min(dsum$ebfmi, na.rm = TRUE))
}
diagnostics <- map_dfr(ALL_MODELS, ~ compute_model_diagnostics(fits[[.x]], .x, stan_data[[.x]]))
print(diagnostics, width = Inf)

# ---- 5. Posterior contrasts -------------------------------------------------------------
draws_mat <- function(fit, variable) {
  m <- as_draws_matrix(fit$draws(variables = variable))
  matrix(as.numeric(m), nrow = nrow(m), dimnames = list(NULL, colnames(m)))
}

get_beta_array <- function(fit, predictor_names = PREDICTORS) {
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

contrast_summary <- function(beta, weights, label) {
  m <- matrix(0, dim(beta)[1], dim(beta)[3], dimnames = list(NULL, NONREFERENCE_LEVELS))
  for (p in names(weights)) m <- m + weights[[p]] * beta[, p, ]
  bind_rows(
    map_dfr(colnames(m), ~ bind_cols(tibble(contrast = label, error_type = .x), summ_draws(m[, .x]))),
    bind_cols(tibble(contrast = label, error_type = "Average of the 3 (vs R)"), summ_draws(rowMeans(m)))
  )
}

type_diff_summary <- function(beta, weights, label) {
  m <- matrix(0, dim(beta)[1], dim(beta)[3])
  for (p in names(weights)) m <- m + weights[[p]] * beta[, p, ]
  short <- sub(" vs R$", "", NONREFERENCE_LEVELS)
  map_dfr(combn(seq_along(short), 2, simplify = FALSE), function(ix) {
    bind_cols(tibble(contrast = label, type_difference = paste(short[ix[1]], "-", short[ix[2]])),
              summ_draws(m[, ix[1]] - m[, ix[2]]))
  })
}

# With distance in the model, the group coefficients compare a distractor of
# category k with an R distractor at the same distance from the correct response.
GROUP_CONTRASTS <- list(
  "BD vs Control" = c(dx_bipolar = 1),
  "SZ vs Control" = c(dx_schizophrenic = 1),
  "SZ vs BD"      = c(dx_schizophrenic = 1, dx_bipolar = -1)
)

contrasts_all <- map_dfr(ALL_MODELS, function(l) {
  beta <- get_beta_array(fits[[l]])
  imap_dfr(GROUP_CONTRASTS, ~ contrast_summary(beta, .x, .y)) %>% mutate(model = l, .before = 1)
})
type_differences <- map_dfr(ALL_MODELS, function(l) {
  beta <- get_beta_array(fits[[l]])
  imap_dfr(GROUP_CONTRASTS, ~ type_diff_summary(beta, .x, .y)) %>% mutate(model = l, .before = 1)
})

# Distance slopes: log-odds of choosing a distractor per one unit of distance,
# at equal error category (OR < 1: distractors further from the correct
# response are chosen less often).
distance_slopes <- map_dfr(model_specs$label[model_specs$use_distance], function(l) {
  g <- as.numeric(draws_mat(fits[[l]], "gamma"))
  out <- if (stan_data[[l]]$use_interaction == 1) {
    dl <- draws_mat(fits[[l]], "delta")
    list("Controls"                              = g,
         "Bipolar"                               = g + dl[, "delta[1]"],
         "Schizophrenia"                         = g + dl[, "delta[2]"],
         "Bipolar - controls (interaction)"      = dl[, "delta[1]"],
         "Schizophrenia - controls (interaction)" = dl[, "delta[2]"],
         "Schizophrenia - bipolar"               = dl[, "delta[2]"] - dl[, "delta[1]"])
  } else {
    list("All groups (common slope)" = g)
  }
  imap_dfr(out, ~ bind_cols(tibble(model = l, distance_var = model_specs$distance_var[model_specs$label == l],
                                   term = .y), summ_draws(.x)))
})
print(distance_slopes, width = Inf)

# ---- 6. PSIS-LOO --------------------------------------------------------------------------
# Chosen distractor: the distractor-level models of each set are fitted to the same
# observations with the same outcome, so their expected log predictive densities
# can be compared; the main set also compares the two distance measures. Model 1
# predicts the error category, not the distractor, and is not in this comparison:
# model1dist_ref is Model 1 on the distractor-level outcome (with the number of
# distractors of each category as an offset), so it is the reference.
# Error category: the outcome of Model 1. The distractor-level models give the
# probability of a category as the sum of the probabilities of its distractors
# (log_lik_cat), so they can be compared with Model 1 on the same 1,099 errors.
loo_list <- map(set_names(model_specs$label), ~ fits[[.x]]$loo(variables = "log_lik", cores = MCMC$parallel_chains))
CATEGORY_MODELS <- c("model1dist_model1", model_specs$label[model_specs$data == "all"])
loo_cat_list <- map(set_names(CATEGORY_MODELS), function(l)
  fits[[l]]$loo(variables = if (l == "model1dist_model1") "log_lik" else "log_lik_cat", cores = MCMC$parallel_chains))
loo_compare_table <- function(labels, lst = loo_list) {
  cmp <- loo::loo_compare(lst[labels])
  cols <- c("elpd_diff", "se_diff", "elpd_loo", "se_elpd_loo", "p_loo", "looic")
  cmp_df <- as.data.frame(unclass(cmp))                 # matrix (loo < 2.9) or data frame
  tibble(model = rownames(cmp_df), !!!lapply(setNames(cols, cols), function(cl) as.numeric(cmp_df[[cl]]))) %>%
    mutate(pareto_k_above_0.7 = map_int(model, ~ sum(lst[[.x]]$diagnostics$pareto_k > 0.7)))
}
distance_loo <- bind_rows(
  mutate(loo_compare_table(c("model1dist_ref", "model1dist_additive", "model1dist_interaction",
                             "model1dist_additive_features", "model1dist_interaction_features")),
         set = "main", .before = 1),
  mutate(loo_compare_table(c("model1dist_ref_noric", "model1dist_interaction_noric")),
         set = "without r.ic items", .before = 1),
  mutate(loo_compare_table(CATEGORY_MODELS, loo_cat_list), set = "error category (with Model 1)", .before = 1)
)
print(distance_loo, width = Inf)

# ---- 7. Appendix tables (Tables H1 and H2) and exports --------------------------------------
fmt2   <- function(x) gsub("-", "\u2212", formatC(ifelse(abs(x) < 0.005, 0, x), format = "f", digits = 2))
fmt1   <- function(x) formatC(x, format = "f", digits = 1)
fmt_ci <- function(lo, hi) sprintf("[%s, %s]", fmt2(lo), fmt2(hi))
fmt_p  <- function(p) sub("^0", "", formatC(p, format = "f", digits = 2))
MODEL_NAMES <- c(
  model1dist_model1               = "Model 1 (error category), refitted to the same errors",
  model1dist_ref                  = "Reference: Model 1 at the distractor level, without distance",
  model1dist_additive             = "Additive: + distance",
  model1dist_interaction          = "Interaction: + distance x diagnosis",
  model1dist_additive_features    = "Additive, features changed: + number of features changed",
  model1dist_interaction_features = "Interaction, features changed: + features changed x diagnosis",
  model1dist_ref_noric            = "Reference, without the r.ic items",
  model1dist_interaction_noric    = "Interaction, without the r.ic items")
CONTRAST_NAMES <- c("BD vs Control" = "Bipolar vs controls", "SZ vs Control" = "Schizophrenia vs controls",
                    "SZ vs BD" = "Schizophrenia vs bipolar")
table_rows <- function(d, panel, group, term) {
  transmute(d, panel = !!panel, group = {{ group }}, term = {{ term }},
            estimate = fmt2(median), cri = fmt_ci(q2.5, q97.5),
            OR = fmt2(OR_median), OR_cri = fmt_ci(OR_q2.5, OR_q97.5), P_gt_0 = fmt_p(P_gt_0))
}
table_H1 <- bind_rows(
  contrasts_all %>% filter(model == "model1dist_interaction") %>%
    table_rows("A. Group contrasts at equal distance, interaction model", CONTRAST_NAMES[contrast], error_type),
  type_differences %>% filter(model == "model1dist_interaction") %>%
    table_rows("B. Differences between error types within each contrast, interaction model",
               CONTRAST_NAMES[contrast], type_difference),
  distance_slopes %>% filter(model != "model1dist_ref_noric") %>%
    table_rows("C. Distance slopes: log-odds of choosing a distractor per unit of distance (or per feature changed)",
               MODEL_NAMES[model], term),
  contrasts_all %>% filter(contrast == "SZ vs Control", error_type %in% c("D vs R", "Average of the 3 (vs R)")) %>%
    table_rows("D. Schizophrenia vs controls in each model", MODEL_NAMES[model], error_type)
)
table_H2 <- distance_loo %>%
  transmute(set, model = MODEL_NAMES[model],
            `elpd (SE)`            = sprintf("%s (%s)", fmt2(elpd_loo), fmt1(se_elpd_loo)),
            `elpd difference (SE)` = sprintf("%s (%s)", fmt2(elpd_diff), fmt1(se_diff)),
            p_loo                  = fmt1(p_loo),
            `Pareto k > 0.7`       = pareto_k_above_0.7)
print(table_H1, n = Inf, width = Inf)
print(table_H2, n = Inf, width = Inf)

write_xlsx(
  list(diagnostics        = diagnostics,
       group_contrasts    = contrasts_all,
       type_differences   = type_differences,
       distance_slopes    = distance_slopes,
       loo                = distance_loo,
       table_H1           = table_H1,
       table_H2           = table_H2,
       ric_items          = tibble(matriks12_item = CONTROL_ITEMS[RIC_ITEMS], matriks88_item = PATIENT_ITEMS[RIC_ITEMS])),
  file.path(OUTPUT_DIR, "distractor_distance_models.xlsx")
)
cat("\nScript 09 done: outputs/distractor_distance_descriptives.xlsx, outputs/distractor_distance_models.xlsx\n")
