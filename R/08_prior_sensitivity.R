###############################################################################
# 08_prior_sensitivity.R
# Prior sensitivity of the partial-pooling priors (Appendix F).
# Models 1 and 2 are refitted on the same data with alternative priors on the
# mean (mu_beta) and spread (tau_beta) of each predictor's coefficients across
# the three contrasts:
#   main                 : mu_beta ~ Normal(0, 0.5), tau_beta ~ Exponential(1)
#                          (paper; draws read from model1_fit.rds / model2_fit.rds)
#   mu_sd_1              : mu_beta ~ Normal(0, 1),   tau_beta ~ Exponential(1)
#   tau_rate_0.5         : mu_beta ~ Normal(0, 0.5), tau_beta ~ Exponential(0.5)
#   mu_sd_1_tau_rate_0.5 : mu_beta ~ Normal(0, 1),   tau_beta ~ Exponential(0.5)
# All other priors and the MCMC settings are those of the main analyses.
# For each setting, the contrasts reported in the paper are recomputed on the
# posterior draws as in script 04.
# Input : outputs/model_stan_data.rds, outputs/model1_fit.rds, outputs/model2_fit.rds
#         (script 03)
# Output: outputs/prior_sensitivity.xlsx (Table F)
###############################################################################

if (!exists("MATRIKS_SETUP_DONE")) source(file.path("R", "00_setup.R"))
require_packages(BAYES_PACKAGES)
suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })

NONREFERENCE_LEVELS <- c("D vs R", "IC vs R", "WP vs R")
MODEL1_PREDICTORS   <- c("dx_bipolar", "dx_schizophrenic")
MODEL2_PREDICTORS   <- c("dx_schizophrenic", "matrici_attentive_z", "span_forward_z", "span_backward_z",
                         "dx_x_attentive", "dx_x_span_forward", "dx_x_span_backward")
SENS_STAN_FILE      <- file.path("stan", "multinomial_exchangeable_prior_sensitivity.stan")

PRIOR_SETTINGS <- tibble(
  setting        = c("main", "mu_sd_1", "tau_rate_0.5", "mu_sd_1_tau_rate_0.5"),
  prior_mu_sd    = c(0.5, 1, 0.5, 1),
  prior_tau_rate = c(1, 1, 0.5, 0.5)
)

stan_data <- readRDS(file.path(OUTPUT_DIR, "model_stan_data.rds"))
main_fits <- list(model1 = readRDS(file.path(OUTPUT_DIR, "model1_fit.rds")),
                  model2 = readRDS(file.path(OUTPUT_DIR, "model2_fit.rds")))

# ---- Helpers (as in script 04) ---------------------------------------------------
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
         P_gt_0 = mean(x > 0), cri95_excludes_0 = qs[1] > 0 | qs[2] < 0)
}

# Each contrast is a draws x 3 matrix (D, IC, WP vs R); the average of the three
# is added as in script 04.
contrast_by_type <- function(mats) {
  imap_dfr(mats, function(m, lab) {
    bind_rows(
      map_dfr(colnames(m), ~ bind_cols(tibble(contrast = lab, error_type = .x), summ_draws(m[, .x]))),
      bind_cols(tibble(contrast = lab, error_type = "Average of the 3 (vs R)"), summ_draws(rowMeans(m)))
    )
  })
}

model1_contrasts <- function(fit) {
  b  <- get_beta_array(fit, MODEL1_PREDICTORS)
  bd <- b[, "dx_bipolar", ]; sz <- b[, "dx_schizophrenic", ]
  contrast_by_type(list("BD vs Control" = bd, "SZ vs Control" = sz, "SZ vs BD" = sz - bd))
}

model2_contrasts <- function(fit) {
  b <- get_beta_array(fit, MODEL2_PREDICTORS)
  cov_map <- c("attention" = "matrici_attentive_z", "forward span" = "span_forward_z",
               "backward span" = "span_backward_z")
  int_map <- c("attention" = "dx_x_attentive", "forward span" = "dx_x_span_forward",
               "backward span" = "dx_x_span_backward")
  slopes <- list("SZ vs BD at mean covariates" = b[, "dx_schizophrenic", ])
  for (nm in names(cov_map)) {
    slopes[[paste0(nm, ": slope in BD")]]           <- b[, cov_map[[nm]], ]
    slopes[[paste0(nm, ": slope in SZ")]]           <- b[, cov_map[[nm]], ] + b[, int_map[[nm]], ]
    slopes[[paste0(nm, ": SZ - BD (interaction)")]] <- b[, int_map[[nm]], ]
  }
  contrast_by_type(slopes)
}

pooling_summary <- function(fit, predictor_names) {
  bind_rows(
    map_dfr(seq_along(predictor_names), function(p) {
      bind_rows(bind_cols(tibble(term = predictor_names[p], parameter = "mu_beta"),
                          summ_draws(draws_mat(fit, sprintf("mu_beta[%d]", p))[, 1])),
                bind_cols(tibble(term = predictor_names[p], parameter = "tau_beta"),
                          summ_draws(draws_mat(fit, sprintf("tau_beta[%d]", p))[, 1])))
    })
  ) %>% select(term, parameter, median, q2.5, q97.5)
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

fit_diagnostics <- function(fit) {
  summ <- fit$summary(variables = c("alpha", "mu_beta", "tau_beta", "beta", "sigma_u", "sigma_v"))
  dsum <- fit$diagnostic_summary(quiet = TRUE)
  tibble(max_rhat = max(summ$rhat, na.rm = TRUE),
         min_bulk_ess = min(summ$ess_bulk, na.rm = TRUE),
         min_tail_ess = min(summ$ess_tail, na.rm = TRUE),
         divergences = sum(dsum$num_divergent),
         treedepth_hits = sum(dsum$num_max_treedepth),
         min_ebfmi = min(dsum$ebfmi, na.rm = TRUE))
}

# ---- Fits --------------------------------------------------------------------------
# The fitted objects are not saved (about 200-270 MB each); only their summaries.
sens_model <- cmdstan_model(SENS_STAN_FILE)

results <- pmap(PRIOR_SETTINGS, function(setting, prior_mu_sd, prior_tau_rate) {
  map(c(model1 = "model1", model2 = "model2"), function(m) {
    if (setting == "main") {
      fit <- main_fits[[m]]
    } else {
      message("\n---- ", m, ", prior setting ", setting, " ----")
      fit <- sens_model$sample(
        data = c(stan_data[[m]], list(prior_mu_sd = prior_mu_sd, prior_tau_rate = prior_tau_rate)),
        seed = SEED, chains = MCMC$chains, parallel_chains = MCMC$parallel_chains,
        iter_warmup = MCMC$iter_warmup, iter_sampling = MCMC$iter_sampling,
        adapt_delta = MCMC$adapt_delta, max_treedepth = MCMC$max_treedepth,
        refresh = 0, output_dir = check_stan_csv_dir(STAN_CSV_DIR),
        output_basename = paste0(m, "_prior_", setting)
      )
    }
    predictors <- if (m == "model1") MODEL1_PREDICTORS else MODEL2_PREDICTORS
    list(contrasts   = if (m == "model1") model1_contrasts(fit) else model2_contrasts(fit),
         pooling     = pooling_summary(fit, predictors),
         diagnostics = fit_diagnostics(fit))
  })
})
names(results) <- PRIOR_SETTINGS$setting

collect <- function(m, part) {
  imap_dfr(results, ~ mutate(.x[[m]][[part]], setting = .y, .before = 1)) %>%
    left_join(PRIOR_SETTINGS, by = "setting") %>%
    relocate(prior_mu_sd, prior_tau_rate, .after = setting)
}

# Compact comparison: "OR [95% CrI], P(beta > 0)" for each setting side by side
fmt_or <- function(d) {
  d %>%
    mutate(value = sprintf("%.2f [%.2f, %.2f], %.2f", OR_median, OR_q2.5, OR_q97.5, P_gt_0)) %>%
    select(contrast, error_type, setting, value) %>%
    pivot_wider(names_from = setting, values_from = value)
}

model1_long <- collect("model1", "contrasts")
model2_long <- collect("model2", "contrasts")
diagnostics <- bind_rows(collect("model1", "diagnostics") %>% mutate(model = "model1", .before = 1),
                         collect("model2", "diagnostics") %>% mutate(model = "model2", .before = 1))
print(fmt_or(model1_long), n = Inf, width = Inf)
print(fmt_or(model2_long), n = Inf, width = Inf)
print(diagnostics, width = Inf)

write_xlsx(
  list(settings          = PRIOR_SETTINGS,
       model1_comparison = fmt_or(model1_long),
       model2_comparison = fmt_or(model2_long),
       model1_contrasts  = model1_long,
       model2_contrasts  = model2_long,
       pooling_mu_tau    = bind_rows(collect("model1", "pooling") %>% mutate(model = "model1", .before = 1),
                                     collect("model2", "pooling") %>% mutate(model = "model2", .before = 1)),
       diagnostics       = diagnostics),
  file.path(OUTPUT_DIR, "prior_sensitivity.xlsx")
)
