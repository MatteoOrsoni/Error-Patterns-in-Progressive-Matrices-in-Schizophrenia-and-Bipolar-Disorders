###############################################################################
# 04_posterior_contrasts.R
# Contrasts computed on each posterior draw (preserving the posterior
# correlation between coefficients; see the paper's Methods):
#   Model 1: BD vs Control, SZ vs Control and the direct SZ - BD contrast for
#            each error type (D, IC, WP vs R) and averaged across the three;
#            differences between error types (D - IC, D - WP, IC - WP);
#            model-implied proportions of each error type by group.
#   Model 2: SZ - BD at mean covariate levels; slope of each cognitive
#            predictor in BD (main effect), in SZ (main effect + interaction)
#            and their difference (interaction).
# Input : outputs/model1_fit.rds, outputs/model2_fit.rds (from script 03)
# Output: outputs/model_contrasts.xlsx (paper: Tables A and B, Figures 4-5)
###############################################################################

if (!exists("MATRIKS_SETUP_DONE")) source(file.path("R", "00_setup.R"))
require_packages(BAYES_PACKAGES)
suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })

OUTCOME_LEVELS      <- c("R", "D", "IC", "WP")
NONREFERENCE_LEVELS <- c("D vs R", "IC vs R", "WP vs R")
MODEL1_PREDICTORS   <- c("dx_bipolar", "dx_schizophrenic")
MODEL2_PREDICTORS   <- c("dx_schizophrenic", "matrici_attentive_z", "span_forward_z", "span_backward_z",
                         "dx_x_attentive", "dx_x_span_forward", "dx_x_span_backward")

model1_fit <- readRDS(file.path(OUTPUT_DIR, "model1_fit.rds"))
model2_fit <- readRDS(file.path(OUTPUT_DIR, "model2_fit.rds"))

# ---- Helpers ------------------------------------------------------------------
draws_mat <- function(fit, variable) {
  m <- as_draws_matrix(fit$draws(variables = variable))
  matrix(as.numeric(m), nrow = nrow(m), dimnames = list(NULL, colnames(m)))
}

# beta as array [draw, predictor, contrast]
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

# For each contrast (a draws x 3 matrix, columns D/IC/WP vs R):
#   by_type  : effect on each error type and the average of the three;
#   type_diff: differences between error types (contrast x error-type interaction).
contrast_tables <- function(mats) {
  by_type <- imap_dfr(mats, function(m, lab) {
    bind_rows(
      map_dfr(colnames(m), ~ bind_cols(tibble(contrast = lab, error_type = .x), summ_draws(m[, .x]))),
      bind_cols(tibble(contrast = lab, error_type = "Average of the 3 (vs R)"), summ_draws(rowMeans(m)))
    )
  })
  short <- sub(" vs R$", "", colnames(mats[[1]]))
  type_diff <- imap_dfr(mats, function(m, lab) {
    map_dfr(combn(seq_along(short), 2, simplify = FALSE), function(ix) {
      bind_cols(tibble(contrast = lab, type_difference = paste(short[ix[1]], "-", short[ix[2]])),
                summ_draws(m[, ix[1]] - m[, ix[2]]))
    })
  })
  list(by_type = by_type, type_diff = type_diff)
}

# ---- Model 1: group contrasts --------------------------------------------------
b1 <- get_beta_array(model1_fit, MODEL1_PREDICTORS)
bd <- b1[, "dx_bipolar", ]
sz <- b1[, "dx_schizophrenic", ]
m1 <- contrast_tables(list("BD vs Control" = bd, "SZ vs Control" = sz, "SZ vs BD" = sz - bd))

print(select(m1$by_type, contrast, error_type, median, q2.5, q97.5, OR_median, OR_q2.5, OR_q97.5, P_gt_0),
      n = Inf, width = Inf)
print(select(m1$type_diff, contrast, type_difference, median, q2.5, q97.5, P_gt_0), n = Inf, width = Inf)

# ---- Model 1: model-implied proportions of each error type by group --------------
# "typical" : average participant and item (random effects = 0)
# "marginal": averaged over participants and items (M simulated random-effect
#             draws per posterior draw from the posterior random-effect SDs)
model_probs <- function(fit, predictor_names, x_list, M = 200, seed = 1) {
  set.seed(seed)
  alpha <- draws_mat(fit, "alpha"); su <- draws_mat(fit, "sigma_u"); sv <- draws_mat(fit, "sigma_v")
  beta  <- get_beta_array(fit, predictor_names)
  S <- nrow(alpha); K1 <- ncol(alpha)
  map(x_list, function(x) {
    lin <- matrix(NA_real_, S, K1)
    for (k in seq_len(K1)) lin[, k] <- alpha[, k] + as.numeric(matrix(beta[, , k], nrow = S) %*% x)
    e_typ <- cbind(1, exp(lin)); p_typ <- e_typ / rowSums(e_typ)
    p_mar <- matrix(0, S, K1 + 1)
    for (m in seq_len(M)) {
      eta <- lin + matrix(rnorm(S * K1), S, K1) * su + matrix(rnorm(S * K1), S, K1) * sv
      e   <- cbind(1, exp(eta)); p_mar <- p_mar + e / rowSums(e)
    }
    p_mar <- p_mar / M
    colnames(p_typ) <- colnames(p_mar) <- OUTCOME_LEVELS
    list(typical = p_typ, marginal = p_mar)
  })
}

summ_pct <- function(x) {
  qs <- unname(quantile(x, c(0.025, 0.975)))
  tibble(median = median(x), q2.5 = qs[1], q97.5 = qs[2], P_gt_0 = mean(x > 0))
}

x1 <- list(Control = c(0, 0), BD = c(1, 0), SZ = c(0, 1))       # order = MODEL1_PREDICTORS
probs1 <- model_probs(model1_fit, MODEL1_PREDICTORS, x1)

prob_table <- map_dfr(c("typical", "marginal"), function(sc) {
  map_dfr(names(x1), function(g) {
    map_dfr(OUTCOME_LEVELS, function(o) {
      bind_cols(tibble(scale = sc, group = g, error_type = o),
                select(summ_pct(100 * probs1[[g]][[sc]][, o]), -P_gt_0))
    })
  })
})

diff_pairs <- list("BD - Control" = c("BD", "Control"), "SZ - Control" = c("SZ", "Control"),
                   "SZ - BD" = c("SZ", "BD"))
prob_diff <- map_dfr(c("typical", "marginal"), function(sc) {
  imap_dfr(diff_pairs, function(gg, lab) {
    map_dfr(OUTCOME_LEVELS, function(o) {
      bind_cols(tibble(scale = sc, contrast = lab, error_type = o),
                summ_pct(100 * (probs1[[gg[1]]][[sc]][, o] - probs1[[gg[2]]][[sc]][, o])))
    })
  })
})
print(filter(prob_table, scale == "marginal"), n = Inf, width = Inf)
print(filter(prob_diff,  scale == "marginal"), n = Inf, width = Inf)

# ---- Model 2: diagnosis effect and simple slopes ---------------------------------
b2 <- get_beta_array(model2_fit, MODEL2_PREDICTORS)
cov_map <- c("attention" = "matrici_attentive_z", "forward span" = "span_forward_z",
             "backward span" = "span_backward_z")
int_map <- c("attention" = "dx_x_attentive", "forward span" = "dx_x_span_forward",
             "backward span" = "dx_x_span_backward")

slopes <- list("SZ vs BD at mean covariates" = b2[, "dx_schizophrenic", ])
for (nm in names(cov_map)) {
  slopes[[paste0(nm, ": slope in BD")]]           <- b2[, cov_map[[nm]], ]
  slopes[[paste0(nm, ": slope in SZ")]]           <- b2[, cov_map[[nm]], ] + b2[, int_map[[nm]], ]
  slopes[[paste0(nm, ": SZ - BD (interaction)")]] <- b2[, int_map[[nm]], ]
}
m2 <- contrast_tables(slopes)
print(select(m2$by_type, contrast, error_type, median, q2.5, q97.5, OR_median, OR_q2.5, OR_q97.5, P_gt_0),
      n = Inf, width = Inf)

# ---- Export --------------------------------------------------------------------
write_xlsx(
  list(model1_group_contrasts   = m1$by_type,
       model1_type_dissociation = m1$type_diff,
       model1_probabilities     = prob_table,
       model1_prob_differences  = prob_diff,
       model2_slopes_contrasts  = m2$by_type,
       model2_type_dissociation = m2$type_diff),
  file.path(OUTPUT_DIR, "model_contrasts.xlsx")
)
