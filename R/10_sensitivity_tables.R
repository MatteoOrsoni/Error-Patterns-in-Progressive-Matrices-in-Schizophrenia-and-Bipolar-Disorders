###############################################################################
# 10_sensitivity_tables.R
# Main analysis and sensitivity analyses side by side (odds ratios [95% CrI] and
# posterior probability of a positive coefficient):
#   Table B (Appendix B): Models 1-2 without item 34, the item with a different
#            distractor composition            (folder outputs_sensitivity_item34/)
#   Table I (Appendix I): Model 1 with "r.ic" responses counted as incomplete
#            correlate (IC) instead of repetition (R) errors
#                                              (folder outputs_sensitivity_ric_ic/)
# Input : model_contrasts.xlsx, multinomial_model_results.xlsx and
#         model_stan_data.rds in outputs/ (run_all.R) and in the two sensitivity
#         folders (run_sensitivity_analyses.R)
# Output: outputs/sensitivity_tables.xlsx
###############################################################################

if (!exists("MATRIKS_SETUP_DONE")) source(file.path("R", "00_setup.R"))
suppressPackageStartupMessages(library(readxl))

MAIN_DIR <- Sys.getenv("MATRIKS_OUTPUT_DIR", unset = "outputs")
SENSITIVITY_DIRS <- c(table_B = "outputs_sensitivity_item34", table_I = "outputs_sensitivity_ric_ic")
OUTCOME_LEVELS <- c("R", "D", "IC", "WP")

fmt_ci  <- function(m, lo, hi) sprintf("%.2f [%.2f, %.2f]", m, lo, hi)
fmt_pos <- function(p) ifelse(p > .99, ">.99", ifelse(p < .01, "<.01", sub("^0", "", sprintf("%.2f", p))))

# One row per reported quantity of one analysis (folder).
collect_estimates <- function(dir, include_model2) {
  xl <- function(file, sheet) read_xlsx(file.path(dir, file), sheet = sheet)
  contrast_rows <- function(d, panel) {
    transmute(d, panel = panel, parameter = paste0(contrast, ", ", error_type),
              estimate = fmt_ci(OR_median, OR_q2.5, OR_q97.5), p_positive = fmt_pos(P_gt_0))
  }
  sd_rows <- function(d, label) {
    transmute(d, panel = "Model 1: intercepts and random-intercept SDs (log-odds)",
              parameter = paste0(label, ", ", comparison),
              estimate = fmt_ci(median, q2.5, q97.5), p_positive = "")
  }
  stan_data <- readRDS(file.path(dir, "model_stan_data.rds"))$model1
  group <- ifelse(stan_data$X[, "dx_bipolar"] == 1, "Bipolar",
                  ifelse(stan_data$X[, "dx_schizophrenic"] == 1, "Schizophrenia", "Control"))
  counts <- as.data.frame(table(group = group, type = factor(OUTCOME_LEVELS[stan_data$y], OUTCOME_LEVELS)))
  out <- bind_rows(
    transmute(counts, panel = "Model 1: errors by type (n)", parameter = paste0(group, ", ", type),
              estimate = as.character(Freq), p_positive = ""),
    contrast_rows(xl("model_contrasts.xlsx", "model1_group_contrasts"),
                  "Model 1: group contrasts (OR vs repetition errors)"),
    if (include_model2) {
      xl("model_contrasts.xlsx", "model2_slopes_contrasts") %>%
        filter(contrast == "SZ vs BD at mean covariates" | grepl("(interaction)", contrast, fixed = TRUE)) %>%
        contrast_rows("Model 2: diagnosis and diagnosis x cognitive predictors (OR vs repetition errors)")
    },
    transmute(xl("multinomial_model_results.xlsx", "model1_intercepts"),
              panel = "Model 1: intercepts and random-intercept SDs (log-odds)",
              parameter = paste0("Intercept (controls), ", comparison),
              estimate = fmt_ci(median, q2.5, q97.5), p_positive = ""),
    sd_rows(xl("multinomial_model_results.xlsx", "model1_random_sd"), "SD participants"),
    sd_rows(xl("multinomial_model_results.xlsx", "model1_item_sd"),   "SD items")
  )
  as_tibble(out)
}

side_by_side <- function(sens_dir, label, include_model2) {
  main <- collect_estimates(MAIN_DIR, include_model2)
  sens <- collect_estimates(sens_dir, include_model2)
  full_join(rename(main, main_estimate = estimate, main_p_positive = p_positive),
            rename(sens, sensitivity_estimate = estimate, sensitivity_p_positive = p_positive),
            by = c("panel", "parameter")) %>%
    rename("Main analysis: OR or estimate [95% CrI]" = main_estimate,
           "Main analysis: P(beta > 0)" = main_p_positive,
           !!paste0(label, ": OR or estimate [95% CrI]") := sensitivity_estimate,
           !!paste0(label, ": P(beta > 0)") := sensitivity_p_positive)
}

tables <- list()
spec <- list(table_B = list(label = "Without item 34", model2 = TRUE),
             table_I = list(label = "r.ic counted as IC", model2 = FALSE))
for (nm in names(SENSITIVITY_DIRS)) {
  d <- SENSITIVITY_DIRS[[nm]]
  if (!file.exists(file.path(d, "model_contrasts.xlsx"))) {
    message("Skipped ", nm, ": run run_sensitivity_analyses.R first (no results in ", d, "/).")
    next
  }
  tables[[nm]] <- side_by_side(d, spec[[nm]]$label, spec[[nm]]$model2)
  print(tables[[nm]], n = Inf, width = Inf)
}
if (length(tables) > 0) write_xlsx(tables, file.path(MAIN_DIR, "sensitivity_tables.xlsx"))
