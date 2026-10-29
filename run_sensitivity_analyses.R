###############################################################################
# run_sensitivity_analyses.R
# Sensitivity analyses of the Bayesian models, each written to its own folder:
#   Appendix B: Models 1-2 without item 34, the item with a different
#               distractor composition            -> outputs_sensitivity_item34/
#   Appendix I: Models 1-2 with "r.ic" responses counted as incomplete correlate
#               (IC) instead of repetition (R) errors -> outputs_sensitivity_ric_ic/
#               (the paper reports Model 1)
# Then R/10_sensitivity_tables.R puts the main and sensitivity estimates side by
# side (outputs/sensitivity_tables.xlsx: Tables B and I).
# Run from the repository root, after run_all.R:  source("run_sensitivity_analyses.R")
###############################################################################

run_sensitivity <- function(settings, output_dir) {
  message("\n==================== Sensitivity analysis: ", output_dir, " ====================")
  do.call(Sys.setenv, c(as.list(settings), MATRIKS_OUTPUT_DIR = output_dir))
  on.exit(Sys.unsetenv(c(names(settings), "MATRIKS_OUTPUT_DIR")), add = TRUE)
  # The settings are read by 00_setup.R, sourced here into a separate environment
  # so that the main analysis settings in the global environment are not changed.
  setup <- new.env(parent = globalenv())
  sys.source(file.path("R", "00_setup.R"), envir = setup)
  for (s in c("R/03_matching_and_models.R", "R/04_posterior_contrasts.R")) {
    source(s, local = new.env(parent = setup))
  }
  invisible(output_dir)
}

run_sensitivity(c(MATRIKS_EXCLUDE_ITEMS = "34"), "outputs_sensitivity_item34")
run_sensitivity(c(MATRIKS_RIC_POLICY = "IC"),    "outputs_sensitivity_ric_ic")

source(file.path("R", "10_sensitivity_tables.R"), local = new.env(parent = globalenv()))
