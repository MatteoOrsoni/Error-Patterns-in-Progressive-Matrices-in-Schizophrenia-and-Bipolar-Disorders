###############################################################################
# run_all.R
# Runs the complete analysis pipeline of
#   "Error Patterns in Progressive Matrices in Schizophrenia and Bipolar Disorder"
# Run from the repository root:  source("run_all.R")
# Each script is run in its own environment, so it only depends on the shared
# setup (R/00_setup.R) and on files written by earlier scripts.
# The sensitivity analyses of Appendices B and I are run by
# run_sensitivity_analyses.R, after this script.
###############################################################################

source(file.path("R", "00_setup.R"))

scripts <- c(
  "R/01_descriptives.R",            # Tables 1-2
  "R/02_accuracy_errors.R",         # ANOVAs and EMMs (Figure 2), error-type percentages (Figure 3)
  "R/03_matching_and_models.R",     # matching (Table 4), Models 1-2 (Stan), Table E
  "R/04_posterior_contrasts.R",     # posterior contrasts (Appendix A, Tables A1-A2)
  "R/05_forest_plots.R",            # Figures 4-5
  "R/06_symptom_severity.R",        # symptom-severity regression (Figure 6)
  "R/07_unsystematic_responding.R", # response times and random-choice mixture models (Appendices C-D)
  "R/08_prior_sensitivity.R",       # Models 1-2 under alternative pooling priors (Appendix F)
  "R/09_distractor_distance.R"      # distance of the distractors from the correct response (Appendix H)
)

for (s in scripts) {
  message("\n==================== ", s, " ====================")
  source(s, local = new.env(parent = globalenv()))
}

write_session_info()
message("\nDone. All outputs are in '", OUTPUT_DIR, "/' (R and package versions in session_info.txt).")
