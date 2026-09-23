###############################################################################
# run_all.R
# Runs the complete analysis pipeline of
#   "Error Patterns in Progressive Matrices in Schizophrenia and Bipolar Disorders"
# Run from the repository root:  source("run_all.R")
# Each script is run in its own environment, so it only depends on the shared
# setup (R/00_setup.R) and on files written by earlier scripts.
###############################################################################

source(file.path("R", "00_setup.R"))

scripts <- c(
  "R/01_descriptives.R",           # Tables 1-2
  "R/02_accuracy_errors.R",        # ANOVAs, EMMs, Figures 2-3
  "R/03_matching_and_models.R",    # matching (Figure 1), Models 1-2 (Stan)
  "R/04_posterior_contrasts.R",    # posterior contrasts (Tables A-B)
  "R/05_forest_plots.R",           # Figures 4-5
  "R/06_symptom_severity.R"        # symptom-severity regression, Figure 6
)

for (s in scripts) {
  message("\n==================== ", s, " ====================")
  source(s, local = new.env(parent = globalenv()))
}

write_session_info()
message("\nDone. All outputs are in '", OUTPUT_DIR, "/' (R and package versions in session_info.txt).")
