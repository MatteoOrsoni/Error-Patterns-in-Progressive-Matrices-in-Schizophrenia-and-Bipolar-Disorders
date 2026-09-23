###############################################################################
# install_dependencies.R
# Installs the R packages and the CmdStan toolchain required by the pipeline.
# Run once from the repository root:  source("install_dependencies.R")
#
# Versions cited in the paper: cmdstanr 0.9.0, posterior 1.7.0, bayesplot 1.15.0.
# The exact versions used for the published results are listed in
# session_info.txt (see README.md).
###############################################################################

cran_packages <- c(
  "readxl", "writexl",                               # data import/export
  "dplyr", "tidyr", "purrr", "tibble", "stringr",    # data manipulation
  "ggplot2", "patchwork",                            # figures
  "emmeans", "broom",                                # ANOVA / regression summaries
  "posterior", "bayesplot"                           # posterior summaries and diagnostics
)

to_install <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

# cmdstanr is distributed through the Stan R-universe (not CRAN).
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
}

# CmdStan (the Stan compiler and sampler). On Windows this requires Rtools
# matching your R version; check_cmdstan_toolchain(fix = TRUE) sets it up.
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstan_ok <- tryCatch({ cmdstanr::cmdstan_version(); TRUE }, error = function(e) FALSE)
if (!cmdstan_ok) cmdstanr::install_cmdstan(cores = 2)

message("CmdStan version: ", cmdstanr::cmdstan_version())
