# Error Patterns in Progressive Matrices in Schizophrenia and Bipolar Disorders — analysis code

Data and R/Stan code to reproduce the analyses of:

> Error Patterns in Progressive Matrices in Schizophrenia and Bipolar Disorders. 

Repository DOI: 10.5281/zenodo.22917478  ·  Licence: MIT

## Contents

```
run_all.R                     runs the whole pipeline
install_dependencies.R        installs R packages, cmdstanr and CmdStan
R/00_setup.R                  shared settings: paths, seed, MCMC settings, sample definitions, helpers
R/01_descriptives.R           Tables 1-2 (demographic and clinical characteristics)
R/02_accuracy_errors.R        ANOVAs and EMMs (Figure 2), error-type proportions (Figure 3)
R/03_matching_and_models.R    performance matching (Figure 1), Bayesian Models 1-2 (Stan)
R/04_posterior_contrasts.R    posterior contrasts (Tables A-B)
R/05_forest_plots.R           Figures 4-5
R/06_symptom_severity.R       symptom-severity regression (Figure 6)
stan/multinomial_exchangeable.stan   Stan program used for both Bayesian models
data/                         input data (see data/README.md)
tools/simulate_example_data.R simulated data with the same structure (for code checking)
tools/make_minimal_dataset.R  keeps only the variables used by the code (for data sharing)
```

## Requirements

- R ≥ 4.3 and the packages listed in `install_dependencies.R`
  (versions cited in the paper: cmdstanr 0.9.0, posterior 1.7.0, bayesplot 1.15.0).
- CmdStan and a C++ toolchain (on Windows: Rtools matching the R version).
- The exact R, package and CmdStan versions used for the published results are recorded in
  `session_info.txt` [to be added by the authors from their final run].

Installation (once, from the repository root):

```r
source("install_dependencies.R")
```

## Running the analyses

1. Put `data_adult.xlsx` and `data_psy.xlsx` in `data/` (see `data/README.md`).
2. Set the working directory to the repository root and run:

```r
source("run_all.R")
```

All results are written to `outputs/`. Each Bayesian model takes several minutes (Model 1: about
6–7 minutes on a 4-core laptop). Scripts can also be run one at a time, in numerical order (scripts 04–05 use files
written by scripts 03–04).


**Paths with non-ASCII characters (e.g. accented letters):** CmdStan's output files are written to a
temporary folder by default. If R reports that the folder contains non-ASCII characters, set an
ASCII-only folder before running, e.g. `Sys.setenv(MATRIKS_STAN_CSV_DIR = "C:/stan_csv")`.

## Outputs and correspondence with the paper

| File in `outputs/` | Paper |
|---|---|
| `Table1_clinical_groups.xlsx`, `Table2_patients_vs_controls.xlsx` | Tables 1–2 |
| `accuracy_errors_anova.xlsx`, `Figure2_accuracy_EMM.jpg` | ANOVAs, EMMs and contrasts; Figure 2 |
| `error_type_proportions.xlsx`, `Figure3_error_type_proportions.jpg` | Error-type proportions; Figure 3 |
| `Figure1_matched_accuracy_density.png`, `data_and_matching_outputs.xlsx` | Matching; Figure 1; Model 2 VIFs |
| `multinomial_model_results.xlsx`, `diagnostics/` | Model 1–2 coefficients and sampling diagnostics |
| `model_contrasts.xlsx` | Posterior contrasts (Tables A–B) |
| `Figure4_model1_forest.png`, `Figure5_model2_forest.png` | Figures 4–5 |
| `symptom_severity_regression.xlsx`, `Figure6_severity_regression.jpg` | Symptom-severity regression; Figure 6 |
| `model1_fit.rds`, `model2_fit.rds` | Complete posterior draws (cmdstanr objects) |
| `session_info.txt` | R, package and CmdStan versions |

## Implementation notes

- **Sample.** Controls: participants aged ≥ 18 with no recorded diagnosis (rows with gender `B` or
  missing gender are excluded), n = 371. Patients: n = 59 after excluding BSPDC40, BSPDC41 and
  BSPDC45 [reason to be documented by the authors].
- **Items.** All analyses use the 18 items shared by the two MatriKS versions; the correspondence
  between item numbers is defined in `R/00_setup.R`.
- **Accuracy.** In the ANOVAs, accuracy is the number of correct responses on the 18 items (skipped
  items count as not correct). For matching, accuracy is the proportion of correct responses among
  answered (non-skipped) items.
- **Matching.** Greedy nearest-neighbour matching in triplets (one control, one bipolar and one
  schizophrenia participant), without replacement and without caliper, anchored on the smallest
  group (schizophrenia, n = 28): Model 1 therefore includes 28 participants per group. The
  procedure is deterministic.
- **`r.ic` responses.** Classified as repetition errors in all error-type analyses (Bayesian models
  and Figure 3), as set by `RIC_POLICY = "R"` in `R/00_setup.R`; this reproduces the published
  results. [Rationale to be documented by the authors.]
- **Error-type percentages (Figure 3).** Descriptive, full (unmatched) sample: errors pooled within
  each group; 95% CIs by participant-level bootstrap (5,000 resamples).
- **Model 2 covariates.** Attention and Digit Span scores are z-standardised across all patients
  with a valid score on that measure, before excluding patients with missing values on any
  covariate.
- **Symptom-severity regression.** BPRS and HoNOS are mean-centred on the complete-case sample; the
  same centring constants are used in the sensitivity analysis.

## Reproducibility

- Random seeds: `SEED = 42` for Stan (all chains) and for the bootstrap CIs of Figure 3;
  `set.seed(1)` for the simulation of model-implied proportions (`R/04_posterior_contrasts.R`).
  All other steps are deterministic.
- Posterior draws are exactly reproducible only with the same CmdStan version, compiler and
  operating system; otherwise small Monte Carlo differences (typically in the second decimal)
  are expected.
- Generative AI (Claude, Anthropic) was used to assist in writing and debugging the code; all code
  and outputs were checked by the authors (see the paper's Methods).

## Contact

[Corresponding author, e-mail]
