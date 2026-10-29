# Error Patterns in Progressive Matrices in Schizophrenia and Bipolar Disorder: analysis code

R and Stan code to reproduce the analyses of:

> [Authors] ([year]). Error Patterns in Progressive Matrices in Schizophrenia and Bipolar Disorder.
> *Intelligence*. [DOI]

Repository DOI: [Zenodo DOI] · Licence: see [Licence](#licence)

The study compares accuracy and error types (repetition, difference, incomplete correlate and wrong
principle errors) on the 18 items shared by two versions of the MatriKS matrix reasoning test in
patients with schizophrenia (n = 28), patients with bipolar disorder (n = 31) and healthy controls
(n = 371). The error-type analyses are Bayesian multinomial mixed-effects models fitted with Stan.

**The study data are not included in this repository** (see [`data/README.md`](data/README.md) for
data availability and for the structure of the two input files). To check that the code runs, the
repository includes a script that writes simulated data with the same structure (see
[Running on simulated data](#running-on-simulated-data)).

## Contents

```
run_all.R                      runs the main pipeline (scripts 01-09)
run_sensitivity_analyses.R     sensitivity analyses of Appendices B and I (after run_all.R)
install_dependencies.R         installs the R packages, cmdstanr and CmdStan
R/00_setup.R                   shared settings: paths, seed, MCMC settings, sample definitions, helpers
R/01_descriptives.R            Tables 1-2 (demographic and clinical characteristics)
R/02_accuracy_errors.R         ANOVAs and EMMs (Figure 2), error-type percentages (Figure 3)
R/03_matching_and_models.R     performance matching (Table 4), Bayesian Models 1-2, Table E
R/04_posterior_contrasts.R     posterior contrasts (Appendix A, Tables A1-A2)
R/05_forest_plots.R            Figures 4-5
R/06_symptom_severity.R        symptom-severity regression (Figure 6)
R/07_unsystematic_responding.R response times and random-choice mixture models (Appendices C-D)
R/08_prior_sensitivity.R       Models 1-2 under alternative pooling priors (Appendix F)
R/09_distractor_distance.R     distance of the distractors from the correct response (Appendix H)
R/10_sensitivity_tables.R      main and sensitivity estimates side by side (Tables B and I)
stan/multinomial_exchangeable.stan                    Models 1-2 and the response-time models
stan/multinomial_exchangeable_prior_sensitivity.stan  same model, pooling-prior scales as data (script 08)
stan/mixture_random_choice.stan                       random-choice mixture models (script 07)
stan/multinomial_distance.stan                        Model 1 at the level of the chosen distractor (script 09)
data/README.md                 input files: data availability, columns, response labels
data/distractor_distance.csv   coding of the 126 distractors (script 09; Appendices G-H)
tools/simulate_example_data.R  simulated data with the structure of the study files
tools/make_minimal_dataset.R   keeps only the columns used by the code (to prepare the data for sharing)
```

## Requirements

- R ≥ 4.3 and the packages listed in `install_dependencies.R` (versions cited in the paper:
  cmdstanr 0.9.0, posterior 1.7.0, bayesplot 1.15.0).
- CmdStan and a C++ toolchain (on Windows, Rtools matching the R version).
- The R, package and CmdStan versions of the run that produced the published results are listed in
  `session_info.txt` [to be added by the authors from their final run]. Every run writes this file to
  the output folder.

Installation (once, from the repository root):

```r
source("install_dependencies.R")
```

## Running the analyses

1. Put `data_adult.xlsx` and `data_psy.xlsx` in `data/` (see `data/README.md`).
2. Set the working directory to the repository root (the folder that contains `run_all.R`) and run:

```r
source("run_all.R")                    # main analyses, results in outputs/
source("run_sensitivity_analyses.R")   # Appendices B and I (optional, after run_all.R)
```

Each Bayesian model takes several minutes (Model 1: about 6–7 minutes on a 4-core laptop). Script 07
fits several models and is the longest step. The scripts can also be run one at a time, in numerical
order: scripts 04–05 read files written by scripts 03–04; scripts 07 and 09 read
`outputs/matched_participants.csv` (script 03); script 08 reads `outputs/model_stan_data.rds` and the
Model 1–2 fits (script 03).

### Sensitivity analyses (Appendices B and I)

`run_sensitivity_analyses.R` refits Models 1–2 (scripts 03–04) twice, each time in its own folder,
and then puts the estimates next to those of the main analysis (`R/10_sensitivity_tables.R`):

| Analysis | Setting | Folder | Paper |
|---|---|---|---|
| Without item 34, the only item with a different distractor composition | `MATRIKS_EXCLUDE_ITEMS = "34"` | `outputs_sensitivity_item34/` | Appendix B, Table B |
| `r.ic` responses counted as incomplete correlate (IC) instead of repetition (R) errors | `MATRIKS_RIC_POLICY = "IC"` | `outputs_sensitivity_ric_ic/` | Appendix I, Table I (Model 1) |

Both tables are written to `outputs/sensitivity_tables.xlsx` (sheets `table_B` and `table_I`). The
main analysis in `outputs/` is not modified.

### Settings (environment variables)

All settings have defaults that reproduce the paper; they can be changed with `Sys.setenv()` before
running the scripts.

| Variable | Default | Effect |
|---|---|---|
| `MATRIKS_DATA_DIR` | `data` | Folder with `data_adult.xlsx` and `data_psy.xlsx` |
| `MATRIKS_OUTPUT_DIR` | `outputs` | Folder for all results |
| `MATRIKS_STAN_CSV_DIR` | a temporary folder | Folder for CmdStan's CSV files; must be an ASCII-only path (see below) |
| `MATRIKS_QUICK_TEST` | unset | `1` = 2 chains × 150 iterations, to check that the code runs (results do **not** match the paper) |
| `MATRIKS_REUSE_FITS` | unset | `1` = scripts 07 and 09 reuse their saved fits and redo only the post-processing |
| `MATRIKS_EXCLUDE_ITEMS` | empty | MatriKS88 item numbers (34–51) excluded from Models 1–2 (sensitivity analysis) |
| `MATRIKS_RIC_POLICY` | `R` | Error type of `r.ic` responses: `R` (paper), `IC` or `drop` (sensitivity analysis; scripts 07 and 09 require `R`) |

**Paths with non-ASCII characters (e.g. accented letters).** CmdStan's output files are written to
a temporary folder by default. If R reports that the folder contains non-ASCII characters, set an
ASCII-only folder before running, e.g. `Sys.setenv(MATRIKS_STAN_CSV_DIR = "C:/stan_csv")`.

### Running on simulated data

`tools/simulate_example_data.R` writes two files with the columns, types and response labels of the
study files, filled with random values, to `data_example/`. They can be used to check the
installation and the code without the study data; the results say nothing about the study.

```r
source("tools/simulate_example_data.R")
Sys.setenv(MATRIKS_DATA_DIR = "data_example", MATRIKS_OUTPUT_DIR = "outputs_example",
           MATRIKS_QUICK_TEST = "1")
source("run_all.R")
```

## Outputs and correspondence with the paper

| File in `outputs/` | Paper |
|---|---|
| `Table1_clinical_groups.xlsx`, `Table2_patients_vs_controls.xlsx` | Tables 1–2 |
| `accuracy_errors_anova.xlsx`, `Figure2_accuracy_EMM.jpg` | ANOVAs, EMMs and contrasts; Figure 2 |
| `error_type_proportions.xlsx`, `Figure3_error_type_proportions.jpg` | Error-type percentages; Figure 3 |
| `matched_sample_balance.xlsx` | Table 4 (matched sample) |
| `matched_participants.csv`, `data_and_matching_outputs.xlsx` | Matched triplets (read by scripts 07 and 09); matching summaries; Model 2 VIFs |
| `model2_missing_covariates.xlsx` | Table E (Model 2: patients included vs excluded for missing cognitive scores) |
| `multinomial_model_results.xlsx`, `diagnostics/` | Model 1–2 coefficients, random-effect SDs and sampling diagnostics |
| `model_contrasts.xlsx` | Posterior contrasts (Appendix A, Tables A1–A2) |
| `Figure4_model1_forest.png`, `Figure5_model2_forest.png` | Figures 4–5 |
| `symptom_severity_regression.xlsx`, `Figure6_severity_regression.jpg` | Symptom-severity regression; Figure 6 |
| `unsystematic_responding_descriptives.xlsx`, `FigureC_response_times.png` | Response times (Appendix C, Table C1 and Figure C) |
| `unsystematic_responding_models.xlsx` | Model 1 with response-time covariates (Table C2) |
| `random_choice_mixture_models.xlsx` | Random-choice mixture models and PSIS-LOO comparison (Appendix D, Tables D1–D2) |
| `prior_sensitivity.xlsx` | Models 1–2 under alternative pooling priors (Appendix F, Table F) |
| `distractor_distance_descriptives.xlsx` | Distance of the 126 distractors by error category and type; distance of the chosen distractors by group |
| `distractor_distance_models.xlsx` | Model 1 with distractor distance: group contrasts, distance slopes and PSIS-LOO (Appendix H, Tables H1–H2) |
| `sensitivity_tables.xlsx` | Main analysis vs sensitivity analyses (Appendix B, Table B; Appendix I, Table I), written by `run_sensitivity_analyses.R` |
| `model1_fit.rds`, `model2_fit.rds`, `model1rt_*_fit.rds`, `model1mix_ref_fit.rds`, `mixture_*_fit.rds`, `model1dist_*_fit.rds` | Complete posterior draws (cmdstanr objects) |
| `session_info.txt` | R, package and CmdStan versions |

The coding rules of the distractors (Appendix G, Table G) are documented in `data/README.md`
together with `data/distractor_distance.csv`.

## Implementation notes

- **Sample.** Controls: participants aged ≥ 18 with no recorded diagnosis (rows with gender `B` or
  missing gender are excluded), n = 371. Patients: n = 59 after excluding three participants
  (codes listed in `R/00_setup.R`) [reason to be documented by the authors].
- **Items.** All analyses use the 18 items shared by the two MatriKS versions; the correspondence
  between item numbers is defined in `R/00_setup.R`. Both versions were administered as fixed forms.
  For two patients the administration was incomplete (16 and 7 of the 18 items).
- **Accuracy.** In the ANOVAs, accuracy is the number of correct responses on the 18 items (skipped
  items count as not correct). For matching, accuracy is the proportion of correct responses among
  answered (non-skipped) items.
- **Matching.** Optimal matching in triplets (one control, one bipolar and one schizophrenia
  participant), anchored on the smallest group (schizophrenia, n = 28), function
  `match_three_groups()` in `R/00_setup.R`. For controls and for bipolar patients separately, the
  anchors are paired with distinct participants by solving the linear assignment problem
  (`clue::solve_LSAP`) on the squared difference in accuracy; pairs more than 0.10 apart are not
  allowed (caliper, `MATCH_CALIPER`) and ties in accuracy are broken by age (`MATCH_AGE_WEIGHT`).
  Triplets whose accuracy range exceeds the caliper are dropped. All 28 triplets meet the caliper,
  so Model 1 includes 28 participants per group (3 bipolar patients unmatched). The procedure is
  deterministic.
- **`r.ic` responses.** In three items (MatriKS88 items 37, 42 and 44; MatriKS 12+ items 15, 38 and
  40), MatriKS records the choice of either repetition distractor with a single label, `r.ic`. In
  the main analysis these responses are repetition (R) errors (`MATRIKS_RIC_POLICY = "R"`, the
  default). Two sensitivity analyses check this choice: Model 1 with the three items excluded and
  with distractor distance (script 09, Appendix H), and Models 1–2 with the `r.ic` responses counted
  as incomplete correlate (IC) errors (`run_sensitivity_analyses.R`, Appendix I). In script 09 the
  likelihood of an `r.ic` response is the probability of choosing either repetition distractor.
- **Error-type percentages (Figure 3).** Descriptive, full (unmatched) sample: errors pooled within
  each group; 95% CIs by participant-level bootstrap (5,000 resamples).
- **Model 2 covariates.** Attention and Digit Span scores are z-standardised across all patients
  with a valid score on that measure, before excluding patients with missing values on any
  covariate.
- **Response times (script 07).** Item response times (seconds) are centred on the item mean of
  the matched controls, scaled by their SD, and split into a between-participant part (mean speed)
  and a within-participant part; a sensitivity analysis removes responses above the 99th
  percentile of the matched sample. Responses without a recorded time are excluded from the
  response-time models, which use the same Stan program as Models 1–2.
- **Random-choice mixture models (script 07).** Each response is a mixture of the systematic
  process of Model 1 (probability 1 − λ) and of a random choice among the response options of the
  item (probability λ), under which each error type is chosen in proportion to its number of
  distractors; λ differs by group and between participants (`stan/mixture_random_choice.stan`).
  Variants with errors only or with all answered items, with or without diagnosis in the
  systematic part, and with random choice allowed to favour the D distractor are compared with
  PSIS-LOO (`loo`).
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
  and outputs were checked by the authors (see the declaration in the paper).

## Citation

If you use this code, please cite the paper above and this repository ([Zenodo DOI]).

## Licence

[To be chosen by the authors, e.g. MIT for the code. Add the licence text in a `LICENSE` file.]

## Contact

[Corresponding author, e-mail]
