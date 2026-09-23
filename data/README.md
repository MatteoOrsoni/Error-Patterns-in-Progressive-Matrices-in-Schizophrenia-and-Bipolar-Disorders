# Data

Place the two data files in this folder:

| File | Content |
|---|---|
| `data_adult.xlsx` | Healthy participants (MatriKS 12+) |
| `data_psy.xlsx` | Psychiatric inpatients (MatriKS88) |

**Data availability:** [to be completed by the authors: open access under a CC-BY 4.0 licence / restricted access for editors and reviewers, available upon reasonable request, and the reason (ethics approval / informed consent)].

The analysis code uses only the columns listed below. `tools/make_minimal_dataset.R` creates files containing only these columns from the full study files.

## `data_adult.xlsx` (one row per participant)

| Column | Type | Use in the code |
|---|---|---|
| `new_id` | text | Participant identifier |
| `age` | numeric (years) | Inclusion (`age >= 18`); Table 2 |
| `gender` | `F` / `M` / `B` | Table 2; rows with `B` (other) or missing values are excluded. |
| `diagnostic` | text / empty | Rows with any recorded diagnosis are excluded (empty = healthy control) |
| `anni_scolarita` | numeric (years of education), `NULL` = missing | Table 2 |
| `response_mat50_item.<n>` for n = 4, 13, 14, 15, 16, 27, 30, 33, 38–47 | response label (see below) | The 18 items shared with MatriKS88 |

## `data_psy.xlsx` (one row per patient)

| Column | Type | Use in the code |
|---|---|---|
| `external_code` | text (e.g. `BSPDC01`) | Patient identifier; codes BSPDC40, BSPDC41 and BSPDC45 are excluded from all analyses |
| `Diagnosi` | `Bipolar` / `Schizophrenic` | Diagnostic group |
| `age`, `gender`, `anni_scolarita` | as above | Tables 1–2 |
| `BPRS_PRE`, `BPRS_POST` | numeric | BPRS total score at admission / discharge (Table 1; regression) |
| `HONOS_PRE`, `HONOS_POST` | numeric | HoNOS total score at admission / discharge (Table 1; regression) |
| `Matrici_Attentive` | numeric (comma or dot decimal) | Attentive Matrices Test score (Model 2) |
| `SPAN_diretto`, `SPAN_inverso` | numeric | Forward / backward Digit Span (Model 2) |
| `response_mat88_item.<n>` for n = 34–51 | response label | The 18 items shared with MatriKS 12+ (same order as the adult items above) |

[Score type of the cognitive tests (raw or age/education-corrected) to be documented by the authors.]

## Response labels

| Label(s) | Meaning |
|---|---|
| `correct` | Correct response |
| `r.diag`, `r.left`, `r.top` | Repetition (R) error |
| `d.union`, `diff`, `diff1`, `diff2` | Difference (D) error |
| `ic.flip`, `ic.inc`, `ic.neg`, `ic.scale` | Incomplete correlate (IC) error |
| `wp.copy`, `wp.matrix`, `wp1` | Wrong principle (WP) error |
| `r.ic` | Distractor coded `r.ic`; classified as a repetition (R) error in all error-type analyses (see README, "Implementation notes") |
| `skip` | No response (excluded from error analyses) |
