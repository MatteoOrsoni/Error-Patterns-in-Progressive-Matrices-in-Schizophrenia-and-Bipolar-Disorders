# Data

The study data files are **not** included in the repository. To run the analyses, place them in
this folder (or set `MATRIKS_DATA_DIR` to the folder that contains them):

| File | Content |
|---|---|
| `data_adult.xlsx` | Healthy participants (MatriKS 12+); not included |
| `data_psy.xlsx` | Psychiatric inpatients (MatriKS88); not included |
| `distractor_distance.csv` | Distance of each of the 126 distractors from the correct response (included; script 09) |

**Data availability:** [to be completed by the authors: open access under a CC-BY 4.0 licence / restricted access for editors and reviewers, available upon reasonable request, and the reason (ethics approval / informed consent)].

The analysis code uses only the columns listed below. `tools/make_minimal_dataset.R` creates files
containing only these columns from the full study files (in `data_minimal/`).
`tools/simulate_example_data.R` writes files with the same columns and random values (in
`data_example/`), to run the code without the study data.

## `data_adult.xlsx` (one row per participant)

| Column | Type | Use in the code |
|---|---|---|
| `new_id` | text | Participant identifier |
| `age` | numeric (years) | Inclusion (`age >= 18`); Table 2; tie-break in the matching; Table 4 |
| `gender` | `F` / `M` / `B` | Tables 2 and 4; rows with `B` (other) or missing values are excluded |
| `diagnostic` | text / empty | Rows with any recorded diagnosis are excluded (empty = healthy control) |
| `anni_scolarita` | numeric (years of education), `NULL` = missing | Tables 2 and 4 |
| `response_mat50_item.<n>` for n = 4, 13, 14, 15, 16, 27, 30, 33, 38–47 | response label (see below) | The 18 items shared with MatriKS88 |
| `time_mat50_item.<n>`, same n | numeric (seconds, integer) | Response time per item (script 07) |

## `data_psy.xlsx` (one row per patient)

| Column | Type | Use in the code |
|---|---|---|
| `external_code` | text | Patient identifier; the three codes in `EXCLUDED_PATIENT_CODES` (`R/00_setup.R`) are excluded from all analyses |
| `Diagnosi` | `Bipolar` / `Schizophrenic` | Diagnostic group |
| `age`, `gender`, `anni_scolarita` | as above | Tables 1, 2 and 4; matching tie-break (age) |
| `BPRS_PRE`, `BPRS_POST` | numeric | BPRS total score at admission / discharge (Table 1; regression) |
| `HONOS_PRE`, `HONOS_POST` | numeric | HoNOS total score at admission / discharge (Table 1; regression) |
| `Matrici_Attentive` | numeric (comma or dot decimal) | Attentive Matrices Test score (Model 2) |
| `SPAN_diretto`, `SPAN_inverso` | numeric | Forward / backward Digit Span (Model 2) |
| `response_mat88_item.<n>` for n = 34–51 | response label | The 18 items shared with MatriKS 12+ (same order as the adult items above) |
| `time_mat88_item.<n>`, same n | numeric (seconds, integer) | Response time per item (script 07) |

Items not administered are empty cells (or `NULL`/`NA`); this occurs for two patients, whose MatriKS88 administration was incomplete.

[Score type of the cognitive tests (raw or age/education-corrected) to be documented by the authors.]

## Response labels

| Label(s) | Meaning |
|---|---|
| `correct` | Correct response |
| `r.diag`, `r.left`, `r.top` | Repetition (R) error |
| `d.union`, `diff`, `diff1`, `diff2` | Difference (D) error |
| `ic.flip`, `ic.inc`, `ic.neg`, `ic.scale` | Incomplete correlate (IC) error |
| `wp.copy`, `wp.matrix`, `wp1` | Wrong principle (WP) error |
| `r.ic` | Either repetition distractor, in three items (MatriKS 12+ items 15, 38 and 40; MatriKS88 items 37, 42 and 44). Repetition (R) error in the main analysis; incomplete correlate (IC) error in the sensitivity analysis of Appendix I (`MATRIKS_RIC_POLICY`; see the main README, "Implementation notes") |
| `skip` | Item presented but not answered (excluded from error analyses and from the matching accuracy) |

## `distractor_distance.csv` (one row per distractor, 18 items x 7)

Consensus coding of two raters (rules in Appendix G, Table G of the paper). Columns: `matriks12_item`, `matriks88_item`
(the same stimulus in the two versions), `distractor` (the response label recorded by MatriKS for that option),
`category` (R, D, IC, WP), the number of departures from the correct response on each feature (`missing`,
`added`, `shape`, `size`, `fill`, `line`, `orientation`), `distance` (their sum), `features_changed` (number of
features with at least one departure) and `matrix_type` (`slot` or `logical`). In the three items where repetition
errors are recorded with the label `r.ic` (MatriKS 12+ items 15, 38 and 40), both repetition distractors are listed
with their own label (`r.left`, `r.top`, `r.diag`); script 09 treats an `r.ic` response as the choice of one of them.
The independent codes of the two raters and the rationale for each consensus value are in
`modello_codifica_distrattori_finale_v2.xlsx` [to be added to the repository by the authors].
