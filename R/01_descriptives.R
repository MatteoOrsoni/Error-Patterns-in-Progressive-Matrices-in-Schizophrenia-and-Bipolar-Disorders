###############################################################################
# 01_descriptives.R
# Demographic and clinical characteristics (paper: Tables 1, 2 and K1).
#   Table 1: Bipolar vs Schizophrenia (Welch t-tests, chi-square test)
#   Table 2: all psychiatric patients vs healthy controls
#   Table K1 (Appendix K): clinical characteristics of all 62 recruited patients
#            (hospitalisations, days from admission to testing, diagnosis in
#            the schizophrenia group, episode at admission in the BD group)
# Output: outputs/Table1_clinical_groups.xlsx, outputs/Table2_patients_vs_controls.xlsx,
#         outputs/TableK1_clinical_characteristics.xlsx,
#         outputs/TableK2_psychotropic_medication.xlsx
#   Table K2 (Appendix K): psychotropic medication during the index admission,
#            from data/pharmacotherapy.xlsx (one row per patient, N = 62)
###############################################################################

if (!exists("MATRIKS_SETUP_DONE")) source(file.path("R", "00_setup.R"))

controls <- load_controls()
patients <- load_patients() %>%
  mutate(BPRS_diff = BPRS_PRE - BPRS_POST)   # treatment effect (admission - discharge)

cat("Patients analysed:", nrow(patients), "| Controls analysed:", nrow(controls), "\n")
print(table(patients$Diagnosi))

# ---- Formatting helpers -----------------------------------------------------
mean_sd <- function(x) sprintf("%.2f \u00b1 %.2f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))

female_n_pct <- function(g) {
  n_f <- sum(g == "F", na.rm = TRUE)
  sprintf("%d (%d%%)", n_f, round(100 * n_f / sum(!is.na(g))))
}

fmt_t <- function(tt) {
  sprintf("t(%.2f) = %.2f, p %s", tt$parameter, tt$statistic,
          ifelse(tt$p.value < .001, "< .001", paste("=", fmt_p(tt$p.value))))
}

fmt_chi <- function(ct) {
  chi <- unname(ct$statistic)
  chi_txt <- if (chi < .01) sub("^0", "", sprintf("%.3f", chi)) else sprintf("%.2f", chi)
  sprintf("\u03c7\u00b2 = %s, p %s", chi_txt,
          ifelse(ct$p.value < .001, "< .001", paste("=", fmt_p(ct$p.value))))
}

# ---- Table 1: Bipolar vs Schizophrenia (Welch t-tests; chi-square with Yates) ----
bd <- filter(patients, Diagnosi == "Bipolar")
sz <- filter(patients, Diagnosi == "Schizophrenic")

tests_t1 <- list(
  age        = t.test(age            ~ Diagnosi, data = patients),
  education  = t.test(anni_scolarita ~ Diagnosi, data = patients),
  bprs_pre   = t.test(BPRS_PRE       ~ Diagnosi, data = patients),
  bprs_post  = t.test(BPRS_POST      ~ Diagnosi, data = patients),
  bprs_diff  = t.test(BPRS_diff      ~ Diagnosi, data = patients),
  honos_pre  = t.test(HONOS_PRE      ~ Diagnosi, data = patients),
  honos_post = t.test(HONOS_POST     ~ Diagnosi, data = patients)
)
chi_t1 <- chisq.test(table(patients$Diagnosi, patients$gender))
invisible(lapply(c(tests_t1, list(sex = chi_t1)), print))

table1 <- tibble(
  Variable = c("Age (Years)", "Sex (Female, n%)", "Education (Years)", "BPRS Admission",
               "BPRS Discharge", "Treatment Effect", "HoNOS Admission", "HoNOS Discharge"),
  !!sprintf("Bipolar (n = %d)", nrow(bd)) := c(
    mean_sd(bd$age), female_n_pct(bd$gender), mean_sd(bd$anni_scolarita),
    mean_sd(bd$BPRS_PRE), mean_sd(bd$BPRS_POST), mean_sd(bd$BPRS_diff),
    mean_sd(bd$HONOS_PRE), mean_sd(bd$HONOS_POST)),
  !!sprintf("Schizophrenia (n = %d)", nrow(sz)) := c(
    mean_sd(sz$age), female_n_pct(sz$gender), mean_sd(sz$anni_scolarita),
    mean_sd(sz$BPRS_PRE), mean_sd(sz$BPRS_POST), mean_sd(sz$BPRS_diff),
    mean_sd(sz$HONOS_PRE), mean_sd(sz$HONOS_POST)),
  `Statistical Analysis` = c(
    fmt_t(tests_t1$age), fmt_chi(chi_t1), fmt_t(tests_t1$education),
    fmt_t(tests_t1$bprs_pre), fmt_t(tests_t1$bprs_post), fmt_t(tests_t1$bprs_diff),
    fmt_t(tests_t1$honos_pre), fmt_t(tests_t1$honos_post))
)

# ---- Table 2: all patients vs healthy controls ----------------------------------
tests_t2 <- list(
  age       = t.test(patients$age, controls$age),
  education = t.test(patients$anni_scolarita, controls$anni_scolarita)
)
chi_t2 <- chisq.test(rbind(table(patients$gender), table(controls$gender)))
invisible(lapply(c(tests_t2, list(sex = chi_t2)), print))

table2 <- tibble(
  Variable = c("Age (Years)", "Sex (Female, n%)", "Education (Years)"),
  !!sprintf("Psychiatric (n = %d)", nrow(patients)) := c(
    mean_sd(patients$age), female_n_pct(patients$gender), mean_sd(patients$anni_scolarita)),
  !!sprintf("Controls (n = %d)", nrow(controls)) := c(
    mean_sd(controls$age), female_n_pct(controls$gender), mean_sd(controls$anni_scolarita)),
  `Statistical Analysis` = c(fmt_t(tests_t2$age), fmt_chi(chi_t2), fmt_t(tests_t2$education))
)

print(table1, width = Inf)
print(table2, width = Inf)
write_xlsx(table1, file.path(OUTPUT_DIR, "Table1_clinical_groups.xlsx"))
write_xlsx(table2, file.path(OUTPUT_DIR, "Table2_patients_vs_controls.xlsx"))

# ---- Helpers for Tables K1 and K2 ---------------------------------------------------
# Formatting helpers (one decimal; integers shown without decimals in IQR/range).
fmt_num  <- function(v) ifelse(abs(v - round(v)) < 1e-9, sprintf("%d", as.integer(round(v))),
                               sub("0+$", "", sprintf("%.2f", v)))
k1_mean  <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
k1_med   <- function(x) { q <- quantile(x, c(.25, .5, .75), na.rm = TRUE)
                          sprintf("%s (%s\u2013%s)", fmt_num(q[2]), fmt_num(q[1]), fmt_num(q[3])) }
k1_range <- function(x) sprintf("%s\u2013%s", fmt_num(min(x, na.rm = TRUE)), fmt_num(max(x, na.rm = TRUE)))
k1_npct  <- function(k, n) sprintf("%d (%.1f)", k, 100 * k / n)
DASH <- "\u2014"

# ---- Table K1 (Appendix K): clinical characteristics ------------------------------
# Descriptive only, on ALL recruited patients (N = 62), i.e. including the three
# bipolar patients excluded from the analyses (EXCLUDED_PATIENT_CODES), as stated
# in the note to Table K1. The raw file is therefore read again without exclusions.
# Quantiles use R's default (type 7), as in the published table.
K1_COLUMNS <- c("external_code", "Diagnosi", "num_hospitalizations", "admission_to_test_days",
                "schizophrenia_status", "bipolar_status", "bipolar_psychotic_symptoms")

psy_all <- read_excel(PSY_FILE)
names(psy_all) <- trimws(names(psy_all))   # e.g. "bipolar_psychotic_symptoms " (trailing space)

if (!all(K1_COLUMNS %in% names(psy_all))) {
  message("Table K1 skipped: data_psy.xlsx lacks the clinical columns ",
          paste(setdiff(K1_COLUMNS, names(psy_all)), collapse = ", "), ".")
} else {
  # Harmonise the free-text clinical codes of the records (Italian/English
  # spellings, typos, capitalisation) into the categories reported in the paper.
  recode_sz_status <- function(x) {
    x <- str_to_lower(str_squish(x))
    case_when(
      is.na(x)                                    ~ NA_character_,
      str_detect(x, "^schizoaff")                 ~ "Schizoaffective disorder",
      str_detect(x, "paranoid")                   ~ "Paranoid schizophrenia",
      x %in% c("schizophrenic", "schizofrenia")   ~ "Schizophrenia, subtype not specified",
      str_detect(x, "psicosi nas|psychosis nos")  ~ "Psychosis NOS",
      str_detect(x, "schizotip|schizotyp")        ~ "Schizotypal disorder",
      TRUE                                        ~ paste0("UNRECOGNISED: ", x)
    )
  }
  recode_bd_status <- function(x) {
    x <- str_to_lower(str_squish(x))
    case_when(
      is.na(x)                               ~ NA_character_,
      str_detect(x, "^(ipo|hypo)man")        ~ "Hypomanic",
      str_detect(x, "^man")                  ~ "Manic",       # "Maniacal", "Manical"
      str_detect(x, "^mix|^mist")            ~ "Mixed",
      str_detect(x, "^dep")                  ~ "Depressive",
      TRUE                                   ~ paste0("UNRECOGNISED: ", x)
    )
  }
  SZ_LEVELS <- c("Schizoaffective disorder", "Paranoid schizophrenia",
                 "Schizophrenia, subtype not specified", "Psychosis NOS", "Schizotypal disorder")
  BD_LEVELS <- c("Manic", "Hypomanic", "Mixed", "Depressive")

  psy_all <- psy_all %>%
    mutate(
      group     = factor(Diagnosi, levels = c("Schizophrenic", "Bipolar")),
      n_hosp    = as.numeric(num_hospitalizations),
      days_test = as.numeric(admission_to_test_days),
      sz_status = recode_sz_status(schizophrenia_status),
      bd_status = recode_bd_status(bipolar_status),
      bd_psychotic = str_to_lower(str_squish(bipolar_psychotic_symptoms)) %in% c("yes", "si", "s\u00ec")
    )

  # Checks: unknown codes, codes in the wrong group, consistency of the dates.
  unknown <- unique(c(grep("^UNRECOGNISED", psy_all$sz_status, value = TRUE),
                      grep("^UNRECOGNISED", psy_all$bd_status, value = TRUE)))
  if (length(unknown) > 0) stop("Table K1: unrecognised clinical codes: ", paste(unknown, collapse = "; "), call. = FALSE)
  if (any(is.na(psy_all$group))) stop("Table K1: Diagnosi must be Bipolar or Schizophrenic.", call. = FALSE)
  if (any(is.na(psy_all$sz_status[psy_all$group == "Schizophrenic"])) ||
      any(!is.na(psy_all$sz_status[psy_all$group == "Bipolar"])) ||
      any(is.na(psy_all$bd_status[psy_all$group == "Bipolar"])) ||
      any(!is.na(psy_all$bd_status[psy_all$group == "Schizophrenic"]))) {
    warning("Table K1: schizophrenia_status / bipolar_status not consistent with Diagnosi.", call. = FALSE)
  }
  if (all(c("admission_date", "test_administration_date") %in% names(psy_all))) {
    adm  <- as.Date(psy_all$admission_date)
    test <- as.Date(psy_all$test_administration_date)
    days_from_dates <- as.numeric(test - adm)
    if (any(days_from_dates != psy_all$days_test, na.rm = TRUE)) {
      warning("Table K1: admission_to_test_days differs from test_administration_date - admission_date for ",
              paste(psy_all$external_code[which(days_from_dates != psy_all$days_test)], collapse = ", "), call. = FALSE)
    }
    fmt_date <- function(x) sub("^0", "", format(x, "%d %B %Y"))
    cat(sprintf("\nTable K1 note: patients admitted between %s and %s, tested between %s and %s.\n",
                fmt_date(min(adm, na.rm = TRUE)), fmt_date(max(adm, na.rm = TRUE)),
                fmt_date(min(test, na.rm = TRUE)), fmt_date(max(test, na.rm = TRUE))))
  }


  sz_k1 <- filter(psy_all, group == "Schizophrenic")
  bd_k1 <- filter(psy_all, group == "Bipolar")
  groups_k1 <- list(sz_k1, bd_k1, psy_all)

  continuous_block <- function(var, label) {
    tibble(Characteristic = c(label, "Mean (SD)", "Median (IQR)", "Range"),
           SZ  = c("", k1_mean(sz_k1[[var]]),  k1_med(sz_k1[[var]]),  k1_range(sz_k1[[var]])),
           BD  = c("", k1_mean(bd_k1[[var]]),  k1_med(bd_k1[[var]]),  k1_range(bd_k1[[var]])),
           Tot = c("", k1_mean(psy_all[[var]]), k1_med(psy_all[[var]]), k1_range(psy_all[[var]])))
  }
  heading_row <- function(label) tibble(Characteristic = label, SZ = "", BD = "", Tot = "")

  sz_rows <- tibble(Characteristic = SZ_LEVELS,
                    SZ  = map_chr(SZ_LEVELS, ~ k1_npct(sum(sz_k1$sz_status == .x), nrow(sz_k1))),
                    BD  = DASH, Tot = DASH)
  bd_rows <- tibble(Characteristic = BD_LEVELS,
                    SZ  = DASH,
                    BD  = map_chr(BD_LEVELS, ~ k1_npct(sum(bd_k1$bd_status == .x), nrow(bd_k1))),
                    Tot = DASH)
  episode_label <- c(Manic = "manic", Hypomanic = "hypomanic", Mixed = "mixed", Depressive = "depressive")
  psych_rows <- bind_rows(
    tibble(Characteristic = "Psychotic symptoms", SZ = DASH,
           BD = k1_npct(sum(bd_k1$bd_psychotic), nrow(bd_k1)), Tot = DASH),
    map_dfr(BD_LEVELS, function(lv) {
      d <- filter(bd_k1, bd_status == lv)
      tibble(Characteristic = sprintf("in %s episodes, n/N (%%)", episode_label[[lv]]), SZ = DASH,
             BD = sprintf("%d/%d (%.1f)", sum(d$bd_psychotic), nrow(d),
                          if (nrow(d) > 0) 100 * mean(d$bd_psychotic) else NA_real_),
             Tot = DASH)
    })
  )

  tableK1 <- bind_rows(
    heading_row("Hospitalisation"),
    continuous_block("n_hosp",    "Number of psychiatric admissions (a)"),
    continuous_block("days_test", "Days from admission to testing (b)"),
    heading_row("Diagnosis, schizophrenia group, n (%)"),
    sz_rows,
    heading_row("Episode at admission, bipolar disorder group, n (%)"),
    bd_rows,
    psych_rows
  )
  names(tableK1) <- c("Characteristic",
                      sprintf("Schizophrenia (n = %d)", nrow(sz_k1)),
                      sprintf("Bipolar disorder (n = %d)", nrow(bd_k1)),
                      sprintf("Total (N = %d)", nrow(psy_all)))

  print(tableK1, n = Inf, width = Inf)
  write_xlsx(tableK1, file.path(OUTPUT_DIR, "TableK1_clinical_characteristics.xlsx"))
}

# ---- Table K2 (Appendix K): psychotropic medication ---------------------------------
# Descriptive only, on all 62 recruited patients, from data/pharmacotherapy.xlsx
# (one row per patient; Yes/No indicators per drug class; codebook in the file).
# Drugs administered at any time during the index admission (medication on the
# day of testing was not available). Diagnostic group taken from data_psy.xlsx.
PHARMA_FILE <- file.path(DATA_DIR, "pharmacotherapy.xlsx")
K2_ROWS <- tribble(
  ~section,                          ~column,                     ~label,
  "Antipsychotics",                  "any_antipsychotic",         "Any antipsychotic",
  "Antipsychotics",                  "two_plus_antipsychotics",   "Two or more antipsychotics during admission",
  "Antipsychotics",                  "fga",                       "First-generation antipsychotic",
  "Antipsychotics",                  "sga",                       "Second-generation antipsychotic",
  "Antipsychotics",                  "clozapine",                 "Clozapine",
  "Antipsychotics",                  "short_acting_im_iv",        "Short-acting IM/IV psychotropic",
  "Antipsychotics",                  "lai",                       "Long-acting injectable antipsychotic (LAI)",
  "Mood stabilisers",                "any_mood_stabiliser",       "Any mood stabiliser",
  "Mood stabilisers",                "lithium",                   "Lithium",
  "Mood stabilisers",                "valproate",                 "Valproate",
  "Mood stabilisers",                "lamotrigine",               "Lamotrigine",
  "Mood stabilisers",                "two_plus_mood_stabilisers", "Two or more mood stabilisers",
  "Other psychotropics",             "any_other_antiepileptic",   "Other antiepileptics (lacosamide, topiramate, primidone, cenobamate)",
  "Other psychotropics",             "gabapentinoid",             "Gabapentinoids (pregabalin, gabapentin)",
  "Other psychotropics",             "any_benzodiazepine_z",      "Benzodiazepines / Z-drugs",
  "Other psychotropics",             "any_antidepressant",        "Antidepressants",
  "Other psychotropics",             "anticholinergic",           "Anticholinergics (biperiden)",
  "Other psychotropics",             "other_cns_drug",            "Other CNS-active drugs (hydroxyzine, memantine, methadone)",
  "Medication burden and sedation",  "any_sedating_drug",         "At least one markedly sedating drug (c)",
  "Medication burden and sedation",  "bdz_only_sedative",         "Benzodiazepine / Z-drug as the only sedating drug",
  "Medication burden and sedation",  "complex_polypharmacy",      "Complex polypharmacy (\u22654 psychotropics)"
)

if (!file.exists(PHARMA_FILE)) {
  message("Table K2 skipped: ", PHARMA_FILE, " not found.")
} else {
  pharma <- read_excel(PHARMA_FILE, sheet = "pharmacotherapy")
  names(pharma) <- trimws(names(pharma))
  check_columns(pharma, c("external_code", "n_psychotropics", K2_ROWS$column), PHARMA_FILE)

  yes_no <- function(x) {
    x <- str_to_lower(str_squish(as.character(x)))
    if (!all(x %in% c("yes", "no"))) stop("Table K2: Yes/No columns contain other values.", call. = FALSE)
    x == "yes"
  }
  pharma <- pharma %>%
    mutate(external_code = as.character(external_code),
           across(all_of(K2_ROWS$column), yes_no),
           n_psychotropics = as.numeric(n_psychotropics))

  # Same patients as data_psy.xlsx, with the group from data_psy.xlsx.
  psy_groups <- read_excel(PSY_FILE) %>% transmute(external_code = as.character(external_code), Diagnosi)
  if (!setequal(pharma$external_code, psy_groups$external_code) || anyDuplicated(pharma$external_code)) {
    stop("Table K2: patient codes in pharmacotherapy.xlsx differ from data_psy.xlsx.", call. = FALSE)
  }
  pharma <- left_join(pharma, psy_groups, by = "external_code")
  if ("group" %in% names(pharma)) {
    grp_check <- c("Schizophrenia" = "Schizophrenic", "Bipolar disorder" = "Bipolar")[pharma$group]
    if (any(grp_check != pharma$Diagnosi, na.rm = TRUE)) {
      warning("Table K2: 'group' in pharmacotherapy.xlsx differs from Diagnosi in data_psy.xlsx for ",
              paste(pharma$external_code[which(grp_check != pharma$Diagnosi)], collapse = ", "), call. = FALSE)
    }
  }
  # Internal consistency of the indicators with the counts, when the counts are present.
  count_checks <- list(any_antipsychotic = "n_antipsychotics", any_other_antiepileptic = "n_other_antiepileptics",
                       any_benzodiazepine_z = "n_benzodiazepines_z", any_antidepressant = "n_antidepressants",
                       any_sedating_drug = "n_sedating_drugs")
  for (flag in names(count_checks)) {
    n_col <- count_checks[[flag]]
    if (n_col %in% names(pharma) && any(pharma[[flag]] != (as.numeric(pharma[[n_col]]) > 0))) {
      warning(sprintf("Table K2: '%s' not consistent with '%s'.", flag, n_col), call. = FALSE)
    }
  }
  if (any(pharma$complex_polypharmacy != (pharma$n_psychotropics >= 4))) {
    warning("Table K2: 'complex_polypharmacy' not consistent with n_psychotropics >= 4.", call. = FALSE)
  }

  sz_k2 <- filter(pharma, Diagnosi == "Schizophrenic")
  bd_k2 <- filter(pharma, Diagnosi == "Bipolar")
  k2_npct <- function(d, col) sprintf("%d (%.1f)", sum(d[[col]]), 100 * mean(d[[col]]))
  k2_med  <- function(x) { q <- quantile(x, c(.25, .5, .75))
                           sprintf("%s (%s\u2013%s)", fmt_num(q[2]), fmt_num(q[1]), fmt_num(q[3])) }

  tableK2 <- map_dfr(unique(K2_ROWS$section), function(sec) {
    rows <- filter(K2_ROWS, section == sec)
    bind_rows(
      tibble(Variable = sec, SZ = "", BD = "", Tot = ""),
      tibble(Variable = rows$label,
             SZ  = map_chr(rows$column, ~ k2_npct(sz_k2, .x)),
             BD  = map_chr(rows$column, ~ k2_npct(bd_k2, .x)),
             Tot = map_chr(rows$column, ~ k2_npct(pharma, .x)))
    )
  }) %>%
    bind_rows(tibble(
      Variable = c("Number of psychotropics, median (IQR)", "Number of psychotropics, range"),
      SZ  = c(k2_med(sz_k2$n_psychotropics),  k1_range(sz_k2$n_psychotropics)),
      BD  = c(k2_med(bd_k2$n_psychotropics),  k1_range(bd_k2$n_psychotropics)),
      Tot = c(k2_med(pharma$n_psychotropics), k1_range(pharma$n_psychotropics))))
  names(tableK2) <- c("Variable",
                      sprintf("Schizophrenia (n = %d)", nrow(sz_k2)),
                      sprintf("Bipolar disorder (n = %d)", nrow(bd_k2)),
                      sprintf("Total (N = %d)", nrow(pharma)))

  print(tableK2, n = Inf, width = Inf)
  write_xlsx(tableK2, file.path(OUTPUT_DIR, "TableK2_psychotropic_medication.xlsx"))
}
