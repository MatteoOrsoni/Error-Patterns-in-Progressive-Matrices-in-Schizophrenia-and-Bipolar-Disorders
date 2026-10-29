###############################################################################
# 01_descriptives.R
# Demographic and clinical characteristics (paper: Tables 1 and 2).
#   Table 1: Bipolar vs Schizophrenia (Welch t-tests, chi-square test)
#   Table 2: all psychiatric patients vs healthy controls
# Output: outputs/Table1_clinical_groups.xlsx, outputs/Table2_patients_vs_controls.xlsx
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
