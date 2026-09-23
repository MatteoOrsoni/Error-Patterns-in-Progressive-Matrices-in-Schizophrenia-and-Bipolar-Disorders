###############################################################################
# 06_symptom_severity.R
# Linear regression of total MatriKS errors (18 common items) on symptom
# severity at admission (BPRS, HoNOS), diagnosis and their interactions,
# in the psychiatric sample (complete cases).
#   - BPRS and HoNOS are mean-centred, so the diagnosis coefficient is the
#     Schizophrenia - Bipolar difference at average severity
#   - influential observations: Cook's distance > 4/n
#   - sensitivity analysis: model refitted without the influential cases
#   - Figure 6: total errors vs BPRS (a) and HoNOS (b), influential cases circled
# Outputs: outputs/symptom_severity_regression.xlsx,
#          outputs/Figure6_severity_regression.jpg
###############################################################################

if (!exists("MATRIKS_SETUP_DONE")) source(file.path("R", "00_setup.R"))

patients_scored <- build_scored_dataset()$patients

data_severity <- patients_scored %>%
  select(external_code, Diagnosi, BPRS_PRE, HONOS_PRE, total_errors) %>%
  mutate(BPRS_PRE  = as.numeric(as.character(BPRS_PRE)),
         HONOS_PRE = as.integer(as.character(HONOS_PRE)),
         Diagnosi  = factor(Diagnosi, levels = c("Bipolar", "Schizophrenic"))) %>%
  drop_na(BPRS_PRE, HONOS_PRE, total_errors) %>%
  mutate(BPRS_c  = BPRS_PRE  - mean(BPRS_PRE),
         HONOS_c = HONOS_PRE - mean(HONOS_PRE))

cat("Patients with complete data:", nrow(data_severity), "\n")

severity_formula <- total_errors ~ BPRS_c * Diagnosi + HONOS_c * Diagnosi

# ---- Full model ------------------------------------------------------------------
model_severity <- lm(severity_formula, data = data_severity)
print(summary(model_severity))

# ---- Influential observations (Cook's distance > 4/n) -----------------------------
cooks_d         <- cooks.distance(model_severity)
threshold       <- 4 / length(cooks_d)
influential_idx <- unname(which(cooks_d > threshold))
influential     <- data_severity[influential_idx, c("external_code", "Diagnosi", "BPRS_PRE",
                                                    "HONOS_PRE", "total_errors")] %>%
  mutate(cooks_d = cooks_d[influential_idx])
cat("\nInfluential observations (Cook's D > 4/n =", round(threshold, 4), "):\n")
print(influential)

# ---- Sensitivity analysis without influential observations --------------------------
# Centring constants are those of the full sample (n with complete data), so the
# diagnosis coefficient is evaluated at the same severity level in both models.
data_severity_clean  <- filter(data_severity, !(row_number() %in% influential_idx))
model_severity_clean <- lm(severity_formula, data = data_severity_clean)
print(summary(model_severity_clean))

write_xlsx(
  list(full_coefficients        = tidy(model_severity),
       full_fit                 = glance(model_severity),
       influential_observations = influential,
       sensitivity_coefficients = tidy(model_severity_clean),
       sensitivity_fit          = glance(model_severity_clean)),
  file.path(OUTPUT_DIR, "symptom_severity_regression.xlsx")
)

# ---- Figure 6 --------------------------------------------------------------------
# Raw (uncentred) scores on the x-axis; one least-squares line per group.
make_regression_plot <- function(data, x_var, x_label, highlight_codes) {
  ggplot(data, aes(x = .data[[x_var]], y = total_errors, color = Diagnosi)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_point(data = filter(data, external_code %in% highlight_codes),
               aes(x = .data[[x_var]], y = total_errors),
               size = 5, shape = 1, stroke = 1.5, color = "red") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1.2) +
    scale_color_manual(values = DX_COLORS) +
    labs(x = x_label, y = "Total Errors", color = "Diagnostic Group") +
    DX_THEME
}

figure6 <- (make_regression_plot(data_severity, "BPRS_PRE",  "BPRS (Score)",  influential$external_code) |
              make_regression_plot(data_severity, "HONOS_PRE", "HoNOS (Score)", influential$external_code)) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Relationship Between Clinical Scores and Total Errors",
                  tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5))) &
  theme(legend.position = "bottom", plot.tag = element_text(size = 18, face = "bold"))

ggsave(file.path(OUTPUT_DIR, "Figure6_severity_regression.jpg"), figure6,
       width = 10, height = 5, units = "in", dpi = 300, quality = 95)
