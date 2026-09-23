###############################################################################
# 02_accuracy_errors.R
# Total accuracy and total errors on the 18 common items (full sample:
# 31 Bipolar, 28 Schizophrenia, 371 controls).
#   - one-way ANOVAs, estimated marginal means (EMMs) and pairwise contrasts
#   - Figure 2: accuracy density (a) and EMMs (b)
#   - Figure 3: descriptive proportions of error types by group
# Outputs: outputs/accuracy_errors_anova.xlsx, outputs/error_type_proportions.xlsx,
#          outputs/Figure2_accuracy_EMM.jpg, outputs/Figure3_error_type_proportions.jpg
###############################################################################

if (!exists("MATRIKS_SETUP_DONE")) source(file.path("R", "00_setup.R"))

data_general <- build_scored_dataset()$general
print(table(data_general$Diagnosi))

# ---- ANOVA, EMMs and pairwise contrasts --------------------------------------
# Factor levels: Bipolar, Control, Schizophrenic. Each contrast is labelled with
# the difference it computes (first group minus second group).
contrast_list <- list(
  "Bipolar - Control"       = c( 1, -1,  0),
  "Schizophrenic - Control" = c( 0, -1,  1),
  "Bipolar - Schizophrenic" = c( 1,  0, -1)
)
stopifnot(identical(levels(data_general$Diagnosi), c("Bipolar", "Control", "Schizophrenic")))

run_anova_emmeans <- function(outcome, data, contrasts) {
  fit <- aov(as.formula(paste(outcome, "~ Diagnosi")), data = data)
  emm <- emmeans(fit, ~ Diagnosi)
  con <- summary(contrast(emm, method = contrasts), infer = c(TRUE, TRUE))
  cat("\n== ANOVA:", outcome, "==\n");     print(summary(fit))
  cat("\n== EMMs:", outcome, "==\n");      print(summary(emm))
  cat("\n== Contrasts:", outcome, "==\n"); print(con)
  list(anova = fit, emmeans = emm, emm_tidy = tidy(emm),
       anova_tidy = tidy(fit), contrasts = as_tibble(as.data.frame(con)))
}

res_accuracy <- run_anova_emmeans("total_matriks", data_general, contrast_list)
res_errors   <- run_anova_emmeans("total_errors",  data_general, contrast_list)

write_xlsx(
  list(
    accuracy_anova     = res_accuracy$anova_tidy,
    accuracy_emmeans   = res_accuracy$emm_tidy,
    accuracy_contrasts = res_accuracy$contrasts,
    errors_anova       = res_errors$anova_tidy,
    errors_emmeans     = res_errors$emm_tidy,
    errors_contrasts   = res_errors$contrasts
  ),
  file.path(OUTPUT_DIR, "accuracy_errors_anova.xlsx")
)

# ---- Figure 2: accuracy density (a) and EMMs +/- 1 SE (b) ---------------------
make_density_plot <- function(data, x_var, x_label) {
  ggplot(data, aes(x = .data[[x_var]], color = Diagnosi, fill = Diagnosi)) +
    geom_density(alpha = 0.5) +
    scale_color_manual(values = DX_COLORS) +
    scale_fill_manual(values  = DX_COLORS) +
    labs(x = x_label, y = "Density", color = "Conditions", fill = "Conditions") +
    DX_THEME
}

make_means_plot <- function(df_means, y_label) {
  ggplot(df_means, aes(x = Diagnosi, y = estimate, color = Diagnosi, fill = Diagnosi)) +
    geom_point(size = 4, shape = 21, stroke = 1.5) +
    geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.2) +
    scale_color_manual(values = DX_COLORS) +
    scale_fill_manual(values  = DX_COLORS) +
    labs(x = NULL, y = y_label) +
    guides(color = "none", fill = "none") +      # a single legend (from panel a) for the figure
    DX_THEME +
    coord_flip()
}

figure2 <- (make_density_plot(data_general, "total_matriks", "Total Accuracy") |
              make_means_plot(res_accuracy$emm_tidy, "Total Accuracy")) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title      = "Distribution of Accuracy and Estimated Marginal Means by Condition",
    tag_levels = "a",
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  ) &
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 16),
    plot.tag        = element_text(size = 16, face = "bold")
  )

ggsave(file.path(OUTPUT_DIR, "Figure2_accuracy_EMM.jpg"), figure2,
       width = 13, height = 8, units = "in", dpi = 300, quality = 95)

# ---- Figure 3: descriptive proportions of error types by group ------------------
# Descriptive only (full, unmatched sample): percentage of each error type among
# all errors pooled within each group. "r.ic" responses follow RIC_POLICY, as in
# the Bayesian models (default "R": counted as repetition errors), so the four
# percentages sum to 100%. 95% CIs: percentile bootstrap resampling participants
# within each group, because errors of the same participant are not independent.
error_types  <- c("incomplete_correlate", "repetition", "difference", "wrong_principle")
error_labels <- c(incomplete_correlate = "Incomplete Correlate", repetition = "Repetition",
                  difference = "Difference", wrong_principle = "Wrong Principle")
N_BOOT <- 5000

error_counts <- data_general %>%
  mutate(repetition           = repetition + if (RIC_POLICY == "R") r_ic else 0L,
         incomplete_correlate = incomplete_correlate + if (RIC_POLICY == "IC") r_ic else 0L,
         total_errors         = total_errors - if (RIC_POLICY == "drop") r_ic else 0L)
stopifnot(all(rowSums(error_counts[, error_types]) == error_counts$total_errors))

pooled_pct <- function(counts, totals) 100 * colSums(counts) / sum(totals)

set.seed(SEED)
error_summary <- error_counts %>%
  group_by(Diagnosi) %>%
  group_modify(function(g, key) {
    counts <- as.matrix(g[, error_types]); totals <- g$total_errors
    boot <- replicate(N_BOOT, {
      i <- sample.int(nrow(g), replace = TRUE)
      pooled_pct(counts[i, , drop = FALSE], totals[i])
    })
    tibble(error_type = error_types,
           n_errors   = unname(colSums(counts)),
           pct        = unname(pooled_pct(counts, totals)),
           ci_lower   = unname(apply(boot, 1, quantile, 0.025)),
           ci_upper   = unname(apply(boot, 1, quantile, 0.975)))
  }) %>%
  ungroup() %>%
  mutate(across(c(pct, ci_lower, ci_upper), ~ round(.x, 2)))

write_xlsx(error_summary, file.path(OUTPUT_DIR, "error_type_proportions.xlsx"))
print(error_summary, n = Inf)

error_long <- error_summary %>%
  mutate(error_type = recode(error_type, !!!error_labels)) %>%
  rename(proportion = pct)

figure3 <- ggplot(error_long, aes(x = error_type, y = proportion, fill = Diagnosi)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = 0.8), width = 0.2, color = "black") +
  scale_fill_manual(values = DX_COLORS) +
  labs(title = "Error Proportion by Condition", x = "Error Type",
       y = "Proportion (%)", fill = "Diagnostic Group") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title      = element_text(size = 18, face = "bold"),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y     = element_text(size = 16),
    legend.position = "bottom"
  )

ggsave(file.path(OUTPUT_DIR, "Figure3_error_type_proportions.jpg"), figure3,
       width = 10, height = 5, units = "in", dpi = 300, quality = 95)
