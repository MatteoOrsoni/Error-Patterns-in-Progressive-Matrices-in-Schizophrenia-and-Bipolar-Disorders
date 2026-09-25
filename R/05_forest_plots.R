###############################################################################
# 05_forest_plots.R
# Forest plots of the posterior contrasts (odds ratios vs repetition errors).
#   Figure 4: Model 1 - BD vs controls, SZ vs controls, SZ vs BD
#   Figure 5: Model 2 - diagnosis effect and attention slopes (the forward and
#             backward span contrasts are reported in Table B only)
# Points = posterior medians; bars = 95% credible intervals; dashed line = OR 1.
# Input : outputs/model_contrasts.xlsx (from script 04)
# Output: outputs/Figure4_model1_forest.png, outputs/Figure5_model2_forest.png
###############################################################################

if (!exists("MATRIKS_SETUP_DONE")) source(file.path("R", "00_setup.R"))

CONTRASTS_FILE <- file.path(OUTPUT_DIR, "model_contrasts.xlsx")
ERROR_TYPES    <- c("D vs R", "IC vs R", "WP vs R")
DIFF_COLOR     <- "#4D4D4D"      # grey for differences between groups

fmt_prob <- function(p) sub("^0", "", sprintf("%.2f", p))   # .93 instead of 0.93

# blocks: named vector (contrast name in the xlsx -> label shown in the figure)
# colours: named vector (contrast name -> colour)
# geom: horizontal layout on the OR (log) scale
make_forest_plot <- function(sheet, blocks, colours, geom, title, header_size) {
  df <- read_excel(CONTRASTS_FILE, sheet = sheet) %>%
    filter(contrast %in% names(blocks), error_type %in% ERROR_TYPES) %>%
    mutate(g_idx  = match(contrast, names(blocks)),
           t_idx  = match(error_type, ERROR_TYPES),
           lab_or = sprintf("%.2f [%.2f, %.2f]", OR_median, OR_q2.5, OR_q97.5),
           lab_p  = fmt_prob(P_gt_0))

  block_h <- length(ERROR_TYPES) + 1 + 1.2           
  df  <- mutate(df, y = -((g_idx - 1) * block_h) - t_idx)
  hdr <- tibble(y = -((seq_along(blocks) - 1) * block_h), label = unname(blocks))
  ref <- tibble(y_top = hdr$y - 0.5, y_bot = hdr$y - length(ERROR_TYPES) - 0.5)
  y_top <- 1.2; y_bot <- min(df$y) - 0.7
  h_center <- (mean(log10(geom$x_range)) - log10(geom$x_min)) / (log10(geom$x_max) - log10(geom$x_min))

  p <- ggplot(df) +
    geom_segment(data = ref, aes(x = 1, xend = 1, y = y_top, yend = y_bot),
                 linetype = "dashed", colour = "grey45", linewidth = 0.6) +
    geom_segment(aes(x = OR_q2.5, xend = OR_q97.5, y = y, yend = y, colour = contrast),
                 linewidth = 1.6, lineend = "round") +
    geom_point(aes(x = OR_median, y = y, colour = contrast), shape = 16, size = 4.6) +
    geom_text(data = hdr, aes(x = geom$lab_x, y = y, label = label),
              hjust = 0, size = header_size, fontface = "bold") +
    geom_text(aes(x = geom$lab_x * geom$row_offset, y = y, label = error_type),
              hjust = 0, size = 5.6, colour = "grey25") +
    geom_text(aes(x = geom$x_or, y = y, label = lab_or), hjust = 0, size = 5.6, colour = "grey15") +
    geom_text(aes(x = geom$x_p,  y = y, label = lab_p),  hjust = 1, size = 5.6, colour = "grey15") +
    annotate("text", x = geom$x_or, y = y_top, label = "OR [95% CrI]", hjust = 0, size = 5.4, colour = "grey35") +
    annotate("text", x = geom$x_p, y = y_top, label = 'paste("P(", beta, " > 0)")', parse = TRUE,
             hjust = 1, size = 5.4, colour = "grey35") +
    annotate("segment", x = 0.93, xend = 0.80, y = y_top, yend = y_top, colour = "grey35",
             linewidth = 0.5, arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
    annotate("text", x = 0.78, y = y_top, label = "relatively more R errors",
             hjust = 1, size = 5.0, colour = "grey35") +
    annotate("segment", x = 1.08, xend = 1.25, y = y_top, yend = y_top, colour = "grey35",
             linewidth = 0.5, arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
    annotate("text", x = 1.28, y = y_top, label = "relatively more D, IC, WP",
             hjust = 0, size = 5.0, colour = "grey35") +
    scale_colour_manual(values = colours, guide = "none") +
    scale_x_log10(breaks = geom$breaks, labels = as.character(geom$breaks),
                  limits = c(geom$x_min, geom$x_max), expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_bot, y_top + 0.6), expand = c(0, 0)) +
    labs(title = title, x = "Odds ratio vs repetition errors (log scale)", y = NULL) +
    DX_THEME +
    theme(plot.title         = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.title.x       = element_text(size = 16, face = "bold", hjust = h_center,
                                            margin = margin(t = 10)),
          axis.text.x        = element_text(size = 16),
          axis.text.y        = element_blank(),
          axis.ticks.x       = element_line(colour = "grey50"),
          axis.line.x        = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.4),
          plot.margin        = margin(10, 15, 10, 10))
  list(plot = p, height = 7.2 * (y_top - y_bot) / 15.4) 
}

# ---- Figure 4: Model 1 ---------------------------------------------------------
fig4 <- make_forest_plot(
  sheet  = "model1_group_contrasts",
  blocks = c("BD vs Control" = "Bipolar vs controls",
             "SZ vs Control" = "Schizophrenia vs controls"
            ),
  colours = c("BD vs Control" = unname(DX_COLORS["Bipolar"]),
              "SZ vs Control" = unname(DX_COLORS["Schizophrenic"])
             ),
  geom = list(x_min = 0.10, x_range = c(0.5, 3.3), x_or = 4.6, x_max = 16, x_p = 15.6,
              lab_x = 0.12, row_offset = 1.35, breaks = c(0.5, 1, 2, 3)),
  title = "Model 1: group contrasts on error type (relative to repetition errors)",
  header_size = 6.2
)
ggsave(file.path(OUTPUT_DIR, "Figure4_model1_forest.png"), fig4$plot,
       width = 12, height = 7.2, dpi = 300, bg = "white")

# ---- Figure 5: Model 2 (diagnosis and attention) --------------------------------
m2_blocks <- c(
  "SZ vs BD at mean covariates"      = "Diagnosis (SZ vs BD, mean covariates)",
  "attention: slope in BD"           = "Attention \u2014 slope in BD",
  "attention: slope in SZ"           = "Attention \u2014 slope in SZ",
  "attention: SZ - BD (interaction)" = "Attention \u2014 SZ \u2212 BD (interaction)"
)
# Colour by role: slope in BD = green, slope in SZ = yellow, differences = grey.
m2_colours <- setNames(
  ifelse(grepl("slope in BD", names(m2_blocks)), DX_COLORS["Bipolar"],
         ifelse(grepl("slope in SZ", names(m2_blocks)), DX_COLORS["Schizophrenic"], DIFF_COLOR)),
  names(m2_blocks))

fig5 <- make_forest_plot(
  sheet   = "model2_slopes_contrasts",
  blocks  = m2_blocks,
  colours = m2_colours,
  geom = list(x_min = 0.08, x_range = c(0.25, 4), x_or = 5.2, x_max = 18, x_p = 17.6,
              lab_x = 0.09, row_offset = 1.9, breaks = c(0.25, 0.5, 1, 2, 4)),
  title = "Model 2: contrasts on error type (relative to repetition errors)",
  header_size = 5.8
)
ggsave(file.path(OUTPUT_DIR, "Figure5_model2_forest.png"), fig5$plot,
       width = 12, height = fig5$height, dpi = 300, bg = "white")
