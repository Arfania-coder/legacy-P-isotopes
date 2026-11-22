################################################################################
# 05_figure3.R
# Figure 3: measured d18O_P of the four fractions relative to equilibrium, by
# slope position and tillage.
#
# Points: one per fraction and composite sample, coloured by fraction and shaped
#         by sampling month.
# Boxes:  equilibrium d18O_P (Eq. 1) of each composite sample over the range of
#         daily mean soil temperature at 5-20 cm in its sampling month (Table S2;
#         calculated in 03_isotopes.R). The line inside a box is its midpoint.
#
# The figure is drawn at its final size (174 mm wide) with lettering of 8 pt or
# more and without a title, as the journal asks; the title belongs in the caption.
#
# Output: output/figures/Figure3_d18Op_vs_equilibrium.png  (600 dpi)
#         output/figures/Figure3_d18Op_vs_equilibrium.pdf  (vector)
################################################################################

cat("\n== 05  Figure 3 ==\n")

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop('Figure 3 needs the ggplot2 package: install.packages("ggplot2")', call. = FALSE)
}
library(ggplot2)

FIGURE_DIR <- file.path("output", "figures")
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

# The y-axis label contains a per mille sign. If the session is not running in a
# UTF-8 locale (for example Rscript on a server with LANG=C), switch the character
# type to UTF-8 so that the sign is drawn correctly.
if (!isTRUE(l10n_info()[["UTF-8"]]) && .Platform$OS.type == "unix") {
  for (loc in c("C.UTF-8", "en_US.UTF-8")) {
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  }
}

# ---- Data for the figure ------------------------------------------------------
group_levels <- c("Top CT", "Toe CT", "Top NT", "Toe NT")
fig3 <- dat
fig3$group <- factor(paste(fig3$position, fig3$tillage), levels = group_levels)

# One box per composite sample (12 boxes).
bands <- aggregate(cbind(band_low, band_high) ~ group + month, data = fig3, FUN = mean)
bands$band_mid <- (bands$band_low + bands$band_high) / 2
stopifnot(nrow(bands) == 12)

# ---- Appearance ---------------------------------------------------------------
fraction_colours <- c(H2O = "#5793C7", NaHCO3 = "#9C6FB6", NaOH = "#79C377", HNO3 = "#EE8E6B")
fraction_labels  <- c(expression(H[2] * "O-P"), expression(NaHCO[3] * "-P"),
                      expression("NaOH-P"), expression(HNO[3] * "-P"))
month_shapes <- c(Dec = 16, Jan = 17, Feb = 15)
month_fills  <- c(Dec = "#DBF0EC", Jan = "#9AD3CB", Feb = "#6AAFA6")
month_labels <- c(Dec = "December", Jan = "January", Feb = "February")

# Left-aligned legend labels (the theme element changed name in ggplot2 3.5.0).
legend_text <- if (utils::packageVersion("ggplot2") >= "3.5.0") {
  theme(legend.text = element_text(size = 8, hjust = 0))
} else {
  theme(legend.text = element_text(size = 8), legend.text.align = 0)
}

figure3 <- ggplot() +
  geom_crossbar(
    data = bands,
    aes(x = group, y = band_mid, ymin = band_low, ymax = band_high, fill = month),
    position = position_dodge(width = 0.7), width = 0.7, colour = "grey25", linewidth = 0.3
  ) +
  geom_point(
    data = fig3,
    aes(x = group, y = d18Op, colour = fraction, shape = month, group = fraction),
    position = position_dodge(width = 0.72), size = 2.1
  ) +
  scale_x_discrete(labels = function(x) sub(" ", "\n", x)) +
  scale_y_continuous(breaks = seq(8, 18, 2)) +
  scale_colour_manual(values = fraction_colours, labels = fraction_labels,
                      name = "Hedley P fraction") +
  scale_shape_manual(values = month_shapes, labels = month_labels, name = "Sampling month") +
  scale_fill_manual(values = month_fills, labels = month_labels,
                    name = expression("Equilibrium " * delta^18 * O[P] * " (soil-T range)")) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(shape = 18, size = 3.4)),
    shape  = guide_legend(order = 2, override.aes = list(colour = "grey25", size = 2.4)),
    fill   = guide_legend(order = 3)
  ) +
  labs(x = "Slope position and tillage system",
       y = expression(delta^18 * O[P] ~ "(\u2030 VSMOW)")) +
  theme_classic(base_size = 9) +
  theme(
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 8.5),
    legend.key.height = unit(4.4, "mm"),
    legend.key.width = unit(5, "mm"),
    legend.spacing.y = unit(2, "mm")
  ) +
  legend_text

# ---- Save ---------------------------------------------------------------------
width_mm  <- 174
height_mm <- 112
png_file <- file.path(FIGURE_DIR, "Figure3_d18Op_vs_equilibrium.png")
ggsave(png_file, figure3, width = width_mm, height = height_mm, units = "mm", dpi = 600)
cat("  wrote", png_file, "\n")

pdf_file <- file.path(FIGURE_DIR, "Figure3_d18Op_vs_equilibrium.pdf")
if (capabilities("cairo")) {
  ggsave(pdf_file, figure3, width = width_mm, height = height_mm, units = "mm",
         device = grDevices::cairo_pdf)
  cat("  wrote", pdf_file, "\n")
} else {
  cat("  PDF not written: this R installation has no cairo graphics device\n")
}
