################################################################################
# 06_figure3.R
# Figure 3: measured d18O_P of the four fractions relative to equilibrium, by
# slope position and tillage.
#
# Points: one per fraction and composite sample; fractions are told apart by
#         colour and by their position within each group, sampling months by symbol.
# Boxes:  equilibrium d18O_P (Eq. 1) of each composite sample over the range of
#         daily mean soil temperature at 5-20 cm in its sampling month
#         (calculated in 05_isotopes.R). The line inside a box is its midpoint.
################################################################################

cat("\n== 06  Figure 3 ==\n")

group_levels <- c("Top CT", "Toe CT", "Top NT", "Toe NT")
fig3 <- dat
fig3$group <- factor(paste(fig3$position, fig3$tillage), levels = group_levels)

# One box per composite sample (12 boxes).
bands <- aggregate(cbind(band_low, band_high) ~ group + month, data = fig3, FUN = mean)
bands$band_mid <- (bands$band_low + bands$band_high) / 2
stopifnot(nrow(bands) == 12)

month_shapes <- c(Dec = 16, Jan = 17, Feb = 15)
month_fills  <- c(Dec = "#DBF0EC", Jan = "#9AD3CB", Feb = "#6AAFA6")

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
  scale_colour_manual(values = FRACTION_COLOURS, labels = fraction_expressions(),
                      name = "Hedley P fraction") +
  scale_shape_manual(values = month_shapes, labels = MONTH_LABELS, name = "Sampling month") +
  scale_fill_manual(values = month_fills, labels = MONTH_LABELS,
                    name = expression("Equilibrium " * delta^18 * O[P] * " (soil-T range)")) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(shape = 18, size = 3.4)),
    shape  = guide_legend(order = 2, override.aes = list(colour = "grey25", size = 2.4)),
    fill   = guide_legend(order = 3)
  ) +
  labs(x = "Slope position and tillage system",
       y = bquote(delta^18 * O[P] ~ "(" * .(PER_MIL) ~ "VSMOW)")) +
  theme_paper()

save_figure("Figure3_d18Op_vs_equilibrium", figure3, height_mm = 112)
