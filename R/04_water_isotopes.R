################################################################################
# 04_water_isotopes.R
# Isotopic composition of precipitation and soil water: local meteoric water
# line, soil-water line and deuterium excess (Section 3.2, Figure 2).
#
# Precipitation: 22 samples. Samples whose deuterium excess (d2H - 8 d18O) lies
# more than 1.5 interquartile ranges outside the quartiles are not used for the
# meteoric water line or shown in Figure 2.
# Soil water: one value per composite sample (n = 12).
################################################################################

cat("\n== 04  Water isotopes (Figure 2, Section 3.2) ==\n")

# ---- Precipitation ---------------------------------------------------------------
quartiles <- quantile(precip$d_excess, c(0.25, 0.75))
fence <- 1.5 * diff(quartiles)
precip$excluded <- precip$d_excess < quartiles[[1]] - fence | precip$d_excess > quartiles[[2]] + fence
precip_used <- precip[!precip$excluded, ]
lmwl <- fit_line(precip_used$d18O, precip_used$d2H)

# ---- Soil water -------------------------------------------------------------------
soil_water <- data.frame(
  sample = samples$sample, tillage = samples$tillage, position = samples$position, month = samples$month,
  d18O = samples$d18Ow_soil_water, d2H = samples$d2Hw_soil_water
)
soil_water$d_excess <- soil_water$d2H - 8 * soil_water$d18O
soil_line <- fit_line(soil_water$d18O, soil_water$d2H)

# ---- Deuterium excess of the two water sources ------------------------------------
water_points <- rbind(
  data.frame(source = "Precipitation", d18O = precip_used$d18O, d2H = precip_used$d2H,
             d_excess = precip_used$d_excess),
  data.frame(source = "Soil water", d18O = soil_water$d18O, d2H = soil_water$d2H,
             d_excess = soil_water$d_excess)
)
water_points$source <- factor(water_points$source, levels = names(WATER_COLOURS))
# Location: Welch t-test and Wilcoxon rank-sum test. Spread: Brown-Forsythe test
# (analysis of variance of the absolute deviations from the group medians).
welch    <- t.test(d_excess ~ source, data = water_points)
wilcoxon <- wilcox.test(d_excess ~ source, data = water_points, exact = FALSE)
spread   <- anova(lm(abs(d_excess - ave(d_excess, source, FUN = median)) ~ source, data = water_points))
d_excess_tests <- data.frame(
  comparison = c("means (Welch t-test)", "location (Wilcoxon rank-sum test)", "spread (Brown-Forsythe test)"),
  statistic = c(unname(welch$statistic), unname(wilcoxon$statistic), spread[["F value"]][1]),
  p = c(welch$p.value, wilcoxon$p.value, spread[["Pr(>F)"]][1])
)

# ---- Effect of the water reference on equilibrium --------------------------------
# Equilibrium of each composite sample recalculated with mean precipitation d18O
# instead of its soil water: with the mean of the samples of Figure 2, and with the
# mean of the samples collected in the three sampling months.
in_sampling_months <- precip_used$sampling_date >= as.Date("2023-12-01") &
  precip_used$sampling_date <= as.Date("2024-02-29")
precip_reference <- list("precipitation samples of Figure 2" = precip_used$d18O,
                         "precipitation of the sampling months" = precip_used$d18O[in_sampling_months])
reference_effect <- do.call(rbind, lapply(names(precip_reference), function(r) {
  lowering <- samples$d18Op_eq - eq_chang_blake(mean(precip_reference[[r]]), samples$soil_temp_C)
  data.frame(water_reference = r, n = length(precip_reference[[r]]), mean_d18O = mean(precip_reference[[r]]),
             lowering_min = min(lowering), lowering_mean = mean(lowering), lowering_max = max(lowering),
             stringsAsFactors = FALSE)
}))

# ---- Tables ------------------------------------------------------------------------
summarise_d_excess <- function(x) {
  data.frame(n = length(x), median = median(x), mean = mean(x), sd = sd(x), min = min(x), max = max(x))
}
water_lines <- rbind(cbind(source = "Precipitation (local meteoric water line)", lmwl),
                     cbind(source = "Soil water", soil_line))
d_excess_summary <- rbind(cbind(source = "Precipitation", summarise_d_excess(precip_used$d_excess)),
                          cbind(source = "Soil water", summarise_d_excess(soil_water$d_excess)))
print(cbind(water_lines[1:2], round(water_lines[3:6], 3)), row.names = FALSE, right = FALSE)
write_table(water_lines, "Sect3.2_water_lines.csv")
print(cbind(d_excess_summary[1:2], round(d_excess_summary[3:7], 2)), row.names = FALSE, right = FALSE)
print(cbind(d_excess_tests[1], round(d_excess_tests[2:3], 3)), row.names = FALSE, right = FALSE)
write_table(d_excess_summary, "Sect3.2_d_excess_summary.csv")
write_table(d_excess_tests, "Sect3.2_d_excess_tests.csv")
write_table(reference_effect, "Sect3.2_precipitation_as_reference.csv")
write_table(precip[c("sampling_date", "d18O", "d2H", "d_excess", "excluded")], "Figure2_precipitation.csv")
write_table(soil_water, "Figure2_soil_water.csv")

# ---- Values quoted in the manuscript ------------------------------------------------
report("Sect. 3.2 / Fig. 2", "local meteoric water line, slope", "8.3", fmt_num(lmwl$slope))
soil_note <- paste("The manuscript values count the January sample of the CT top slope three times, because",
                   "its four fraction rows carry three different soil-water d18O values; there are 12 composite samples.")
report("Sect. 3.2 / Fig. 2", "soil-water line, slope", "2.5", fmt_num(soil_line$slope), soil_note)
report("Sect. 3.2 / Fig. 2", "soil-water line, number of samples", "14", as.character(soil_line$n), soil_note)
report("Sect. 3.2", "largest lowering of equilibrium d18O_P with precipitation as the water reference (per mil)",
       "about 2", paste("about", fmt_num(reference_effect$lowering_max[1], 0)))
season_span <- function(d) paste(fmt_date(range(d)), collapse = " to ")
first_season <- precip$sampling_date < as.Date("2023-09-01")
report("Sect. 2.5", "precipitation sampling periods", "January to June, 2022-2023 and 2023-2024 seasons",
       paste(season_span(precip$sampling_date[first_season]), season_span(precip$sampling_date[!first_season]), sep = "; "),
       note = "The first sample of the data file is from 8 March 2023 and the second season starts on 4 December 2023.")

# ---- Figure 2 ------------------------------------------------------------------------
line_labels <- c(
  sprintf("Precipitation (slope %s, n = %d)", fmt_num(lmwl$slope), lmwl$n),
  sprintf("Soil water (slope %s, n = %d)", fmt_num(soil_line$slope), soil_line$n)
)
figure2a <- ggplot(water_points, aes(x = d18O, y = d2H, colour = source, shape = source)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = WATER_COLOURS, labels = line_labels, name = NULL) +
  scale_shape_manual(values = c(16, 17), labels = line_labels, name = NULL) +
  scale_x_continuous(labels = label_minus) +
  scale_y_continuous(labels = label_minus) +
  labs(x = bquote(delta^18 * O ~ "(" * .(PER_MIL) ~ "VSMOW)"),
       y = bquote(delta^2 * H ~ "(" * .(PER_MIL) ~ "VSMOW)")) +
  theme_paper() +
  legend_inside(0.02, 1)

figure2b <- ggplot(water_points, aes(x = source, y = d_excess)) +
  geom_hline(yintercept = 10, colour = "grey45", linewidth = 0.3, linetype = "dashed") +
  geom_boxplot(aes(fill = source), coef = Inf, width = 0.5, linewidth = 0.35, colour = INK, alpha = 0.25) +
  geom_point(aes(colour = source, shape = source), size = 1.4,
             position = position_jitter(width = 0.13, height = 0, seed = 1)) +
  scale_colour_manual(values = WATER_COLOURS) +
  scale_fill_manual(values = WATER_COLOURS) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = NULL, y = bquote("d-excess (" * .(PER_MIL) * ")")) +
  theme_paper() +
  theme(legend.position = "none")

save_figure("Figure2_water_isotopes", list(figure2a, figure2b), height_mm = 95, ncol = 2, widths = c(1.65, 1))
