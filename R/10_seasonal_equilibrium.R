################################################################################
# 10_seasonal_equilibrium.R
# Equilibrium d18O_P by season, from precipitation d18O and the soil temperature
# of the station record (Supplementary Text S2, Figure S1).
#
# For each season
#   - water:        d18O of the precipitation samples collected in the months of
#                   the season (all 22 samples of the two collection periods);
#   - temperature:  daily mean soil temperature at 20 cm in those months, 2023 and
#                   2024 pooled;
#   - point:        Eq. 1 at the mean d18O and the mean temperature;
#   - bar:          from Eq. 1 at the lower quartile of d18O and the upper quartile
#                   of temperature, to Eq. 1 at the upper quartile of d18O and the
#                   lower quartile of temperature.
# The seasons overlap (winter is part of the dormant season, and spring shares two
# months with the growing season).
#
# The calculation is repeated without the two precipitation samples that are not
# used for the meteoric water line (04_water_isotopes.R); that result is written
# to the table but not plotted.
#
# Uses periods (02_station_climate.R), precip$excluded (04_water_isotopes.R) and
# by_fraction (05_isotopes.R).
################################################################################

cat("\n== 10  Seasonal equilibrium (Figure S1, Text S2) ==\n")

seasons <- data.frame(
  season = c("Winter", "Dormant", "Spring", "Growing"),
  months = c("Dec-Feb", "Oct-Mar", "Mar-May", "Apr-Sep"),
  period = c("Winter (Dec-Feb), 2023-2024 pooled", "Dormant season (Oct-Mar), 2023-2024 pooled",
             "Spring (Mar-May), 2023-2024 pooled", "Growing season (Apr-Sep), 2023-2024 pooled"),
  stringsAsFactors = FALSE
)
season_months <- list(Winter = c(12, 1, 2), Dormant = c(10, 11, 12, 1, 2, 3), Spring = 3:5, Growing = 4:9)

seasonal_equilibrium <- function(precipitation, samples_used) {
  do.call(rbind, lapply(seq_len(nrow(seasons)), function(k) {
    water <- precipitation$d18O[precipitation$month_number %in% season_months[[seasons$season[k]]]]
    temp  <- station$soil_temp_20cm_C[periods[[seasons$period[k]]]]
    temp  <- temp[!is.na(temp)]
    w <- quantile(water, c(0.25, 0.75), names = FALSE)
    t <- quantile(temp,  c(0.25, 0.75), names = FALSE)
    data.frame(
      precipitation_samples = samples_used, season = seasons$season[k], months = seasons$months[k],
      n_precipitation = length(water), n_days = length(temp),
      d18O_precip_mean = mean(water), d18O_precip_q25 = w[1], d18O_precip_q75 = w[2],
      soil_temp_20cm_mean = mean(temp), soil_temp_20cm_q25 = t[1], soil_temp_20cm_q75 = t[2],
      eq_mean = eq_chang_blake(mean(water), mean(temp)),
      eq_low  = eq_chang_blake(w[1], t[2]),
      eq_high = eq_chang_blake(w[2], t[1]),
      stringsAsFactors = FALSE
    )
  }))
}
seasonal_eq      <- seasonal_equilibrium(precip, "all 22")
seasonal_eq_kept <- seasonal_equilibrium(precip[!precip$excluded, ], "without the 2 samples excluded in Figure 2")
shown <- rbind(seasonal_eq, seasonal_eq_kept)
print(cbind(shown[1:4], round(shown[c("eq_mean", "eq_low", "eq_high")], 2)), row.names = FALSE, right = FALSE)
write_table(shown, "FigureS1_seasonal_equilibrium.csv")

# ---- Values quoted in the manuscript ------------------------------------------------
growing <- seasonal_eq[seasonal_eq$season == "Growing", ]
fraction_mean <- setNames(by_fraction$d18Op_mean, FRACTIONS)
inside <- seasonal_eq$season[seasonal_eq$eq_low <= fraction_mean[["HNO3"]] &
                               seasonal_eq$eq_high >= fraction_mean[["HNO3"]]]
report("Text S2 / Sect. 4.2", "growing-season equilibrium from precipitation, mean (per mil)", "9.0",
       fmt_num(growing$eq_mean))
report("Text S2 / Sect. 4.2", "growing-season equilibrium, range (per mil)", "7.1 to 11.0",
       paste(fmt_num(growing$eq_low), "to", fmt_num(growing$eq_high)),
       note = "The bar of Figure S1 runs from 6.5 to 11.0; its lower end is 6.1 to 6.5 with the temperature at 5, 10 or 20 cm.")
report("Fig. S1 caption", "seasons whose range includes the HNO3-P mean", "Spring, Growing",
       paste(inside, collapse = ", "))
report("Fig. S1 caption", "labile fraction means above the range of every season", "yes",
       if (min(fraction_mean[LABILE]) > max(seasonal_eq$eq_high)) "yes" else "no")

# ---- Figure S1 ------------------------------------------------------------------------
season_labels <- paste0(seasonal_eq$season, "\n(", sub("-", EN_DASH, seasonal_eq$months, fixed = TRUE), ")")
seasonal_eq$label <- factor(season_labels, levels = season_labels)
fraction_lines <- data.frame(fraction = factor(FRACTIONS, levels = FRACTIONS), d18Op = by_fraction$d18Op_mean)

figure_s1 <- ggplot(seasonal_eq, aes(x = label)) +
  geom_segment(data = fraction_lines, aes(x = 0.5, xend = 4.5, y = d18Op, yend = d18Op, colour = fraction),
               linewidth = 0.4, linetype = "dashed") +
  geom_text(data = fraction_lines, aes(x = 4.58, y = d18Op, label = FRACTION_LABELS[as.character(fraction)]),
            parse = TRUE, hjust = 0, size = 8 / .pt, colour = INK) +
  geom_linerange(aes(ymin = eq_low, ymax = eq_high), colour = INK, linewidth = 0.5) +
  geom_point(aes(y = eq_mean), colour = INK, size = 2.2) +
  scale_colour_manual(values = FRACTION_COLOURS, guide = "none") +
  scale_x_discrete(limits = season_labels) +
  scale_y_continuous(breaks = seq(6, 18, 2)) +
  coord_cartesian(xlim = c(0.5, 5.35), ylim = c(6, 17.5), expand = FALSE) +
  labs(x = NULL, y = bquote("Equilibrium or measured" ~ delta^18 * O[P] ~ "(" * .(PER_MIL) ~ "VSMOW)")) +
  theme_paper() +
  theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25))

save_figure("FigureS1_seasonal_equilibrium", figure_s1, height_mm = 95, width_mm = 129)
