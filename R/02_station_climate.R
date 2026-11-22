################################################################################
# 02_station_climate.R
# Daily record of the Pullman 2 NW station, 2023-2024: soil-temperature summary
# (Supplementary Table S2), the station values quoted in Sections 2.1 and 2.5,
# and Supplementary Figure S2.
#
# Creates:
#   soil_T          minimum, mean and maximum of the daily mean soil temperature
#                   by period and depth (Table S2, unrounded)
#   month_T_range   range of the daily means at 5-20 cm in each sampling month
#                   (used for the equilibrium bands of Figure 3)
################################################################################

cat("\n== 02  Station climate (Table S2, Figure S2) ==\n")

in_months <- function(months) station$month_number %in% months
periods <- list(
  "December 2023"  = station$year == 2023 & station$month_number == 12,
  "January 2024"   = station$year == 2024 & station$month_number == 1,
  "February 2024"  = station$year == 2024 & station$month_number == 2,
  "Sampling months combined (Dec 2023-Feb 2024)" =
    station$date >= as.Date("2023-12-01") & station$date <= as.Date("2024-02-29"),
  "Winter (Dec-Feb), 2023-2024 pooled"         = in_months(c(12, 1, 2)),
  "Dormant season (Oct-Mar), 2023-2024 pooled" = in_months(c(10, 11, 12, 1, 2, 3)),
  "Spring (Mar-May), 2023-2024 pooled"         = in_months(3:5),
  "Growing season (Apr-Sep), 2023-2024 pooled" = in_months(4:9),
  "Summer (Jun-Aug), 2023-2024 pooled"         = in_months(6:8),
  "Annual 2023" = station$year == 2023,
  "Annual 2024" = station$year == 2024,
  "Full record (Jan 2023-Dec 2024)" = rep(TRUE, nrow(station))
)

soil_T <- do.call(rbind, lapply(names(periods), function(p) {
  do.call(rbind, lapply(names(SOIL_DEPTHS), function(depth) {
    x <- station[[SOIL_DEPTHS[[depth]]]][periods[[p]]]
    data.frame(period = p, depth_cm = as.integer(depth),
               min_C = min(x, na.rm = TRUE), mean_C = mean(x, na.rm = TRUE), max_C = max(x, na.rm = TRUE),
               days_with_data = sum(!is.na(x)), days_total = length(x), stringsAsFactors = FALSE)
  }))
}))
soil_period <- function(p) soil_T[soil_T$period == p, ]

# Range of the daily means over the three depths in each sampling month.
month_period <- c(Dec = "December 2023", Jan = "January 2024", Feb = "February 2024")
month_T_range <- t(vapply(month_period, function(p) {
  r <- soil_period(p)
  c(T_min = min(r$min_C), T_max = max(r$max_C))
}, numeric(2)))

# ---- Table S2 -------------------------------------------------------------------
table_s2 <- do.call(rbind, lapply(names(periods), function(p) {
  r <- soil_period(p)
  cell <- function(depth) {
    k <- r[r$depth_cm == depth, ]
    paste(fmt_num(k$min_C), fmt_num(k$mean_C), fmt_num(k$max_C), sep = " / ")
  }
  data.frame(period = p, cm5 = cell(5), cm10 = cell(10), cm20 = cell(20),
             days = sprintf("%d of %d", r$days_with_data[1], r$days_total[1]), stringsAsFactors = FALSE)
}))
names(table_s2) <- c("Period", "5 cm: min / mean / max (C)", "10 cm: min / mean / max (C)",
                     "20 cm: min / mean / max (C)", "Days with data (of total)")
write_table(table_s2, "TableS2_soil_temperature.csv")
write_table(soil_T, "TableS2_soil_temperature_unrounded.csv")

table_s2_reported <- rbind(
  c("December 2023", "-0.7 / 0.4 / 2.6", "0.0 / 0.9 / 2.6", "1.0 / 1.6 / 2.4", "9 of 31"),
  c("January 2024", "-2.1 / 0.6 / 5.5", "-1.1 / 0.9 / 5.1", "-0.2 / 1.3 / 4.4", "31 of 31"),
  c("February 2024", "0.2 / 3.0 / 5.7", "0.8 / 3.3 / 5.6", "1.5 / 3.5 / 5.2", "29 of 29"),
  c("Sampling months combined (Dec 2023-Feb 2024)", "-2.1 / 1.6 / 5.7", "-1.1 / 1.9 / 5.6", "-0.2 / 2.2 / 5.2", "69 of 91"),
  c("Winter (Dec-Feb), 2023-2024 pooled", "-2.3 / 1.2 / 5.7", "-1.1 / 1.5 / 5.6", "-0.2 / 1.9 / 5.2", "149 of 181"),
  c("Dormant season (Oct-Mar), 2023-2024 pooled", "-2.3 / 3.9 / 15.9", "-1.1 / 4.2 / 15.4", "-0.2 / 4.8 / 15.6", "333 of 365"),
  c("Spring (Mar-May), 2023-2024 pooled", "-0.7 / 9.9 / 24.0", "-0.3 / 9.5 / 22.2", "0.2 / 9.1 / 20.1", "184 of 184"),
  c("Growing season (Apr-Sep), 2023-2024 pooled", "2.7 / 20.1 / 33.1", "3.1 / 19.3 / 30.2", "3.4 / 18.6 / 27.9", "366 of 366"),
  c("Summer (Jun-Aug), 2023-2024 pooled", "12.8 / 25.4 / 33.1", "13.2 / 24.2 / 30.2", "13.6 / 23.0 / 27.9", "184 of 184"),
  c("Annual 2023", "-2.3 / 12.9 / 32.8", "-0.9 / 12.7 / 30.2", "-0.2 / 12.4 / 27.9", "333 of 365"),
  c("Annual 2024", "-2.1 / 11.9 / 33.1", "-1.1 / 11.6 / 29.8", "-0.2 / 11.7 / 27.7", "366 of 366"),
  c("Full record (Jan 2023-Dec 2024)", "-2.3 / 12.4 / 33.1", "-1.1 / 12.1 / 30.2", "-0.2 / 12.0 / 27.9", "699 of 731")
)
stopifnot(identical(table_s2_reported[, 1], table_s2$Period))
column_names <- c("5 cm", "10 cm", "20 cm", "days with data")
for (i in seq_len(nrow(table_s2))) {
  for (j in 1:4) {
    note <- if (i == 1 && j == 3) {
      "The largest December value at 20 cm is exactly 2.45 C (36.41 F on 30 December 2023); halves are rounded up here."
    } else ""
    report("Table S2", paste0(table_s2$Period[i], " | ", column_names[j]),
           table_s2_reported[i, j + 1], table_s2[i, j + 1], note)
  }
}

# ---- Station values quoted in the text ------------------------------------------
annual_precip <- tapply(station$precip_mm, station$year, sum)
wet_season <- in_months(c(11, 12, 1, 2, 3, 4))
wet_share <- 100 * tapply(station$precip_mm[wet_season], station$year[wet_season], sum) / annual_precip
report("Sect. 2.1", "annual precipitation, 2023 and 2024 (mm)", "369, 534",
       paste(fmt_num(annual_precip, 0), collapse = ", "))
report("Sect. 2.1", "share of annual precipitation in November-April, 2023 and 2024 (%)", "70, 70",
       paste(fmt_num(wet_share, 0), collapse = ", "))

soil_columns <- unname(SOIL_DEPTHS)
all_soil <- unlist(station[soil_columns])
sampling_months <- periods[["Sampling months combined (Dec 2023-Feb 2024)"]]
sm <- soil_period("Sampling months combined (Dec 2023-Feb 2024)")
gs <- soil_period("Growing season (Apr-Sep), 2023-2024 pooled")
report("Sect. 2.1 / 2.5", "daily mean soil temperature at 5-20 cm, two-year range (C)", "-2.3 to 33.1",
       fmt_range(all_soil[!is.na(all_soil)]))
report("Sect. 2.1 / 2.5", "daily mean soil temperature at 5-20 cm, sampling months (C)", "-2.1 to 5.7",
       paste(fmt_num(min(sm$min_C)), "to", fmt_num(max(sm$max_C))))
report("Sect. 2.5", "sampling-month mean soil temperature at 5, 10 and 20 cm (C)", "1.6, 1.9, 2.2",
       paste(fmt_num(sm$mean_C), collapse = ", "))
cold_days <- station$soil_temp_5cm_C[sampling_months]
report("Sect. 2.5 / Fig. S2", "days with the 5-cm soil temperature below 0 C in the sampling months", "20 of 69",
       sprintf("%d of %d", sum(cold_days < 0, na.rm = TRUE), sum(!is.na(cold_days))))
report("Sect. 2.5 / Fig. S2", "growing-season mean soil temperature, range over depths (C)", "18.6 to 20.1",
       fmt_range(gs$mean_C))
on_date <- function(d) unlist(station[station$date == as.Date(d), soil_columns])
report("Table S2 caption", "station soil temperature on 18 January 2024 at 5, 10, 20 cm (C)", "-1.5, -0.9, -0.2",
       paste(fmt_num(on_date("2024-01-18")), collapse = ", "))
report("Table S2 caption", "station soil temperature on 22 February 2024 at 5, 10, 20 cm (C)", "4.6, 4.4, 3.9",
       paste(fmt_num(on_date("2024-02-22")), collapse = ", "))
report("Table S2 caption", "station soil-temperature record on 15 December 2023", "none",
       if (all(is.na(on_date("2023-12-15")))) "none" else "present")
missing_days <- station$date[is.na(station$soil_temp_10cm_C)]
gap <- function(d) paste(fmt_date(range(d)), collapse = " to ")
report("Table S2 caption", "days without a soil-temperature record",
       "1 Jan 2023 to 10 Jan 2023; 2 Dec 2023 to 23 Dec 2023",
       paste(gap(missing_days[missing_days < as.Date("2023-07-01")]),
             gap(missing_days[missing_days > as.Date("2023-07-01")]), sep = "; "))

# ---- Figure S2 --------------------------------------------------------------------
depth_labels <- c("5" = "Soil, 5 cm", "10" = "Soil, 10 cm", "20" = "Soil, 20 cm")
soil_long <- do.call(rbind, lapply(names(SOIL_DEPTHS), function(depth) {
  data.frame(date = station$date, temp_C = station[[SOIL_DEPTHS[[depth]]]],
             depth = factor(depth_labels[[depth]], levels = unname(depth_labels)))
}))
sampling_dates  <- sort(unique(dat$sampling_date))
sampling_window <- data.frame(xmin = as.Date("2023-12-01"), xmax = as.Date("2024-02-29"))
quarter_starts  <- seq(as.Date("2023-01-01"), as.Date("2025-01-01"), by = "3 months")

figure_s2 <- ggplot() +
  geom_rect(data = sampling_window,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = "Sampling months")) +
  geom_ribbon(data = station,
              aes(x = date, ymin = air_tmin_C, ymax = air_tmax_C, fill = "Daily air temperature range")) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.25, linetype = "dotted") +
  geom_line(data = soil_long, aes(x = date, y = temp_C, colour = depth), linewidth = 0.4, na.rm = TRUE) +
  geom_vline(xintercept = sampling_dates, colour = INK, linewidth = 0.3, linetype = "dashed") +
  annotate("text", x = as.Date("2024-03-08"), y = 37.5, label = "Sampling dates", hjust = 0, size = 2.8,
           colour = INK) +
  scale_x_date(breaks = quarter_starts, labels = function(d) paste0(month.abb[as.integer(format(d, "%m"))], "\n",
                                                                     format(d, "%Y")),
               expand = c(0.005, 0)) +
  scale_y_continuous(breaks = seq(-30, 40, 10), labels = label_minus) +
  scale_colour_manual(values = c("Soil, 5 cm" = "#86B6EF", "Soil, 10 cm" = "#2A78D6", "Soil, 20 cm" = "#0D366B"),
                      name = NULL) +
  scale_fill_manual(values = c("Daily air temperature range" = "grey82", "Sampling months" = "#F6EBD2"),
                    name = NULL) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2, override.aes = list(linewidth = 0.9))) +
  labs(x = NULL, y = expression("Temperature (" * degree * "C)")) +
  theme_paper() +
  theme(legend.position = "bottom", legend.box = "horizontal",
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25))

save_figure("FigureS2_station_temperature", figure_s2, height_mm = 95)
