################################################################################
# 09_equilibrium_sensitivity.R
# Sensitivity of the equilibrium reference to the phosphate-water calibration and
# to extrapolation below 3 degrees C (Section 2.5, Supplementary Text S1, Table S3).
################################################################################

cat("\n== 09  Sensitivity of the equilibrium reference (Text S1, Table S3) ==\n")

# ---- Equilibrium of every observation under each calibration ------------------
calibrations <- c("Chang and Blake (2015), Eq. 1", names(LINEAR_CALIBRATIONS))
eq_by_calibration <- function(name, d18Ow, temp_C) {
  if (name == "Chang and Blake (2015), Eq. 1") {
    eq_chang_blake(d18Ow, temp_C)
  } else {
    k <- LINEAR_CALIBRATIONS[[name]]
    eq_linear(d18Ow, temp_C, k[["A"]], k[["B"]])
  }
}

table_s3 <- do.call(rbind, lapply(calibrations, function(name) {
  eq  <- eq_by_calibration(name, dat$d18Ow_soil_water, dat$soil_temp_C)
  off <- dat$d18Op - eq
  mean_off <- tapply(off, dat$fraction, mean)
  data.frame(
    calibration = name,
    eq_at_2C_minus14 = eq_by_calibration(name, -14, 2),
    eq_mean_48 = mean(eq),
    offset_H2O = mean_off[["H2O"]], offset_NaHCO3 = mean_off[["NaHCO3"]],
    offset_NaOH = mean_off[["NaOH"]], offset_HNO3 = mean_off[["HNO3"]],
    labile_above = sum(off[dat$fraction %in% LABILE] > 0),
    HNO3_below = sum(off[dat$fraction == "HNO3"] < 0),
    stringsAsFactors = FALSE
  )
}))
offset_cols <- c("offset_H2O", "offset_NaHCO3", "offset_NaOH", "offset_HNO3")
print(cbind(table_s3[1], round(table_s3[2:7], 2), table_s3[8:9]), row.names = FALSE, right = FALSE)
write_table(table_s3, "TableS3_calibration_sensitivity.csv")

# Largest change of a fraction offset relative to the calibration of the main text.
shift_from_eq1 <- sweep(as.matrix(table_s3[-1, offset_cols]), 2, unlist(table_s3[1, offset_cols]))

# ---- Extrapolation below the calibrated range (3-37 degrees C) ----------------
# Equilibrium at the coldest sampling-date temperature (0.2 C) relative to the
# lower limit of the calibration (3.0 C), and the mean slope over that interval.
extrapolation <- eq_chang_blake(-14, 0.2) - eq_chang_blake(-14, 3.0)
slope_0.2_to_3 <- -extrapolation / (3.0 - 0.2)

# ---- Range over the sampling months -------------------------------------------
# Full range of daily mean soil temperature at 5-20 cm (Dec 2023-Feb 2024; Table S2)
# combined with the full range of measured soil-water d18O.
winter <- soil_period("Sampling months combined (Dec 2023-Feb 2024)")
stopifnot(nrow(winter) == 3)
T_range <- c(min(winter$min_C), max(winter$max_C))
w_range <- range(dat$d18Ow_soil_water)
eq_range <- c(eq_chang_blake(w_range[1], T_range[2]), eq_chang_blake(w_range[2], T_range[1]))
labile <- dat[dat$fraction %in% LABILE, ]
inside <- sort(labile$d18Op[labile$d18Op > eq_range[1] & labile$d18Op < eq_range[2]])

write_table(data.frame(
  quantity = c("equilibrium at 0.2 C minus equilibrium at 3.0 C (per mil)",
               "mean slope of Eq. 1 between 0.2 and 3.0 C (per mil per K)",
               "sampling-month soil temperature, minimum (C)", "sampling-month soil temperature, maximum (C)",
               "equilibrium over the sampling months, minimum (per mil)",
               "equilibrium over the sampling months, maximum (per mil)"),
  value = c(extrapolation, slope_0.2_to_3, T_range, eq_range)
), "TextS1_extrapolation_and_range.csv")

# ---- Values quoted in the manuscript ------------------------------------------
s3_reported <- list(
  c("11.9", "12.1", "+3.0", "+3.6", "+4.4", "-2.4", "36 / 12"),
  c("11.4", "11.6", "+3.4", "+4.0", "+4.8", "-2.0", "36 / 12"),
  c("11.6", "11.8", "+3.2", "+3.8", "+4.6", "-2.2", "36 / 12"),
  c("13.7", "13.9", "+1.2", "+1.8", "+2.6", "-4.2", "34 / 12")
)
s3_items <- c("equilibrium at 2 C, d18O_W = -14", "mean equilibrium, 48 obs.", "H2O-P offset",
              "NaHCO3-P offset", "NaOH-P offset", "HNO3-P offset", "labile above / HNO3-P below")
for (i in seq_len(nrow(table_s3))) {
  r <- table_s3[i, ]
  computed <- c(fmt_num(r$eq_at_2C_minus14), fmt_num(r$eq_mean_48),
                fmt_num(unlist(r[offset_cols]), 1, TRUE),
                sprintf("%d / %d", r$labile_above, r$HNO3_below))
  for (j in seq_along(s3_items)) {
    report("Table S3", paste0(r$calibration, " | ", s3_items[j]), s3_reported[[i]][j], computed[j])
  }
}
report("Sect. 2.5 / Text S1", "spread of the four calibrations at 2 C (per mil)", "2.3",
       fmt_num(diff(range(table_s3$eq_at_2C_minus14))),
       note = "13.65 - 11.44 = 2.21; the 2.3 in the text is the difference of the rounded values 13.7 and 11.4.")
report("Text S1", "largest change of a fraction offset relative to Eq. 1 (per mil)", "1.8",
       fmt_num(max(abs(shift_from_eq1))))
report("Sect. 2.5 / Text S1", "mean slope of Eq. 1 between 0.2 and 3.0 C (per mil per K)", "-0.19",
       fmt_num(slope_0.2_to_3, 2))
report("Text S1", "equilibrium at 0.2 C minus equilibrium at 3.0 C (per mil)", "0.54", fmt_num(extrapolation, 2))
report("Sect. 2.5", "sampling-month soil temperature at 5-20 cm (C)", "-2.1 to 5.7", fmt_range(T_range))
report("Text S1", "equilibrium over the sampling-month range (per mil)", "10.7 to 14.1",
       paste(fmt_num(eq_range[1]), "to", fmt_num(eq_range[2])),
       note = "Lower limit is 10.62 with the lowest measured soil-water value (-14.54); 10.7 corresponds to -14.5.")
report("Text S1", "labile observations inside that range", "3 of 36: 12.5, 13.3, 13.7",
       sprintf("%d of %d: %s", length(inside), nrow(labile), paste(fmt_num(inside), collapse = ", ")))
