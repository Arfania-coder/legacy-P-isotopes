################################################################################
# 03_isotopes.R
# Oxygen isotope composition of the four fractions and their departure from
# equilibrium: Table 2 and the statistics of Sections 3.3, 3.4 and 3.5.
#
# Model: fraction + position + tillage + sampling date (additive, type II),
# fitted to measured d18O_P and to the departure from equilibrium (offset).
################################################################################

cat("\n== 03  Isotopic composition of the fractions (Table 2, Sections 3.3-3.5) ==\n")

labile_rows <- dat$fraction %in% LABILE

# ---- Summary by fraction and one-sample t-tests -------------------------------
by_fraction <- do.call(rbind, lapply(FRACTIONS, function(f) {
  d <- dat[dat$fraction == f, ]
  tt <- t.test(d$offset, mu = 0)
  data.frame(
    fraction = f, n = nrow(d),
    d18Op_mean = mean(d$d18Op), d18Op_sd = sd(d$d18Op),
    d18Op_min = min(d$d18Op), d18Op_max = max(d$d18Op),
    offset_mean = mean(d$offset), offset_sd = sd(d$offset),
    offset_min = min(d$offset), offset_max = max(d$offset),
    n_above = sum(d$offset > 0), n_below = sum(d$offset < 0),
    t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value
  )
}))

# ---- Linear models ------------------------------------------------------------
fit_offset <- lm(offset ~ fraction + position + tillage + month, data = dat)
fit_d18Op  <- lm(d18Op  ~ fraction + position + tillage + month, data = dat)
anova_offset <- anova_type2(fit_offset)
anova_d18Op  <- anova_type2(fit_d18Op)

interaction_terms <- c("fraction:position", "fraction:tillage", "fraction:month")
interactions <- rbind(
  cbind(response = "offset", do.call(rbind, lapply(interaction_terms, test_added_term, fit = fit_offset))),
  cbind(response = "d18Op",  do.call(rbind, lapply(interaction_terms, test_added_term, fit = fit_d18Op)))
)

# ---- Tukey-adjusted comparisons of estimated marginal means --------------------
tukey_fraction <- as.data.frame(pairs(emmeans::emmeans(fit_d18Op, ~ fraction), adjust = "tukey"))
tukey_month    <- as.data.frame(pairs(emmeans::emmeans(fit_offset, ~ month), adjust = "tukey"))

p_matrix <- matrix(1, 4, 4, dimnames = list(FRACTIONS, FRACTIONS))
for (i in seq_len(nrow(tukey_fraction))) {
  pair <- strsplit(as.character(tukey_fraction$contrast[i]), " - ")[[1]]
  p_matrix[pair[1], pair[2]] <- p_matrix[pair[2], pair[1]] <- tukey_fraction$p.value[i]
}
by_fraction$letter <- letter_display(setNames(by_fraction$d18Op_mean, FRACTIONS), p_matrix)
tukey_pair <- function(a, b) {
  tukey_fraction[tukey_fraction$contrast %in% c(paste(a, "-", b), paste(b, "-", a)), ]
}

# ---- HNO3-P analysed alone ----------------------------------------------------
anova_HNO3 <- anova_type2(lm(d18Op ~ position + tillage + month, data = dat[dat$fraction == "HNO3", ]))

# ---- Within-month equilibrium band (Figure 3) ----------------------------------
# Range of daily mean soil temperature at 5-20 cm in each sampling month (Table S2),
# combined with the soil-water d18O of each sample.
month_period <- c(Dec = "December 2023", Jan = "January 2024", Feb = "February 2024")
band_T <- t(vapply(month_period, function(p) {
  r <- soil_T[soil_T$period == p, ]
  c(T_min = min(r$min_C), T_max = max(r$max_C))
}, numeric(2)))
dat$band_low  <- eq_chang_blake(dat$d18Ow_soil_water, band_T[as.character(dat$month), "T_max"])
dat$band_high <- eq_chang_blake(dat$d18Ow_soil_water, band_T[as.character(dat$month), "T_min"])

# ---- Pairwise isotopic separation (Figure 5c) ----------------------------------
# Root mean square of the paired per-sample differences in d18O_P between fractions.
# Equilibrium is a property of the composite sample and cancels in the difference
# between two fractions of that sample, so this is also the separation in the
# departure from equilibrium.
d18Op_wide <- as.matrix(samples[paste0("d18Op_", FRACTIONS)])
colnames(d18Op_wide) <- FRACTIONS
separation <- matrix(0, 4, 4, dimnames = list(FRACTIONS, FRACTIONS))
for (i in FRACTIONS) for (j in FRACTIONS) {
  separation[i, j] <- sqrt(mean((d18Op_wide[, i] - d18Op_wide[, j])^2))
}
sep_labile <- separation[LABILE, LABILE][upper.tri(diag(3))]
sep_HNO3   <- separation["HNO3", LABILE]

# ---- Means by sampling date ---------------------------------------------------
by_month <- aggregate(cbind(d18Op, d18Op_eq, offset) ~ fraction + month, data = dat, FUN = mean)
eq_month <- tapply(samples$d18Op_eq, samples$month, mean)

# ---- Table 2 and supporting tables --------------------------------------------
table2 <- data.frame(
  fraction = paste0(FRACTIONS, "-P"),
  d18Op = paste(fmt_num(by_fraction$d18Op_mean), "+/-", fmt_num(by_fraction$d18Op_sd)),
  letter = by_fraction$letter,
  d18Op_range = paste(fmt_num(by_fraction$d18Op_min), "to", fmt_num(by_fraction$d18Op_max)),
  offset = paste(fmt_num(by_fraction$offset_mean, 1, TRUE), "+/-", fmt_num(by_fraction$offset_sd)),
  offset_range = paste(fmt_num(by_fraction$offset_min, 1, TRUE), "to", fmt_num(by_fraction$offset_max, 1, TRUE)),
  relative_to_equilibrium = ifelse(by_fraction$offset_mean > 0,
                                   sprintf("%d of %d above", by_fraction$n_above, by_fraction$n),
                                   sprintf("%d of %d below", by_fraction$n_below, by_fraction$n)),
  t11_p = paste0(fmt_num(by_fraction$t), " (", fmt_p(by_fraction$p), ")"),
  stringsAsFactors = FALSE
)
print(table2, row.names = FALSE, right = FALSE)
write_table(table2, "Table2_d18Op_by_fraction.csv")
write_table(by_fraction, "Table2_d18Op_by_fraction_unrounded.csv")
write_table(rbind(cbind(response = "offset", anova_offset), cbind(response = "d18Op", anova_d18Op),
                  cbind(response = "d18Op of HNO3-P alone", anova_HNO3)),
            "Table2_anova_type2.csv")
write_table(interactions, "Table2_interaction_tests.csv")
write_table(tukey_fraction, "Table2_tukey_fractions_d18Op.csv")
write_table(tukey_month, "Sect3.4_tukey_dates_offset.csv")
write_table(by_month, "Sect3.4_means_by_fraction_and_date.csv")
write_table(data.frame(fraction = FRACTIONS, round(separation, 3), check.names = FALSE),
            "Figure5c_pairwise_separation.csv")
write_table(dat[c("tillage", "position", "month", "fraction", "d18Op", "d18Op_eq", "band_low", "band_high")],
            "Figure3_equilibrium_bands.csv")

# ---- Values quoted in the manuscript ------------------------------------------
table2_reported <- data.frame(
  d18Op = c("15.0 +/- 1.2", "15.7 +/- 1.1", "16.5 +/- 1.3", "9.7 +/- 0.8"),
  d18Op_range = c("12.5 to 16.3", "13.7 to 17.2", "14.2 to 18.4", "7.7 to 10.9"),
  offset = c("+3.0 +/- 1.4", "+3.6 +/- 1.2", "+4.4 +/- 1.2", "-2.4 +/- 0.9"),
  offset_range = c("+0.3 to +4.9", "+2.1 to +5.7", "+2.5 to +6.6", "-3.8 to -1.5"),
  relative_to_equilibrium = c("12 of 12 above", "12 of 12 above", "12 of 12 above", "12 of 12 below"),
  t11_p = c("7.4 (< 0.001)", "10.6 (< 0.001)", "12.4 (< 0.001)", "-9.4 (< 0.001)"),
  stringsAsFactors = FALSE
)
for (i in 1:4) {
  for (col in names(table2_reported)) {
    note <- if (i == 2 && col == "offset_range") {
      "The largest NaHCO3-P offset is 5.649; Table S1 prints it as 5.65."
    } else ""
    report("Table 2", paste0(table2$fraction[i], " | ", col), table2_reported[[col]][i],
           table2[[col]][i], note)
  }
}
report("Table 2", "letters for H2O-P, NaHCO3-P, NaOH-P, HNO3-P", "c, bc, b, a",
       paste(by_fraction$letter, collapse = ", "),
       note = "Same grouping. The script gives 'a' to the greatest mean, as in Table 1.")

term_row <- function(a, term) a[a$term == term, ]
caption <- list(
  c("offset: fraction",  "F(3,40) = 101.8, p < 0.001", fmt_test(term_row(anova_offset, "fraction"))),
  c("offset: position",  "F(1,40) = 1.4, p = 0.25",    fmt_test(term_row(anova_offset, "position"))),
  c("offset: tillage",   "F(1,40) = 2.4, p = 0.13",    fmt_test(term_row(anova_offset, "tillage"))),
  c("offset: date",      "F(2,40) = 5.6, p = 0.007",   fmt_test(term_row(anova_offset, "month"))),
  c("d18O_P: fraction",  "F(3,40) = 119.3, p < 0.001", fmt_test(term_row(anova_d18Op, "fraction"))),
  c("d18O_P: position",  "F(1,40) = 3.1, p = 0.088",   fmt_test(term_row(anova_d18Op, "position"))),
  c("d18O_P: tillage",   "F(1,40) = 0.7, p = 0.41",    fmt_test(term_row(anova_d18Op, "tillage"))),
  c("d18O_P: date",      "F(2,40) = 5.5, p = 0.008",   fmt_test(term_row(anova_d18Op, "month")))
)
for (x in caption) report("Table 2 caption / Sect. 3.4", x[1], x[2], x[3])
eta <- function(a, term, d = 2) fmt_num(a$eta2[a$term == term], d)
report("Table 2 caption", "offset: eta2 of fraction", "0.85", eta(anova_offset, "fraction"))
report("Table 2 caption", "d18O_P: share of variance explained by fraction", "0.87", eta(anova_d18Op, "fraction"))
report("Table 2 caption", "offset: share of variance explained by date", "0.03", eta(anova_offset, "month"))
report("Sect. 4.3", "offset: eta2 of position", "0.004", eta(anova_offset, "position", 3))
report("Sect. 4.3", "offset: eta2 of tillage",  "0.007", eta(anova_offset, "tillage", 3))
report("Table 2 caption", "smallest p among fraction x factor interactions", "0.11", fmt_p(min(interactions$p)))
int_row <- function(resp, term) interactions[interactions$response == resp & interactions$term == term, ]
report("Sect. 3.4", "offset: fraction x position", "F(3,37) = 1.3, p = 0.28", fmt_test(int_row("offset", "fraction:position")))
report("Sect. 3.4 / 4.3", "offset: fraction x tillage", "F(3,37) = 1.8, p = 0.16", fmt_test(int_row("offset", "fraction:tillage")))

report("Sect. 3.3", "labile fractions pooled, mean +/- SD (n = 36)", "15.7 +/- 1.3", fmt_mean_sd(dat$d18Op[labile_rows]))
report("Sect. 3.3", "pooled labile mean minus HNO3-P mean (per mil)", "6.0",
       fmt_num(mean(dat$d18Op[labile_rows]) - mean(dat$d18Op[!labile_rows])))
report("Sect. 3.3", "labile observations above the sampling-date equilibrium", "36 of 36",
       sprintf("%d of %d", sum(dat$offset[labile_rows] > 0), sum(labile_rows)))
report("Sect. 3.3", "smallest labile offset (per mil)", "+0.3", fmt_num(min(dat$offset[labile_rows]), 1, TRUE))
report("Sect. 3.3", "labile observations above the within-month band", "35 of 36",
       sprintf("%d of %d", sum(dat$d18Op[labile_rows] > dat$band_high[labile_rows]), sum(labile_rows)))
report("Sect. 3.3", "HNO3-P observations below the within-month band", "12 of 12",
       sprintf("%d of %d", sum(dat$d18Op[!labile_rows] < dat$band_low[!labile_rows]), sum(!labile_rows)))
report("Sect. 3.2", "within-month equilibrium band (per mil)", "11.1 to 13.6",
       paste(fmt_num(min(dat$band_low)), "to", fmt_num(max(dat$band_high))),
       note = "Lower limit is 11.04 with the monthly temperature ranges of Table S2; quoted as approximate.")
report("Sect. 3.3", "SD of d18O_P, labile fractions (per mil)", "1.1 to 1.3", fmt_range(by_fraction$d18Op_sd[1:3]))
report("Sect. 3.3", "HNO3-P alone: position", "F(1,7) = 0.02, p = 0.90", fmt_test(term_row(anova_HNO3, "position")))
report("Sect. 3.3", "HNO3-P alone: tillage",  "F(1,7) = 0.08, p = 0.78", fmt_test(term_row(anova_HNO3, "tillage")))
report("Sect. 3.3", "HNO3-P alone: date",     "F(2,7) = 1.0, p = 0.42",  fmt_test(term_row(anova_HNO3, "month")))

d_HNO3 <- abs(vapply(LABILE, function(f) tukey_pair(f, "HNO3")$estimate, numeric(1)))
tukey_note <- paste("The manuscript value is from a one-way Tukey HSD by fraction (44 df);",
                    "the comparison of estimated marginal means described in Section 2.6 gives the computed value.")
report("Sect. 3.3", "Tukey: HNO3-P versus labile fractions, differences (per mil)", "5.3 to 6.8", fmt_range(d_HNO3))
report("Sect. 3.3", "Tukey: H2O-P versus NaOH-P, difference (per mil)", "1.5",
       fmt_num(abs(tukey_pair("H2O", "NaOH")$estimate)))
report("Sect. 3.3", "Tukey: H2O-P versus NaOH-P, p", "0.011", fmt_p(tukey_pair("H2O", "NaOH")$p.value), tukey_note)
report("Sect. 3.3", "Tukey: NaHCO3-P versus H2O-P or NaOH-P, smallest p", "0.29",
       fmt_p(min(tukey_pair("H2O", "NaHCO3")$p.value, tukey_pair("NaHCO3", "NaOH")$p.value)), tukey_note)

report("Sect. 3.4", "equilibrium d18O_P, January minus December (per mil)", "0.9",
       fmt_num(eq_month[["Jan"]] - eq_month[["Dec"]]))
jan_shift <- vapply(c("H2O", "NaHCO3"), function(f) {
  m <- by_month[by_month$fraction == f, ]
  m$d18Op[m$month != "Jan"] - m$d18Op[m$month == "Jan"]
}, numeric(2))
report("Sect. 3.4", "H2O-P and NaHCO3-P: December and February minus January (per mil)", "0.9 to 1.6",
       fmt_range(jan_shift))

report("Sect. 3.5 / Fig. 5c", "pairwise separation within the labile domain (per mil)", "1.29 to 1.99",
       fmt_range(sep_labile, 2))
report("Sect. 3.5 / Fig. 5c", "within the labile domain, mean +/- SD", "1.71 +/- 0.37", fmt_mean_sd(sep_labile, 2))
report("Sect. 3.5 / Fig. 5c", "HNO3-P versus labile fractions (per mil)", "5.53 to 6.93", fmt_range(sep_HNO3, 2),
       note = "6.924 with the d18O_P values of Table S1 (two decimals); 6.926 with three-decimal values.")
report("Sect. 3.5 / Fig. 5c", "HNO3-P versus labile fractions, mean +/- SD", "6.21 +/- 0.70", fmt_mean_sd(sep_HNO3, 2))
report("Sect. 3.5", "ratio of the two mean separations", "3.6", fmt_num(mean(sep_HNO3) / mean(sep_labile)))
