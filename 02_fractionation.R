################################################################################
# 02_fractionation.R
# Sequential (Hedley) P fractionation: Table 1 and the statistics of Section 3.1.
#
# Responses: concentration of each fraction, the sum of the three labile
# fractions, and total extracted P (n = 12 composite samples per response).
# Model: position + tillage + sampling date (additive, type II); the
# position x tillage interaction is tested by adding it to that model.
################################################################################

cat("\n== 02  Phosphorus fractionation (Table 1, Section 3.1) ==\n")

responses <- c(P_H2O = "H2O-P", P_NaHCO3 = "NaHCO3-P", P_NaOH = "NaOH-P",
               P_HNO3 = "HNO3-P", P_labile = "Labile sum", P_total = "Total P")
digits <- c(P_H2O = 1, P_NaHCO3 = 1, P_NaOH = 1, P_HNO3 = 0, P_labile = 0, P_total = 0)

fit_fraction <- function(y) {
  fit <- lm(reformulate(c("position", "tillage", "month"), response = y), data = samples)
  rbind(
    anova_type2(fit)[c("term", "df", "df_resid", "F", "p")],
    test_added_term(fit, "position:tillage")
  )
}
anova_P <- lapply(names(responses), fit_fraction)
names(anova_P) <- names(responses)
anova_row <- function(y, term) anova_P[[y]][anova_P[[y]]$term == term, ]

# ---- Factor means (n = 6) with letters ----------------------------------------
factor_means <- function(y, factor_name) {
  x <- split(samples[[y]], samples[[factor_name]])
  m <- vapply(x, mean, numeric(1))
  lett <- if (anova_row(y, factor_name)$p < ALPHA) ifelse(m == max(m), "a", "b") else c("a", "a")
  data.frame(response = y, factor = factor_name, level = names(x),
             n = lengths(x), mean = m, sd = vapply(x, sd, numeric(1)), letter = lett,
             row.names = NULL)
}
means_P <- do.call(rbind, lapply(names(responses), function(y) {
  rbind(factor_means(y, "position"), factor_means(y, "tillage"))
}))

# ---- Position x tillage cell means (n = 3) ------------------------------------
cell_means <- do.call(rbind, lapply(names(responses), function(y) {
  a <- aggregate(samples[[y]], samples[c("tillage", "position")],
                 function(v) c(mean = mean(v), sd = sd(v)))
  data.frame(response = y, tillage = a$tillage, position = a$position,
             mean = a$x[, "mean"], sd = a$x[, "sd"])
}))
toe_gain <- function(y, till) {
  m <- cell_means[cell_means$response == y & cell_means$tillage == till, ]
  100 * (m$mean[m$position == "Toe"] / m$mean[m$position == "Top"] - 1)
}

# ---- Table 1 (lower block) ----------------------------------------------------
mean_cell <- function(y, factor_name, level) {
  r <- means_P[means_P$response == y & means_P$factor == factor_name & means_P$level == level, ]
  paste(fmt_num(r$mean, digits[[y]]), "+/-", fmt_num(r$sd, digits[[y]]), r$letter)
}
anova_cell <- function(y, term, F_digits = NULL) {
  r <- anova_row(y, term)
  paste0(fmt_F(r$F, F_digits), ", ", fmt_p(r$p))
}
table1_cols <- c("P_H2O", "P_NaHCO3", "P_NaOH", "P_HNO3", "P_total")
table1 <- data.frame(
  row = c("Position: Top", "Position: Toe", "Tillage: CT", "Tillage: NT",
          "Position (df 1, 7): F, p", "Tillage (df 1, 7): F, p", "Date (df 2, 7): F, p",
          "Position x tillage (df 1, 6): F, p"),
  stringsAsFactors = FALSE
)
for (y in table1_cols) {
  table1[[responses[[y]]]] <- c(
    mean_cell(y, "position", "Top"), mean_cell(y, "position", "Toe"),
    mean_cell(y, "tillage", "CT"), mean_cell(y, "tillage", "NT"),
    anova_cell(y, "position"), anova_cell(y, "tillage"), anova_cell(y, "month"),
    anova_cell(y, "position:tillage")
  )
}
print(table1, row.names = FALSE, right = FALSE)
write_table(table1, "Table1_factor_means_and_anova.csv")
write_table(do.call(rbind, Map(function(y, a) cbind(response = y, a), names(anova_P), anova_P)),
            "Table1_anova_unrounded.csv")
write_table(cell_means, "Sect3.1_position_by_tillage_means.csv")
write_table(samples[c("tillage", "position", "month", paste0("P_", FRACTIONS), "P_labile", "P_total")],
            "Table1_concentrations_by_sample.csv")

# ---- Values quoted in the manuscript ------------------------------------------
table1_reported <- list(
  P_H2O    = c("16.5 +/- 5.0 a", "21.7 +/- 5.4 a", "21.2 +/- 4.5 a", "17.1 +/- 6.3 a",
               "4.6, 0.069", "2.9, 0.13", "2.7, 0.13", "0.05, 0.83"),
  P_NaHCO3 = c("42.5 +/- 8.1 b", "54.6 +/- 8.0 a", "54.4 +/- 5.5 a", "42.8 +/- 10.3 b",
               "15.0, 0.006", "13.7, 0.008", "0.7, 0.54", "0.7, 0.43"),
  P_NaOH   = c("91.0 +/- 13.4 b", "112.7 +/- 1.5 a", "106.4 +/- 7.0 a", "97.4 +/- 19.2 a",
               "17.4, 0.004", "3.0, 0.13", "0.7, 0.55", "5.1, 0.065"),
  P_HNO3   = c("546 +/- 5 a", "566 +/- 23 a", "546 +/- 8 b", "567 +/- 22 a",
               "5.6, 0.050", "6.4, 0.040", "0.3, 0.75", "15.5, 0.008"),
  P_total  = c("697 +/- 21 b", "755 +/- 23 a", "728 +/- 20 a", "724 +/- 51 a",
               "18.2, 0.004", "0.1, 0.78", "0.7, 0.54", "16.8, 0.006")
)
anova_terms <- c("position", "tillage", "month", "position:tillage")
for (y in table1_cols) {
  for (i in seq_len(nrow(table1))) {
    quoted <- table1_reported[[y]][i]
    computed <- if (i <= 4) {
      table1[[responses[[y]]]][i]
    } else {                       # F compared at the number of decimals printed in Table 1
      anova_cell(y, anova_terms[i - 4], decimals_of(sub(",.*$", "", quoted)))
    }
    report("Table 1", paste0(responses[[y]], " | ", table1$row[i]), quoted, computed)
  }
}

share_HNO3   <- 100 * samples$P_HNO3 / samples$P_total
share_labile <- 100 * samples$P_labile / samples$P_total
report("Sect. 3.1", "total P, range (mg/kg)", "667 to 790", fmt_range(samples$P_total, 0))
report("Sect. 3.1", "total P, mean +/- SD (n = 12)", "726 +/- 37", fmt_mean_sd(samples$P_total, 0))
report("Sect. 3.1", "HNO3-P share of total P, range (%)", "73 to 82", fmt_range(share_HNO3, 0))
report("Sect. 3.1", "HNO3-P share of total P, mean +/- SD (%)", "77 +/- 3", fmt_mean_sd(share_HNO3, 0))
report("Sect. 3.1", "labile share of total P, range (%)", "18 to 27", fmt_range(share_labile, 0))

cm <- function(y, till, pos, d = 0) {
  r <- cell_means[cell_means$response == y & cell_means$tillage == till & cell_means$position == pos, ]
  paste(fmt_num(r$mean, d), "+/-", fmt_num(r$sd, d))
}
report("Sect. 3.1", "total P, NT toe", "768 +/- 19", cm("P_total", "NT", "Toe"))
report("Sect. 3.1", "total P, NT top", "679 +/- 11", cm("P_total", "NT", "Top"))
report("Sect. 3.1", "total P, CT toe", "741 +/- 20", cm("P_total", "CT", "Toe"))
report("Sect. 3.1", "total P, CT top", "714 +/- 5",  cm("P_total", "CT", "Top"))
report("Sect. 3.1", "total P, toe-slope gain under NT (%)", "+13", fmt_num(toe_gain("P_total", "NT"), 0, TRUE))
report("Sect. 3.1", "total P, toe-slope gain under CT (%)", "+4",  fmt_num(toe_gain("P_total", "CT"), 0, TRUE))

by_level <- function(y, factor_name, level, d = 0) fmt_mean_sd(samples[[y]][samples[[factor_name]] == level], d)
report("Sect. 3.1", "labile sum, toe", "189 +/- 11", by_level("P_labile", "position", "Toe"))
report("Sect. 3.1", "labile sum, top", "150 +/- 22", by_level("P_labile", "position", "Top"))
report("Sect. 3.1", "labile sum, CT",  "182 +/- 15", by_level("P_labile", "tillage", "CT"))
report("Sect. 3.1", "labile sum, NT",  "157 +/- 30", by_level("P_labile", "tillage", "NT"))
report("Sect. 3.1", "labile sum, position: F(1,7), p", "36.9, < 0.001", anova_cell("P_labile", "position"))
report("Sect. 3.1", "labile sum, tillage: F(1,7), p",  "14.8, 0.006",   anova_cell("P_labile", "tillage"))
report("Sect. 3.1", "labile sum, position x tillage: F(1,6), p", "7.4, 0.035",
       anova_cell("P_labile", "position:tillage"))
report("Sect. 3.1", "labile sum, toe-slope gain under NT (%)", "+39", fmt_num(toe_gain("P_labile", "NT"), 0, TRUE))
report("Sect. 3.1", "labile sum, toe-slope gain under CT (%)", "+16", fmt_num(toe_gain("P_labile", "CT"), 0, TRUE))

report("Sect. 3.1", "HNO3-P, NT toe", "585 +/- 7", cm("P_HNO3", "NT", "Toe"))
other_cells <- cell_means[cell_means$response == "P_HNO3" &
                            !(cell_means$tillage == "NT" & cell_means$position == "Toe"), "mean"]
report("Sect. 3.1", "HNO3-P, other three position-tillage means", "545 to 548", fmt_range(other_cells, 0))
report("Sect. 3.1", "HNO3-P, toe-slope gain under NT (%)", "+7",   fmt_num(toe_gain("P_HNO3", "NT"), 0, TRUE))
report("Sect. 3.1", "HNO3-P, toe-slope gain under CT (%)", "+0.2", fmt_num(toe_gain("P_HNO3", "CT"), 1, TRUE))
report("Sect. 3.1", "labile share of total P under CT (%)", "25.0 +/- 1.5",
       fmt_mean_sd(share_labile[samples$tillage == "CT"]))
report("Sect. 3.1", "labile share of total P under NT (%)", "21.6 +/- 2.7",
       fmt_mean_sd(share_labile[samples$tillage == "NT"]))
p_date <- vapply(names(responses)[1:4], function(y) anova_row(y, "month")$p, numeric(1))
report("Sect. 3.1", "smallest p for a date effect on a fraction", "0.13", fmt_p(min(p_date)))

# ---- Phosphorus stocks in the 0-15 cm layer -----------------------------------
# kg P per ha = mg/kg x 10^-6 x depth (m) x bulk density (kg/m3) x 10^4 m2/ha
stock <- function(mg_kg, bulk_density_g_cm3) mg_kg * 1e-6 * 0.15 * bulk_density_g_cm3 * 1e3 * 1e4
crop_removal <- 12.98                                  # kg P/ha/yr (Welikhe et al. 2023)
report("Sect. 3.1", "NaHCO3-P, mean +/- SD (n = 12)", "48.6 +/- 9.9", fmt_mean_sd(samples$P_NaHCO3))
report("Sect. 3.1", "labile sum, mean +/- SD (n = 12)", "170 +/- 26", fmt_mean_sd(samples$P_labile, 0))
report("Sect. 3.1", "HNO3-P, mean +/- SD (n = 12)", "556 +/- 20", fmt_mean_sd(samples$P_HNO3, 0),
       note = "The SD of the 12 values in Table 1 is 19.0.")
stocks <- data.frame(
  pool = c("NaHCO3-P", "Labile sum", "HNO3-P"),
  mean_mg_kg = c(mean(samples$P_NaHCO3), mean(samples$P_labile), mean(samples$P_HNO3))
)
stocks$kg_ha_bd1.2 <- stock(stocks$mean_mg_kg, 1.2)
stocks$kg_ha_bd1.3 <- stock(stocks$mean_mg_kg, 1.3)
stocks$years_removal_bd1.2 <- stocks$kg_ha_bd1.2 / crop_removal
stocks$years_removal_bd1.3 <- stocks$kg_ha_bd1.3 / crop_removal
write_table(stocks, "Sect3.1_P_stocks.csv")
nearest5 <- function(x) 5 * round(x / 5)
report("Sect. 3.1", "labile stock, 0-15 cm (kg P/ha), nearest 5", "305 to 330",
       paste(nearest5(stocks$kg_ha_bd1.2[2]), "to", nearest5(stocks$kg_ha_bd1.3[2])))
report("Sect. 3.1", "labile stock in years of crop removal", "23 to 26",
       paste(floor(stocks$years_removal_bd1.2[2]), "to", ceiling(stocks$years_removal_bd1.3[2])))
report("Sect. 3.1", "HNO3-P stock, 0-15 cm (kg P/ha), nearest 10", "1000 to 1080",
       paste(10 * round(stocks$kg_ha_bd1.2[3] / 10), "to", 10 * round(stocks$kg_ha_bd1.3[3] / 10)))
