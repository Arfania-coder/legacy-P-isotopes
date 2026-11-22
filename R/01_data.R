################################################################################
# 01_data.R
# Read the data, check the design, and calculate the equilibrium d18O_P and the
# departure from equilibrium for every observation (Supplementary Table S1).
#
# Creates:
#   dat      48 observations: 4 fractions x 2 positions x 2 tillage systems x 3 dates
#   samples  12 composite samples (one per plot and date), fractions in columns
#   soil_T   soil-temperature summary of the Pullman 2 NW station (Table S2)
################################################################################

cat("\n== 01  Data and equilibrium reference ==\n")

dat <- read.csv(file.path(DATA_DIR, "d18Op_hedley_fractions.csv"), stringsAsFactors = FALSE)
soil_T <- read.csv(file.path(DATA_DIR, "soil_temperature_summary.csv"), stringsAsFactors = FALSE)

dat$fraction <- factor(dat$fraction, levels = FRACTIONS)
dat$position <- factor(dat$position, levels = c("Top", "Toe"))
dat$tillage  <- factor(dat$tillage,  levels = c("CT", "NT"))
dat$month    <- factor(dat$month,    levels = c("Dec", "Jan", "Feb"))
dat$sampling_date <- as.Date(dat$sampling_date)

# ---- Design checks ------------------------------------------------------------
stopifnot(
  nrow(dat) == 48,
  !anyNA(dat[c("P_mg_kg", "d18Op", "d18Ow_soil_water", "soil_temp_C")]),
  all(table(dat$fraction, dat$position, dat$tillage, dat$month) == 1)
)

# ---- Equilibrium (Eq. 1) and departure from equilibrium -----------------------
# Each observation uses the soil-water d18O of its composite sample and the
# 10-cm soil temperature measured in the field on the sampling date.
dat$d18Op_eq <- eq_chang_blake(dat$d18Ow_soil_water, dat$soil_temp_C)
dat$offset   <- dat$d18Op - dat$d18Op_eq          # Delta-d18O_P in the paper

table_s1 <- data.frame(
  tillage = dat$tillage, fraction = dat$fraction, position = dat$position, month = dat$month,
  d18Op = dat$d18Op, d18Ow_soil_water = dat$d18Ow_soil_water, soil_temp_C = dat$soil_temp_C,
  d18Op_eq = fmt_num(dat$d18Op_eq, 2), offset = fmt_num(dat$offset, 2, signed = TRUE)
)
write_table(table_s1, "TableS1_d18Op_equilibrium.csv")

# ---- Composite samples (n = 12) -----------------------------------------------
# Soil water and temperature belong to the composite sample, not to the fraction.
# In one sample (CT, top slope, January) the soil-water value differs slightly
# among the four fraction rows (-13.42, -13.40, -13.38, -13.38 per mil); the
# sample value is taken as the mean of its rows.
samples <- aggregate(
  cbind(d18Ow_soil_water, soil_temp_C, d18Op_eq) ~ tillage + position + month,
  data = dat, FUN = mean
)
sample_key <- function(d) paste(d$tillage, d$position, d$month)
for (f in FRACTIONS) {
  sub <- dat[dat$fraction == f, ]
  idx <- match(sample_key(samples), sample_key(sub))
  samples[[paste0("P_", f)]]      <- sub$P_mg_kg[idx]
  samples[[paste0("d18Op_", f)]]  <- sub$d18Op[idx]
  samples[[paste0("offset_", f)]] <- sub$offset[idx]
}
samples$P_labile <- rowSums(samples[paste0("P_", LABILE)])
samples$P_total  <- rowSums(samples[paste0("P_", FRACTIONS)])
stopifnot(nrow(samples) == 12, !anyNA(samples))

n_soil_water <- nrow(unique(dat[c("tillage", "position", "month", "d18Ow_soil_water")]))
cat(sprintf("  %d observations, %d composite samples, %d distinct soil-water d18O entries\n",
            nrow(dat), nrow(samples), n_soil_water))

# ---- Values quoted in Section 3.2 and Table 2 ---------------------------------
report("Sect. 3.2", "soil-water d18O, range (per mil)", "-14.5 to -12.6", fmt_range(dat$d18Ow_soil_water))
report("Sect. 3.2", "sampling-date soil temperature, range (C)", "0.2 to 2.5", fmt_range(dat$soil_temp_C))
report("Sect. 3.2 / Table 2", "equilibrium d18O_P, mean +/- SD (n = 12)", "12.1 +/- 0.6",
       fmt_mean_sd(samples$d18Op_eq))
report("Sect. 3.2 / Table 2", "equilibrium d18O_P, range", "11.3 to 13.2", fmt_range(samples$d18Op_eq))
