#!/usr/bin/env Rscript
################################################################################
# run_all.R
# Reproduces the statistics, tables and data figures of
#
#   Arfania H, Kayler Z, Strawn D, Brooks E, Feber G, Laan M.
#   Oxygen isotopes in phosphate (d18O_P) reveal functional separation between
#   labile and stable legacy phosphorus in soil. Biogeochemistry (in revision).
#
# Run from the repository root:
#   Rscript run_all.R            (terminal)
#   source("run_all.R")          (R console, working directory = repository root)
################################################################################

if (!file.exists(file.path("data", "d18Op_hedley_fractions.csv"))) {
  stop("Run this script from the repository root (the folder that contains data/ and R/).",
       call. = FALSE)
}

scripts <- c(
  "00_setup.R",                    # packages, constants, helper functions
  "01_data.R",                     # data, equilibrium (Eq. 1), Table S1
  "02_station_climate.R",          # Table S2, Figure S2
  "03_fractionation.R",            # Table 1, Section 3.1
  "04_water_isotopes.R",           # Figure 2, Section 3.2
  "05_isotopes.R",                 # Table 2, Sections 3.3-3.4
  "06_figure3.R",                  # Figure 3
  "07_pca.R",                      # Figure 4
  "08_connectivity.R",             # Figure 5, Section 3.5
  "09_equilibrium_sensitivity.R",  # Text S1, Table S3
  "10_seasonal_equilibrium.R"      # Text S2, Figure S1
)
for (script in scripts) source(file.path("R", script))

# ---- Comparison with the values quoted in the manuscript ----------------------
check <- reported_values()
write_table(check, "manuscript_values_check.csv")
cat(sprintf("\n== Values quoted in the manuscript: %d of %d reproduced ==\n",
            sum(check$match), nrow(check)))
if (any(!check$match)) {
  cat("Not reproduced to the printed precision:\n")
  for (i in which(!check$match)) {
    cat(sprintf("  [%s] %s\n      manuscript: %s | computed: %s\n", check$location[i], check$item[i],
                check$manuscript[i], check$computed[i]))
    if (nzchar(check$note[i])) cat("      note:", check$note[i], "\n")
  }
}

writeLines(capture.output(sessionInfo()), file.path("output", "sessionInfo.txt"))
