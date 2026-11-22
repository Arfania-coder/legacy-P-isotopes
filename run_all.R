#!/usr/bin/env Rscript
################################################################################
# run_all.R
# Reproduces the statistics, tables and Figure 3 of
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

source(file.path("R", "00_setup.R"))
source(file.path("R", "01_data.R"))
source(file.path("R", "02_fractionation.R"))
source(file.path("R", "03_isotopes.R"))
source(file.path("R", "04_equilibrium_sensitivity.R"))
source(file.path("R", "05_figure3.R"))

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
