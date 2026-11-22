################################################################################
# 00_setup.R
# Paths, packages, constants and helper functions shared by all scripts.
#
# Arfania H, Kayler Z, Strawn D, Brooks E, Feber G, Laan M. Oxygen isotopes in
# phosphate (d18O_P) reveal functional separation between labile and stable
# legacy phosphorus in soil. Biogeochemistry (in revision).
################################################################################

# ---- Packages ----------------------------------------------------------------
required_packages <- c("emmeans")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install the missing package(s) first: install.packages(c(",
       paste0('"', missing_packages, '"', collapse = ", "), "))", call. = FALSE)
}

# ---- Paths (relative to the repository root) ----------------------------------
DATA_DIR  <- "data"
TABLE_DIR <- file.path("output", "tables")
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Constants ----------------------------------------------------------------
FRACTIONS <- c("H2O", "NaHCO3", "NaOH", "HNO3")   # extraction order
LABILE    <- c("H2O", "NaHCO3", "NaOH")           # labile legacy domain
ALPHA     <- 0.05                                 # criterion for a difference

# ---- Equilibrium d18O_P -------------------------------------------------------
# Eq. 1 of the main text: Chang and Blake (2015), dissolved phosphate-water,
# calibrated between 3 and 37 degrees C.
#   1000 ln(alpha) = 14.43 * 10^3 / T - 26.54        (T in kelvin)
eq_chang_blake <- function(d18Ow, temp_C) {
  temp_K <- temp_C + 273.15
  exp((14.43e3 / temp_K - 26.54) / 1000) * (d18Ow + 1000) - 1000
}

# Apatite-water calibrations of the form T (deg C) = A - B * (d18O_P - d18O_W),
# used only for the sensitivity analysis (Supplementary Text S1, Table S3).
LINEAR_CALIBRATIONS <- list(
  "Longinelli and Nuti (1973)" = c(A = 111.4, B = 4.30),
  "Lecuyer et al. (2013)"      = c(A = 117.4, B = 4.50),
  "Puceat et al. (2010)"       = c(A = 118.7, B = 4.22)
)
eq_linear <- function(d18Ow, temp_C, A, B) d18Ow + (A - temp_C) / B

# ---- Statistics helpers -------------------------------------------------------
# Type II analysis of variance for an additive linear model: each factor is
# tested after all the others (drop1). eta2 is the factor sum of squares divided
# by the total sum of squares.
anova_type2 <- function(fit) {
  d <- drop1(fit, test = "F")
  y <- model.response(model.frame(fit))
  data.frame(
    term     = rownames(d)[-1],
    df       = d[["Df"]][-1],
    df_resid = df.residual(fit),
    sum_sq   = d[["Sum of Sq"]][-1],
    F        = d[["F value"]][-1],
    p        = d[["Pr(>F)"]][-1],
    eta2     = d[["Sum of Sq"]][-1] / sum((y - mean(y))^2),
    row.names = NULL
  )
}

# F test for one term added to a fitted model (used for interactions).
test_added_term <- function(fit, term) {
  bigger <- update(fit, as.formula(paste(". ~ . +", term)))
  a <- anova(fit, bigger)
  data.frame(term = term, df = a[["Df"]][2], df_resid = a[["Res.Df"]][2],
             F = a[["F"]][2], p = a[["Pr(>F)"]][2])
}

# Compact letter display from a matrix of pairwise p-values. Groups are taken in
# decreasing order of their mean, so "a" marks the greatest mean; groups that
# share a letter are not different at the chosen alpha. (Sets are built from
# neighbours in that order, which is sufficient for the comparisons made here.)
letter_display <- function(means, p_matrix, alpha = ALPHA) {
  ord <- names(sort(means, decreasing = TRUE))
  different <- p_matrix[ord, ord] < alpha
  diag(different) <- FALSE
  sets <- list()
  for (i in seq_along(ord)) {
    j <- i
    while (j < length(ord) && !any(different[i:(j + 1), i:(j + 1)])) j <- j + 1
    candidate <- i:j
    if (!any(vapply(sets, function(s) all(candidate %in% s), logical(1)))) {
      sets[[length(sets) + 1]] <- candidate
    }
  }
  out <- setNames(character(length(ord)), ord)
  for (k in seq_along(sets)) {
    out[ord[sets[[k]]]] <- paste0(out[ord[sets[[k]]]], letters[k])
  }
  out[names(means)]
}

# ---- Formatting helpers -------------------------------------------------------
# Numbers are rounded half away from zero (106.35 -> 106.4), the convention used
# in the tables of the paper; sprintf() alone rounds the binary representation.
fmt_num <- function(x, digits = 1, signed = FALSE) {
  m <- 10^digits
  r <- sign(x) * floor(abs(x) * m + 0.5 + 1e-9) / m
  out <- formatC(abs(r), format = "f", digits = digits)
  zero <- as.numeric(out) == 0
  paste0(ifelse(r < 0 & !zero, "-", ifelse(signed, "+", "")), out)
}
fmt_range <- function(x, digits = 1, signed = FALSE) {
  paste(fmt_num(min(x), digits, signed), "to", fmt_num(max(x), digits, signed))
}
fmt_mean_sd <- function(x, digits = 1, signed = FALSE) {
  paste(fmt_num(mean(x), digits, signed), "+/-", fmt_num(sd(x), digits))
}
fmt_p <- function(p) {
  ifelse(p < 0.001, "< 0.001", ifelse(p < 0.1, fmt_num(p, 3), fmt_num(p, 2)))
}
# F statistics: one decimal; two decimals below 0.1.
fmt_F <- function(f, digits = NULL) {
  if (!is.null(digits)) return(fmt_num(f, digits))
  ifelse(f < 0.1, fmt_num(f, 2), fmt_num(f, 1))
}
# Number of decimals printed in a quoted value such as "0.05".
decimals_of <- function(token) {
  ifelse(grepl(".", token, fixed = TRUE), nchar(sub("^[^.]*\\.", "", token)), 0)
}
fmt_test <- function(row, F_digits = NULL) {
  sprintf("F(%d,%d) = %s, p %s", row$df, row$df_resid, fmt_F(row$F, F_digits),
          ifelse(row$p < 0.001, "< 0.001", paste("=", fmt_p(row$p))))
}

# ---- Record of the values quoted in the manuscript ----------------------------
# Every statistic quoted in the paper is registered with report(); run_all.R
# writes the collected table (output/tables/manuscript_values_check.csv) so that
# script output and manuscript can be compared line by line.
.reported <- new.env()
.reported$rows <- list()
report <- function(location, item, manuscript, computed, note = "") {
  computed <- as.character(computed)
  .reported$rows[[length(.reported$rows) + 1]] <- data.frame(
    location = location, item = item, manuscript = manuscript, computed = computed,
    match = identical(manuscript, computed), note = note, stringsAsFactors = FALSE
  )
  invisible(computed)
}
reported_values <- function() do.call(rbind, .reported$rows)

write_table <- function(x, file) {
  path <- file.path(TABLE_DIR, file)
  write.csv(x, path, row.names = FALSE)
  cat("  wrote", path, "\n")
  invisible(path)
}
