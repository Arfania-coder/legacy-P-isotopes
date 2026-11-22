################################################################################
# 00_setup.R
# Paths, packages, constants and helper functions shared by all scripts.
#
# Arfania H, Kayler Z, Strawn D, Brooks E, Feber G, Laan M. Oxygen isotopes in
# phosphate (d18O_P) reveal functional separation between labile and stable
# legacy phosphorus in soil. Biogeochemistry (in revision).
################################################################################

# ---- Packages ----------------------------------------------------------------
required_packages <- c("emmeans", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install the missing package(s) first: install.packages(c(",
       paste0('"', missing_packages, '"', collapse = ", "), "))", call. = FALSE)
}
library(ggplot2)

# ---- Paths (relative to the repository root) ----------------------------------
DATA_DIR   <- "data"
TABLE_DIR  <- file.path("output", "tables")
FIGURE_DIR <- file.path("output", "figures")
dir.create(TABLE_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Constants ----------------------------------------------------------------
FRACTIONS <- c("H2O", "NaHCO3", "NaOH", "HNO3")   # extraction order
LABILE    <- c("H2O", "NaHCO3", "NaOH")           # labile legacy domain
ALPHA     <- 0.05                                 # criterion for a difference
SOIL_DEPTHS <- c("5" = "soil_temp_5cm_C", "10" = "soil_temp_10cm_C", "20" = "soil_temp_20cm_C")

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

# Share of the variance of y associated with a factor (one-way sum of squares
# divided by the total sum of squares).
variance_share <- function(y, factor_values) {
  a <- anova(lm(y ~ factor(factor_values)))
  a[["Sum Sq"]][1] / sum(a[["Sum Sq"]])
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

# Ordinary least-squares line y = a + b x, with the statistics reported for it.
fit_line <- function(x, y) {
  fit <- lm(y ~ x)
  s <- summary(fit)
  data.frame(n = length(x), slope = coef(fit)[[2]], intercept = coef(fit)[[1]],
             r2 = s$r.squared, p = s$coefficients[2, 4])
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
# Dates as "8 Mar 2023", independent of the language of the R session.
fmt_date <- function(d) {
  d <- as.Date(d)
  paste(as.integer(format(d, "%d")), month.abb[as.integer(format(d, "%m"))], format(d, "%Y"))
}
# Probability as a percentage with its 95% interval: "50 (29-73)".
fmt_pct_interval <- function(mean, lower, upper) {
  sprintf("%s (%s-%s)", fmt_num(100 * mean, 0), fmt_num(100 * lower, 0), fmt_num(100 * upper, 0))
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

# ---- Figures ------------------------------------------------------------------
# One colour and one symbol per fraction in every figure. The colours are from
# the Okabe-Ito set; because two of them are close for red-green colour-blind
# readers, fractions are always distinguished by a second cue as well (symbol,
# position or a direct label).
FRACTION_COLOURS <- c(H2O = "#0072B2", NaHCO3 = "#009E73", NaOH = "#CC79A7", HNO3 = "#D55E00")
FRACTION_SHAPES  <- c(H2O = 16, NaHCO3 = 17, NaOH = 15, HNO3 = 18)
FRACTION_LABELS  <- c(H2O = 'H[2]*"O-P"', NaHCO3 = 'NaHCO[3]*"-P"', NaOH = '"NaOH-P"', HNO3 = 'HNO[3]*"-P"')
fraction_expressions <- function(x = FRACTIONS) parse(text = FRACTION_LABELS[as.character(x)])
MONTH_LABELS <- c(Dec = "December", Jan = "January", Feb = "February")
WATER_COLOURS <- c(Precipitation = "#0072B2", "Soil water" = "#D55E00")
INK <- "grey15"

# The axis labels contain a per mille sign. If the session is not running in a
# UTF-8 locale (for example Rscript on a server with LANG=C), switch the character
# type to UTF-8 so that the sign is drawn correctly.
if (!isTRUE(l10n_info()[["UTF-8"]]) && .Platform$OS.type == "unix") {
  for (loc in c("C.UTF-8", "en_US.UTF-8")) {
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  }
}
PER_MIL <- "\u2030"
MINUS   <- "\u2212"
EN_DASH <- "\u2013"

# Numbers in figures carry a true minus sign instead of a hyphen.
with_minus  <- function(x) sub("-", MINUS, x, fixed = TRUE)
label_minus <- function(x) {
  out <- with_minus(format(x, trim = TRUE, scientific = FALSE))
  out[is.na(x)] <- NA
  out
}

# Theme for all figures: drawn at final size, lettering of 8 pt or more, no
# titles inside the figure (the title belongs in the caption).
theme_paper <- function(base_size = 9) {
  left_aligned <- if (utils::packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.text = element_text(size = 8, hjust = 0))
  } else {
    theme(legend.text = element_text(size = 8), legend.text.align = 0)
  }
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(size = 8, colour = "black"),
      axis.title = element_text(size = 9),
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      legend.title = element_text(size = 8.5),
      legend.key.height = unit(4.2, "mm"),
      legend.key.width = unit(5, "mm"),
      legend.background = element_blank(),
      plot.margin = margin(3, 3, 2, 3, unit = "mm")
    ) +
    left_aligned
}

# Legend inside the plotting area; x and y are fractions of the panel and `just`
# is the corner of the legend placed there. (ggplot2 3.5.0 renamed the setting.)
legend_inside <- function(x, y, just = c(0, 1)) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.position = "inside", legend.position.inside = c(x, y), legend.justification = just)
  } else {
    theme(legend.position = c(x, y), legend.justification = just)
  }
}

# Draw one or more plots on a grid, with bold lower-case panel letters.
draw_panels <- function(plots, nrow = 1, ncol = 1, widths = rep(1, ncol), heights = rep(1, nrow),
                        labels = NULL) {
  grid::grid.newpage()
  layout <- grid::grid.layout(nrow, ncol, widths = grid::unit(widths, "null"),
                              heights = grid::unit(heights, "null"))
  grid::pushViewport(grid::viewport(layout = layout))
  for (i in seq_along(plots)) {
    row <- (i - 1) %/% ncol + 1
    col <- (i - 1) %% ncol + 1
    grid::pushViewport(grid::viewport(layout.pos.row = row, layout.pos.col = col))
    grid::grid.draw(ggplotGrob(plots[[i]]))
    if (!is.null(labels)) {
      grid::grid.text(labels[i], x = grid::unit(1.2, "mm"),
                      y = grid::unit(1, "npc") - grid::unit(1, "mm"), just = c("left", "top"),
                      gp = grid::gpar(fontface = "bold", fontsize = 10))
    }
    grid::popViewport()
  }
  grid::popViewport()
  invisible(NULL)
}

# Save a figure as PNG (600 dpi) and, where a cairo device is available, as PDF.
# `plots` is a ggplot object or a list of them (panels in reading order).
save_figure <- function(name, plots, height_mm, width_mm = 174, nrow = 1, ncol = 1,
                        widths = rep(1, ncol), heights = rep(1, nrow), labels = NULL) {
  if (inherits(plots, "ggplot")) plots <- list(plots)
  if (is.null(labels) && length(plots) > 1) labels <- letters[seq_along(plots)]
  draw <- function() draw_panels(plots, nrow, ncol, widths, heights, labels)

  png_file <- file.path(FIGURE_DIR, paste0(name, ".png"))
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(png_file, width = width_mm, height = height_mm, units = "mm", res = 600)
  } else if (capabilities("cairo")) {
    grDevices::png(png_file, width = width_mm, height = height_mm, units = "mm", res = 600, type = "cairo")
  } else {
    grDevices::png(png_file, width = width_mm, height = height_mm, units = "mm", res = 600)
  }
  draw()
  grDevices::dev.off()
  cat("  wrote", png_file, "\n")

  if (capabilities("cairo")) {
    pdf_file <- file.path(FIGURE_DIR, paste0(name, ".pdf"))
    grDevices::cairo_pdf(pdf_file, width = width_mm / 25.4, height = height_mm / 25.4)
    draw()
    grDevices::dev.off()
    cat("  wrote", pdf_file, "\n")
  }
  invisible(png_file)
}
