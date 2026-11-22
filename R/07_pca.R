################################################################################
# 07_pca.R
# Principal component analysis of the 48 observations (Section 3.4, Figure 4).
#
# Six standardized variables: measured d18O_P, equilibrium d18O_P (Eq. 1), their
# difference, soil-water d18O and d2H, and soil temperature. Equilibrium and
# departure from equilibrium are derived from the other variables, so the
# analysis describes the variation in the data; it is not used to test effects.
################################################################################

cat("\n== 07  Principal component analysis (Figure 4, Section 3.4) ==\n")

pca_variables <- c(d18Op = "d18Op", d18Op_eq = "d18Op_eq", offset = "offset",
                   d18Ow = "d18Ow_soil_water", d2Hw = "d2Hw_soil_water", soil_T = "soil_temp_C")
pca <- prcomp(dat[unname(pca_variables)], center = TRUE, scale. = TRUE)
rownames(pca$rotation) <- names(pca_variables)

# The sign of a principal component is arbitrary. Orient PC1 so that soil-water
# d18O loads negatively and PC2 so that the departure from equilibrium loads
# positively, as in the paper.
flip <- c(-sign(pca$rotation["d18Ow", 1]), sign(pca$rotation["offset", 2]))
pca$rotation[, 1:2] <- sweep(pca$rotation[, 1:2], 2, flip, "*")
pca$x[, 1:2]        <- sweep(pca$x[, 1:2], 2, flip, "*")

variance_pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)
loadings <- sweep(pca$rotation, 2, pca$sdev, "*")        # correlations between variables and components
scores <- data.frame(dat[c("tillage", "position", "month", "fraction")], PC1 = pca$x[, 1], PC2 = pca$x[, 2])

# Share of the variance of the scores associated with each design factor.
factors <- c(date = "month", fraction = "fraction", position = "position", tillage = "tillage")
score_shares <- data.frame(
  factor = names(factors),
  PC1_pct = vapply(factors, function(f) 100 * variance_share(scores$PC1, scores[[f]]), numeric(1)),
  PC2_pct = vapply(factors, function(f) 100 * variance_share(scores$PC2, scores[[f]]), numeric(1)),
  row.names = NULL
)
pc2_by_fraction <- data.frame(
  fraction = FRACTIONS,
  mean = as.vector(tapply(scores$PC2, scores$fraction, mean)),
  sd   = as.vector(tapply(scores$PC2, scores$fraction, sd))
)

cat(sprintf("  PC1 %.1f%%, PC2 %.1f%%\n", variance_pct[1], variance_pct[2]))
write_table(data.frame(component = paste0("PC", seq_along(variance_pct)), variance_pct = variance_pct,
                       cumulative_pct = cumsum(variance_pct)), "Figure4_pca_variance.csv")
write_table(data.frame(variable = rownames(loadings), loadings[, 1:2]), "Figure4_pca_loadings.csv")
write_table(scores, "Figure4_pca_scores.csv")
write_table(score_shares, "Sect3.4_pca_score_variance_by_factor.csv")
write_table(pc2_by_fraction, "Sect3.4_pca_PC2_by_fraction.csv")

# ---- Values quoted in the manuscript ------------------------------------------------
pca_note <- paste("The manuscript values were calculated with the equilibrium column of the first submission;",
                  "the values here use Eq. 1, as in the rest of the revised manuscript.")
check_pca <- function(item, manuscript, computed) {
  report("Sect. 3.4 / Fig. 4", item, manuscript, computed, if (identical(manuscript, computed)) "" else pca_note)
}
share <- function(component, f) score_shares[[paste0(component, "_pct")]][score_shares$factor == f]
check_pca("variance explained by PC1 (%)", "40.9", fmt_num(variance_pct[1]))
check_pca("variance explained by PC2 (%)", "34.1", fmt_num(variance_pct[2]))
check_pca("variance explained by PC1 and PC2 (%)", "75.0", fmt_num(sum(variance_pct[1:2])))
check_pca("PC1 loading, soil-water d18O", "-0.86", fmt_num(loadings["d18Ow", 1], 2))
check_pca("PC1 loading, soil-water d2H", "-0.84", fmt_num(loadings["d2Hw", 1], 2))
check_pca("PC1 loading, equilibrium d18O_P", "-0.70", fmt_num(loadings["d18Op_eq", 1], 2))
check_pca("PC1 loading, soil temperature", "-0.59", fmt_num(loadings["soil_T", 1], 2))
check_pca("PC2 loading, measured d18O_P", "0.86", fmt_num(loadings["d18Op", 2], 2))
check_pca("PC2 loading, departure from equilibrium", "0.95", fmt_num(loadings["offset", 2], 2))
check_pca("share of PC1 score variance: sampling date (%)", "73", fmt_num(share("PC1", "date"), 0))
check_pca("share of PC1 score variance: fraction (%)", "4", fmt_num(share("PC1", "fraction"), 0))
check_pca("share of PC2 score variance: fraction (%)", "68", fmt_num(share("PC2", "fraction"), 0))
check_pca("largest share of a component: landscape position (%)", "2.7",
          fmt_num(max(share("PC1", "position"), share("PC2", "position"))))
check_pca("largest share of a component: tillage (%)", "6.4",
          fmt_num(max(share("PC1", "tillage"), share("PC2", "tillage"))))
hno3 <- pc2_by_fraction[pc2_by_fraction$fraction == "HNO3", ]
check_pca("PC2 score of HNO3-P, mean +/- SD", "-1.98 +/- 0.72", paste(fmt_num(hno3$mean, 2), "+/-", fmt_num(hno3$sd, 2)))
check_pca("PC2 score of the labile fractions, range of means", "+0.35 to +0.99",
          fmt_range(pc2_by_fraction$mean[pc2_by_fraction$fraction %in% LABILE], 2, signed = TRUE))

# ---- Figure 4 -------------------------------------------------------------------------
# Ellipse that contains 95% of a bivariate normal distribution with the mean and
# the covariance of a group of points (a data ellipse, not a confidence region
# for the group mean).
data_ellipse <- function(x, y, level = 0.95, n = 120) {
  radius <- sqrt(2 * qf(level, 2, length(x) - 1))
  angle  <- seq(0, 2 * pi, length.out = n)
  xy <- t(c(mean(x), mean(y)) + radius * t(chol(cov(cbind(x, y)))) %*% rbind(cos(angle), sin(angle)))
  data.frame(PC1 = xy[, 1], PC2 = xy[, 2])
}
group_ellipses <- function(group) {
  out <- do.call(rbind, lapply(levels(scores[[group]]), function(g) {
    cbind(level = g, data_ellipse(scores$PC1[scores[[group]] == g], scores$PC2[scores[[group]] == g]))
  }))
  out[[group]] <- factor(out$level, levels = levels(scores[[group]]))
  out
}
ellipses_fraction <- group_ellipses("fraction")
ellipses_position <- group_ellipses("position")

# Loadings are drawn 2.9 times longer than their value, to fill the range of the scores.
# Each label is placed next to the tip of its arrow: dx and dy are the distance from
# the tip, hjust and vjust the point of the label that is put there.
arrow_scale <- 2.9
arrows <- merge(
  data.frame(variable = rownames(loadings), PC1 = arrow_scale * loadings[, 1], PC2 = arrow_scale * loadings[, 2]),
  data.frame(
    variable = c("d18Op", "offset", "soil_T", "d2Hw", "d18Ow", "d18Op_eq"),
    label = c("delta^18*O[P]", "Delta*delta^18*O[P]", '"Soil T"', "delta^2*H[W]", "delta^18*O[W]",
              'delta^18*O[P*","*eq]'),
    dx    = c(-0.10, 0.05, -0.15,  0.45, -0.12,  0.00),
    dy    = c( 0.12, 0.12,  0.14, -0.36,  0.02, -0.20),
    hjust = c(1, 0, 1, 0.5, 1, 0.5),
    vjust = c(0, 0, 0, 0.5, 0.5, 1)
  ),
  by = "variable"
)

axis_label <- function(k) sprintf("PC%d (%s%%)", k, fmt_num(variance_pct[k]))
limits_x <- range(ellipses_fraction$PC1, ellipses_position$PC1, scores$PC1) + c(-0.2, 0.2)
limits_y <- range(ellipses_fraction$PC2, ellipses_position$PC2, scores$PC2) + c(-0.2, 0.2)

biplot_axes <- function() {
  list(
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25),
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.25)
  )
}
biplot_arrows <- function() {
  list(
    geom_segment(data = arrows, aes(x = 0, y = 0, xend = PC1, yend = PC2), inherit.aes = FALSE,
                 colour = INK, linewidth = 0.3, arrow = arrow(length = unit(1.3, "mm"), type = "closed")),
    geom_text(data = arrows, aes(x = PC1 + dx, y = PC2 + dy, label = label, hjust = hjust, vjust = vjust),
              inherit.aes = FALSE, parse = TRUE, size = 8 / .pt, colour = INK),
    scale_x_continuous(labels = label_minus),
    scale_y_continuous(labels = label_minus),
    coord_cartesian(xlim = limits_x, ylim = limits_y, expand = FALSE),
    labs(x = axis_label(1), y = axis_label(2)),
    theme_paper(),
    theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0),
          legend.spacing.y = unit(0.5, "mm"))
  )
}

figure4a <- ggplot(scores, aes(x = PC1, y = PC2)) +
  biplot_axes() +
  geom_polygon(data = ellipses_fraction, aes(group = fraction, fill = fraction, colour = fraction), alpha = 0.07,
               linewidth = 0.3, show.legend = FALSE) +
  geom_point(aes(colour = fraction, shape = fraction), size = 1.7) +
  scale_colour_manual(values = FRACTION_COLOURS, labels = fraction_expressions(), name = NULL) +
  scale_fill_manual(values = FRACTION_COLOURS, guide = "none") +
  scale_shape_manual(values = FRACTION_SHAPES, labels = fraction_expressions(), name = NULL) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  biplot_arrows()

position_labels <- c(Top = "Top slope", Toe = "Toe slope")
tillage_labels  <- c(CT = "Conventional tillage", NT = "No-till")
figure4b <- ggplot(scores, aes(x = PC1, y = PC2)) +
  biplot_axes() +
  geom_path(data = ellipses_position, aes(group = position, linetype = position), colour = "grey35",
            linewidth = 0.3) +
  geom_point(aes(shape = tillage, fill = position), size = 1.7, colour = INK, stroke = 0.35) +
  scale_shape_manual(values = c(CT = 24, NT = 21), labels = tillage_labels, name = NULL) +
  scale_fill_manual(values = c(Top = "white", Toe = "grey35"), labels = position_labels, name = NULL) +
  scale_linetype_manual(values = c(Top = "dashed", Toe = "solid"), labels = position_labels, name = NULL) +
  guides(shape = guide_legend(order = 1, override.aes = list(fill = "grey70")),
         fill = guide_legend(order = 2, override.aes = list(shape = 21)),
         linetype = guide_legend(order = 2)) +
  biplot_arrows()

save_figure("Figure4_pca", list(figure4a, figure4b), height_mm = 105, ncol = 2)
