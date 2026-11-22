################################################################################
# 08_connectivity.R
# Isotopic connectivity among the four fractions (Section 3.5, Figure 5).
#
# The connectivity index of a pathway (from one fraction to another) is obtained
# from 4000 draws. In each draw
#   1. the 12 composite samples are resampled with replacement;
#   2. the isotopic distance between every pair of fractions is calculated as the
#      root mean square of their paired differences in d18O_P;
#   3. a value is drawn for every pathway from the Beta distribution assigned to
#      its distance class (Eq. 6 of the main text):
#          distance < 2 per mil      Beta(5, 2)
#          2 to < 4 per mil          Beta(3, 4)
#          4 per mil or more         Beta(1, 8)
#   4. the three values leaving a fraction are divided by their sum, so that the
#      pathways from each fraction add up to 1.
# Each pathway is summarised by the mean, the median and the 2.5th and 97.5th
# percentiles of its 4000 draws. The index expresses isotopic similarity; it is
# not a measured flux.
#
# Uses d18Op_wide, separation and by_fraction from 05_isotopes.R.
################################################################################

cat("\n== 08  Isotopic connectivity (Figure 5, Section 3.5) ==\n")

N_DRAWS <- 4000
beta_shape <- function(distance) {
  if (distance < 2) c(5, 2) else if (distance < 4) c(3, 4) else c(1, 8)
}

# The random-number generator and its seed are fixed, so the draws are the same
# in every run.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(2026)

connectivity <- array(NA_real_, dim = c(N_DRAWS, 4, 4), dimnames = list(NULL, FRACTIONS, FRACTIONS))
distance     <- array(0, dim = c(N_DRAWS, 4, 4), dimnames = list(NULL, FRACTIONS, FRACTIONS))
for (b in seq_len(N_DRAWS)) {
  resampled <- d18Op_wide[sample(nrow(d18Op_wide), replace = TRUE), , drop = FALSE]
  for (i in FRACTIONS) {
    others <- setdiff(FRACTIONS, i)
    theta <- vapply(others, function(j) {
      distance[b, i, j] <<- sqrt(mean((resampled[, i] - resampled[, j])^2))
      shape <- beta_shape(distance[b, i, j])
      rbeta(1, shape[1], shape[2])
    }, numeric(1))
    connectivity[b, i, others] <- theta / sum(theta)
  }
}

# ---- Summary of the 12 pathways --------------------------------------------------
pathways <- expand.grid(to = FRACTIONS, from = FRACTIONS, stringsAsFactors = FALSE)[c("from", "to")]
pathways <- pathways[pathways$from != pathways$to, ]
pathway_summary <- do.call(rbind, lapply(seq_len(nrow(pathways)), function(k) {
  p <- connectivity[, pathways$from[k], pathways$to[k]]
  data.frame(from = pathways$from[k], to = pathways$to[k], mean = mean(p), median = median(p),
             lower = quantile(p, 0.025, names = FALSE), upper = quantile(p, 0.975, names = FALSE),
             stringsAsFactors = FALSE)
}))
pathway <- function(from, to) pathway_summary[pathway_summary$from == from & pathway_summary$to == to, ]

# Share of the resamples in which a pair of fractions falls in each distance class.
pair_index <- which(upper.tri(diag(4)), arr.ind = TRUE)
distance_classes <- do.call(rbind, lapply(seq_len(nrow(pair_index)), function(k) {
  i <- FRACTIONS[pair_index[k, 1]]
  j <- FRACTIONS[pair_index[k, 2]]
  d <- distance[, i, j]
  data.frame(fraction_1 = i, fraction_2 = j, separation_12_samples = separation[i, j],
             share_below_2 = mean(d < 2), share_2_to_4 = mean(d >= 2 & d < 4), share_4_or_more = mean(d >= 4),
             stringsAsFactors = FALSE)
}))

shown <- data.frame(pathway = paste(pathway_summary$from, "->", pathway_summary$to),
                    "mean (95% interval), %" = fmt_pct_interval(pathway_summary$mean, pathway_summary$lower,
                                                                 pathway_summary$upper),
                    check.names = FALSE)
print(shown, row.names = FALSE, right = FALSE)
write_table(pathway_summary, "Figure5_connectivity_summary.csv")
write_table(distance_classes, "Figure5_distance_classes.csv")
draws_long <- do.call(rbind, lapply(seq_len(nrow(pathways)), function(k) {
  data.frame(draw = seq_len(N_DRAWS), from = pathways$from[k], to = pathways$to[k],
             prob = round(connectivity[, pathways$from[k], pathways$to[k]], 6), stringsAsFactors = FALSE)
}))
draws_long$order <- rep(seq_len(nrow(pathways)), each = N_DRAWS)
draws_long <- draws_long[order(draws_long$draw, draws_long$order), c("draw", "from", "to", "prob")]
write_table(draws_long, "Figure5_connectivity_draws.csv")

# ---- Values quoted in the manuscript ------------------------------------------------
pct <- function(from, to) {
  r <- pathway(from, to)
  fmt_pct_interval(r$mean, r$lower, r$upper)
}
from_HNO3 <- as.vector(connectivity[, "HNO3", LABILE])
connectivity_note <- "Figure 5 of the manuscript shows the values computed here; the values in the text of Section 3.5 differ from the figure."
report("Sect. 3.5", "NaHCO3-P to H2O-P: mean (95% interval), %", "48 (25-72)", pct("NaHCO3", "H2O"), connectivity_note)
report("Sect. 3.5", "H2O-P to NaHCO3-P: mean (95% interval), %", "46 (22-70)", pct("H2O", "NaHCO3"), connectivity_note)
report("Sect. 3.5", "NaOH-P to HNO3-P: mean (95% interval), %", "26 (8-48)", pct("NaOH", "HNO3"), connectivity_note)
report("Sect. 3.5", "HNO3-P to a labile fraction: mean (95% interval), %", "20 (4-42)",
       fmt_pct_interval(mean(from_HNO3), quantile(from_HNO3, 0.025), quantile(from_HNO3, 0.975)),
       paste(connectivity_note, "HNO3-P is 4 per mil or more from every labile fraction in every resample, so its",
             "three pathways are drawn from the same distribution and average one third after division by their sum."))

# ---- Figure 5 -------------------------------------------------------------------------
pathway_label <- function(from, to) paste(FRACTION_LABELS[from], "%->%", FRACTION_LABELS[to])
tint <- function(colour, amount = 0.2) {          # the colour mixed with white
  m <- 255 - amount * (255 - col2rgb(colour))
  rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
}

# (a) Network. Coordinates are in millimetres at final size. The nodes follow the
#     extraction order clockwise from the top left; each arrow lies to the right of
#     its direction of travel, and its width is proportional to the mean probability.
node_radius <- 8.5
nodes <- data.frame(fraction = FRACTIONS, x = c(20, 60, 60, 20), y = c(62, 62, 22, 22),
                    offset = by_fraction$offset_mean, stringsAsFactors = FALSE)
arrow_width <- function(p) 0.3 + 3.0 * p
arrow_polygon <- function(x0, y0, x1, y1, width) {
  len <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
  ux <- (x1 - x0) / len
  uy <- (y1 - y0) / len
  head_width  <- width + 1.6
  head_length <- 0.95 * head_width
  bx <- x1 - head_length * ux
  by <- y1 - head_length * uy
  half <- c(width, width, head_width, 0, -head_width, -width, -width) / 2
  along_x <- c(x0, bx, bx, x1, bx, bx, x0)
  along_y <- c(y0, by, by, y1, by, by, y0)
  data.frame(x = along_x - uy * half, y = along_y + ux * half)
}
network_arrows <- do.call(rbind, lapply(seq_len(nrow(pathway_summary)), function(k) {
  a <- nodes[nodes$fraction == pathway_summary$from[k], ]
  b <- nodes[nodes$fraction == pathway_summary$to[k], ]
  len <- sqrt((b$x - a$x)^2 + (b$y - a$y)^2)
  ux <- (b$x - a$x) / len
  uy <- (b$y - a$y) / len
  shift <- 2.0                                           # to the right of the direction of travel
  clear <- sqrt(node_radius^2 - shift^2) + 1.0           # leaves 1 mm between arrow and node
  cbind(pathway = k, from = pathway_summary$from[k],
        arrow_polygon(a$x + clear * ux + shift * uy, a$y + clear * uy - shift * ux,
                      b$x - clear * ux + shift * uy, b$y - clear * uy - shift * ux,
                      arrow_width(pathway_summary$mean[k])))
}))
circle <- function(x, y, r, n = 120) {
  angle <- seq(0, 2 * pi, length.out = n)
  data.frame(x = x + r * cos(angle), y = y + r * sin(angle))
}
node_circles <- do.call(rbind, lapply(seq_len(nrow(nodes)), function(k) {
  cbind(fraction = nodes$fraction[k], circle(nodes$x[k], nodes$y[k], node_radius))
}))
key <- data.frame(p = c(0.1, 0.3, 0.5), x0 = c(33, 49, 65))
key_arrows <- do.call(rbind, lapply(seq_len(nrow(key)), function(k) {
  cbind(id = k, arrow_polygon(key$x0[k], 3, key$x0[k] + 7, 3, arrow_width(key$p[k])))
}))

figure5a <- ggplot() +
  geom_polygon(data = network_arrows, aes(x = x, y = y, group = pathway, fill = from),
               colour = "white", linewidth = 0.25) +
  geom_polygon(data = node_circles, aes(x = x, y = y, group = fraction, colour = fraction),
               fill = rep(tint(FRACTION_COLOURS[nodes$fraction]), each = 120), linewidth = 0.9) +
  geom_text(data = nodes, aes(x = x, y = y + 2.1, label = FRACTION_LABELS[fraction]), parse = TRUE,
            size = 8 / .pt, colour = INK) +
  geom_text(data = nodes, aes(x = x, y = y - 2.4, label = paste0(with_minus(fmt_num(offset, 1, signed = TRUE)), PER_MIL)),
            size = 8 / .pt, colour = INK) +
  geom_polygon(data = key_arrows, aes(x = x, y = y, group = id), fill = "grey35") +
  annotate("text", x = 2, y = 3, label = "Mean probability", hjust = 0, size = 8 / .pt, colour = INK) +
  annotate("text", x = key$x0 + 8, y = 3, label = paste0(100 * key$p, "%"), hjust = 0, size = 8 / .pt,
           colour = INK) +
  scale_fill_manual(values = FRACTION_COLOURS, guide = "none") +
  scale_colour_manual(values = FRACTION_COLOURS, guide = "none") +
  coord_fixed(xlim = c(0, 80), ylim = c(-1, 74), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(6, 2, 2, 2, unit = "mm"))

# (b) Distribution of the draws for the four pathways that follow the extraction
#     order. The sample is reflected at 0 and 1 so that the density estimate does
#     not fall off artificially at the limits of a probability.
density_unit <- function(x, n = 256) {
  d <- density(c(x, -x, 2 - x), bw = bw.nrd0(x), from = 0, to = 1, n = n)
  data.frame(prob = d$x, density = 3 * d$y)
}
sequence_pairs <- data.frame(from = FRACTIONS, to = c(FRACTIONS[-1], FRACTIONS[1]), row = 4:1,
                             stringsAsFactors = FALSE)
ridges <- do.call(rbind, lapply(seq_len(nrow(sequence_pairs)), function(k) {
  cbind(sequence_pairs[k, ], density_unit(connectivity[, sequence_pairs$from[k], sequence_pairs$to[k]]),
        row.names = NULL)
}))
ridge_scale <- 0.92 / max(ridges$density)   # the tallest ridge almost fills its row
ridges$height <- ridges$row + ridge_scale * ridges$density
ridge_marks <- do.call(rbind, lapply(seq_len(nrow(sequence_pairs)), function(k) {
  r <- ridges[ridges$row == sequence_pairs$row[k], ]
  s <- pathway(sequence_pairs$from[k], sequence_pairs$to[k])
  at <- c(s$lower, s$median, s$upper)
  data.frame(row = sequence_pairs$row[k], prob = at, top = approx(r$prob, r$height, at)$y,
             mark = c("limit", "median", "limit"))
}))

figure5b <- ggplot(ridges, aes(x = 100 * prob)) +
  geom_ribbon(aes(ymin = row, ymax = height, group = row, fill = from), alpha = 0.75) +
  geom_line(aes(y = height, group = row), colour = INK, linewidth = 0.3) +
  geom_segment(data = ridge_marks, aes(x = 100 * prob, xend = 100 * prob, y = row, yend = top, linewidth = mark),
               colour = INK) +
  scale_fill_manual(values = FRACTION_COLOURS, guide = "none") +
  scale_linewidth_manual(values = c(limit = 0.3, median = 0.6), guide = "none") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  scale_y_continuous(breaks = sequence_pairs$row,
                     labels = parse(text = pathway_label(sequence_pairs$from, sequence_pairs$to)),
                     limits = c(1, max(ridges$height) + 0.08), expand = c(0, 0)) +
  labs(x = "Probability (%)", y = NULL) +
  theme_paper() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(vjust = 0, margin = margin(r = 2, unit = "mm")),
        plot.margin = margin(8, 6, 2, 2, unit = "mm"))

# (c) Isotopic separation of each pair of fractions (all 12 composite samples).
pair_cells <- data.frame(
  row = factor(distance_classes$fraction_2, levels = rev(FRACTIONS[-1])),
  column = factor(distance_classes$fraction_1, levels = FRACTIONS[-4]),
  separation = distance_classes$separation_12_samples
)
figure5c <- ggplot(pair_cells, aes(x = column, y = row, fill = separation)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = fmt_num(separation, 2), colour = separation > 4), size = 8.5 / .pt) +
  scale_fill_gradient(low = "#DCEAFB", high = "#0D366B", limits = c(0, 7), breaks = seq(0, 6, 2),
                      name = bquote("Separation (" * .(PER_MIL) * ")")) +
  scale_colour_manual(values = c("FALSE" = INK, "TRUE" = "white"), guide = "none") +
  scale_x_discrete(labels = fraction_expressions(FRACTIONS[-4]), position = "bottom") +
  scale_y_discrete(labels = fraction_expressions(rev(FRACTIONS[-1]))) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_paper() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.key.width = unit(6.4, "mm"), legend.key.height = unit(2.5, "mm"),
        plot.margin = margin(8, 4, 2, 4, unit = "mm"))

# (d) All pathways, in decreasing order of the mean. Colour and symbol show the
#     fraction the pathway starts from, as in the other panels.
ranked <- pathway_summary[order(pathway_summary$mean), ]
ranked$label <- factor(seq_len(nrow(ranked)))
ranked$from  <- factor(ranked$from, levels = FRACTIONS)
figure5d <- ggplot(ranked, aes(y = label, colour = from)) +
  geom_vline(xintercept = 50, colour = "grey45", linewidth = 0.3, linetype = "dashed") +
  geom_segment(aes(x = 100 * lower, xend = 100 * upper, yend = label), linewidth = 0.6) +
  geom_point(aes(x = 100 * mean, shape = from), size = 2.2) +
  scale_colour_manual(values = FRACTION_COLOURS, guide = "none") +
  scale_shape_manual(values = FRACTION_SHAPES, guide = "none") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  scale_y_discrete(labels = parse(text = pathway_label(as.character(ranked$from), ranked$to))) +
  labs(x = "Probability (%): mean and 95% interval", y = NULL) +
  theme_paper() +
  theme(plot.margin = margin(8, 6, 2, 2, unit = "mm"))

save_figure("Figure5_connectivity", list(figure5a, figure5b, figure5c, figure5d), height_mm = 172,
            nrow = 2, ncol = 2, widths = c(1, 1.12), heights = c(1, 1.05))
