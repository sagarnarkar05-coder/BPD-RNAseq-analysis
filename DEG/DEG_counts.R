rm(list = ls())
gc()

# Package and file paths
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

results_dir <- "/home/sagar007/RNAseq/DEG/results"
figure_dir <- "/home/sagar007/RNAseq/figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

down_colour <- "#3C5488"
up_colour <- "#E64B35"

# Comparisons
comparisons <- data.frame(
  Section = c(rep("Disease contrasts", 4), rep("Longitudinal contrasts", 6)),
  Comparison = c(
    "Overall", "Cord", "Day 14", "Day 28",
    "BPD: Day 14 vs cord", "BPD: Day 28 vs day 14", "BPD: Day 28 vs cord",
    "nonBPD: Day 14 vs cord", "nonBPD: Day 28 vs day 14", "nonBPD: Day 28 vs cord"
  ),
  Analysis = c(
    "Overall_BPD_vs_nonBPD", "Cord_BPD_vs_nonBPD",
    "Day14_BPD_vs_nonBPD", "Day28_BPD_vs_nonBPD",
    "BPD_Day14_vs_Cord_paired", "BPD_Day28_vs_Day14_paired",
    "BPD_Day28_vs_Cord_paired", "nonBPD_Day14_vs_Cord_paired",
    "nonBPD_Day28_vs_Day14_paired", "nonBPD_Day28_vs_Cord_paired"
  ),
  File = c(
    "DEG_Overall_BPD_vs_nonBPD_significant.csv",
    "DEG_Cord_BPD_vs_nonBPD_significant.csv",
    "DEG_Day14_BPD_vs_nonBPD_significant.csv",
    "DEG_Day28_BPD_vs_nonBPD_significant.csv",
    "DEG_BPD_paired_Day14_vs_Cord_significant.csv",
    "DEG_BPD_paired_Day28_vs_Day14_significant.csv",
    "DEG_BPD_paired_Day28_vs_Cord_significant.csv",
    "DEG_nonBPD_paired_Day14_vs_Cord_significant.csv",
    "DEG_nonBPD_paired_Day28_vs_Day14_significant.csv",
    "DEG_nonBPD_paired_Day28_vs_Cord_significant.csv"
  )
)

# Count upregulated and downregulated genes
count_degs <- function(file) {
  results <- read.csv(file.path(results_dir, file), check.names = FALSE)
  if (!"log2FoldChange" %in% names(results)) stop("Missing log2FoldChange column: ", file)
  c(Up = sum(results$log2FoldChange > 0), Down = sum(results$log2FoldChange < 0))
}

deg_counts <- t(vapply(comparisons$File, count_degs, numeric(2)))
comparisons$Up <- as.integer(deg_counts[, "Up"])
comparisons$Down <- as.integer(deg_counts[, "Down"])
comparisons$Total <- comparisons$Up + comparisons$Down

# Check totals against the DEG summary
summary_file <- file.path(results_dir, "DEG_analysis_summary.csv")
if (file.exists(summary_file)) {
  deg_summary <- read.csv(summary_file)
  expected <- deg_summary$Significant_DEGs[match(comparisons$Analysis, deg_summary$Analysis)]
  if (anyNA(expected) || any(comparisons$Total != expected)) stop("DEG totals do not match DEG_analysis_summary.csv")
}

plot_data <- rbind(
  data.frame(comparisons[c("Section", "Comparison")], Direction = "Downregulated", Count = comparisons$Down),
  data.frame(comparisons[c("Section", "Comparison")], Direction = "Upregulated", Count = comparisons$Up)
)
plot_data$Direction <- factor(plot_data$Direction, levels = c("Downregulated", "Upregulated"))
plot_data$Label <- format(plot_data$Count, big.mark = ",", scientific = FALSE)

# Plot function
make_panel <- function(section, x_limit, legend = FALSE) {
  data <- plot_data[plot_data$Section == section, ]
  order <- comparisons$Comparison[comparisons$Section == section]
  data$Comparison <- factor(data$Comparison, levels = rev(order))

  p <- ggplot(data, aes(Count, Comparison, fill = Direction)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62) +
    geom_text(
      aes(label = Label, colour = Direction),
      position = position_dodge(width = 0.72), hjust = -0.12,
      size = 3.4, fontface = "bold", show.legend = FALSE
    ) +
    scale_fill_manual(values = c(Downregulated = down_colour, Upregulated = up_colour)) +
    scale_colour_manual(values = c(Downregulated = down_colour, Upregulated = up_colour)) +
    scale_x_continuous(
      limits = c(0, x_limit), breaks = pretty(c(0, x_limit), 7),
      labels = function(x) format(x, big.mark = ",", scientific = FALSE),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(title = section, x = "Number of significant DEGs", y = NULL, fill = NULL) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_text(colour = "black"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey88"),
      plot.margin = margin(8, 35, 8, 8),
      legend.position = if (legend) "top" else "none",
      legend.justification = "left"
    )
  p
}

disease_plot <- make_panel("Disease contrasts", 450, TRUE)
longitudinal_plot <- make_panel("Longitudinal contrasts", 2250)

# Combine panels
draw_figure <- function() {
  grid.newpage()
  layout <- grid.layout(3, 1, heights = unit(c(0.5, 3.2, 4.5), "null"))
  pushViewport(viewport(layout = layout))
  grid.text(
    "Significant DEG counts across all comparisons",
    gp = gpar(fontsize = 19, fontface = "bold"),
    vp = viewport(layout.pos.row = 1)
  )
  print(disease_plot, vp = viewport(layout.pos.row = 2))
  print(longitudinal_plot, vp = viewport(layout.pos.row = 3))
  popViewport()
}

# Save outputs
write.csv(comparisons, file.path(figure_dir, "Figure_1E_DEG_counts.csv"), row.names = FALSE)

png(file.path(figure_dir, "Figure_1E_DEG_counts.png"), 3600, 2700, res = 300, bg = "white")
draw_figure()
dev.off()

pdf(file.path(figure_dir, "Figure_1E_DEG_counts.pdf"), 12, 9, useDingbats = FALSE)
draw_figure()
dev.off()

print(comparisons[c("Comparison", "Up", "Down", "Total")], row.names = FALSE)

