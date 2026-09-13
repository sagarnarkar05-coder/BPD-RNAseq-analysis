library(ggplot2)
library(dplyr)

# -----------------------------------------
# Read CIBERSORTx results and metadata
# -----------------------------------------

cell <- read.delim(
  "../CIBERSORTx_Results.txt",
  check.names = FALSE
)

meta <- read.delim(
  "../Sample_metadata.txt",
  check.names = FALSE
)

# Join sample information with cell fractions
data <- merge(
  meta,
  cell,
  by.x = "Sample",
  by.y = "Mixture"
)

# Correct time-point order
data$Time_point <- factor(
  data$Time_point,
  levels = c("Cord", "Day14", "Day28")
)

# Correct group order
data$Group <- factor(
  data$Group,
  levels = c("nonBPD", "BPD")
)


# -----------------------------------------
# Function to make one trajectory plot
# -----------------------------------------

make_trajectory <- function(cell_type, title, filename) {

  # Mean and standard error
  summary <- data %>%
    group_by(Group, Time_point) %>%
    summarise(
      Mean = mean(.data[[cell_type]], na.rm = TRUE),
      SEM = sd(.data[[cell_type]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  p <- ggplot(
    summary,
    aes(
      x = Time_point,
      y = Mean,
      group = Group
    )
  ) +

    # Lines
    geom_line(
      aes(linetype = Group),
      colour = "black",
      linewidth = 0.8
    ) +

    # Error bars
    geom_errorbar(
      aes(
        ymin = Mean - SEM,
        ymax = Mean + SEM
      ),
      width = 0.07,
      linewidth = 0.6,
      colour = "black"
    ) +

    # Points
    geom_point(
      aes(shape = Group),
      size = 3,
      colour = "black"
    ) +

    # Match your existing style
    scale_shape_manual(
      values = c(
        "nonBPD" = 16,
        "BPD" = 17
      )
    ) +

    scale_linetype_manual(
      values = c(
        "nonBPD" = "solid",
        "BPD" = "dashed"
      )
    ) +

    labs(
      title = title,
      x = NULL,
      y = "Mean estimated fraction",
      shape = "Group",
      linetype = "Group"
    ) +

    # Add space so error bars are never clipped
    scale_y_continuous(
      expand = expansion(
        mult = c(0.08, 0.15)
      )
    ) +

    theme_classic(base_size = 13) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 17
      ),

      axis.title.y = element_text(
        face = "bold"
      ),

      axis.text = element_text(
        colour = "black"
      ),

      legend.title = element_text(
        face = "bold"
      ),

      legend.position = "right"
    )

  # Save directly in CELL/figures
  ggsave(
    filename,
    p,
    width = 5.5,
    height = 4.5,
    dpi = 600,
    bg = "white"
  )

  print(p)
}


# -----------------------------------------
# Make the three plots
# -----------------------------------------

make_trajectory(
  "Neutrophils",
  "Neutrophils trajectory",
  "Neutrophils_trajectory.png"
)

make_trajectory(
  "T cells CD4 naive",
  "T cells CD4 naive trajectory",
  "T_cells_CD4_naive_trajectory.png"
)

make_trajectory(
  "B cells naive",
  "B cells naive trajectory",
  "B_cells_naive_trajectory.png"
)

cat("\nTrajectory plots completed.\n")

