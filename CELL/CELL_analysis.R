suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
})

# File paths and settings
project_dir <- "/home/sagar007/RNAseq/CELL"
results_file <- file.path(project_dir, "CIBERSORTx_Results.txt")
metadata_file <- file.path(project_dir, "Sample_metadata.txt")
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")

fdr_cutoff <- 0.05
deconvolution_p_cutoff <- 0.05
sparse_detection_fraction <- 0.20

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(results_file)) stop("CIBERSORTx result file not found: ", results_file)
if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)

# Read and match input files
res <- read.delim(results_file, check.names = FALSE)
meta <- read.delim(metadata_file, check.names = FALSE)
required <- c("Sample", "Infant_ID", "Time_point", "Group")

if (!"Mixture" %in% names(res)) stop("The result file has no Mixture column.")
if (length(setdiff(required, names(meta)))) stop("Required metadata columns are missing.")
if (anyDuplicated(res$Mixture) || anyDuplicated(meta$Sample)) stop("Duplicate sample IDs found.")

rownames(res) <- res$Mixture
rownames(meta) <- meta$Sample
if (!setequal(rownames(res), rownames(meta))) stop("CIBERSORTx results and metadata do not match.")

meta <- meta[rownames(res), , drop = FALSE]
meta$Infant_ID <- factor(meta$Infant_ID)
meta$Group <- factor(meta$Group, levels = c("nonBPD", "BPD"))
meta$Time_point <- factor(meta$Time_point, levels = c("Cord", "Day14", "Day28"))
if (anyNA(meta[, required])) stop("Required metadata contain missing or unexpected values.")

# CIBERSORTx sample and cell-type QC
p_col <- intersect(c("P-value", "P.value", "P value"), names(res))
if (length(p_col) != 1) stop("Could not identify exactly one deconvolution P-value column.")

metric_cols <- intersect(c("P-value", "P.value", "P value", "Correlation", "RMSE", "Absolute score"), names(res))
cell_cols <- setdiff(names(res), c("Mixture", metric_cols))
if (length(cell_cols) < 2) stop("Fewer than two cell-fraction columns found.")

fractions_all <- as.matrix(res[, cell_cols, drop = FALSE])
storage.mode(fractions_all) <- "double"
if (anyNA(fractions_all) || any(!is.finite(fractions_all))) stop("Fractions contain invalid values.")
if (any(fractions_all < -1e-12) || any(fractions_all > 1 + 1e-12)) stop("Relative fractions must lie between zero and one.")
fractions_all[abs(fractions_all) < 1e-15] <- 0

fraction_sum <- rowSums(fractions_all)
if (any(abs(fraction_sum - 1) > 1e-6)) stop("Fractions do not sum to one; check relative mode.")

deconv_p <- as.numeric(res[[p_col]])
qc_pass <- is.finite(deconv_p) & deconv_p < deconvolution_p_cutoff
sample_qc <- data.frame(Sample = rownames(res), Fraction_sum = fraction_sum, QC_pass = qc_pass, meta[, c("Infant_ID", "Time_point", "Group")], check.names = FALSE)
for (column in metric_cols) sample_qc[[column]] <- res[[column]]
write.csv(sample_qc, file.path(results_dir, "CIBERSORTx_sample_QC.csv"), row.names = FALSE)

if (!all(qc_pass)) warning(sum(!qc_pass), " sample(s) failed CIBERSORTx P < ", deconvolution_p_cutoff, " and will be excluded.")
if (sum(qc_pass) < 10) stop("Fewer than 10 samples passed deconvolution QC.")

fractions <- fractions_all[qc_pass, , drop = FALSE]
meta <- droplevels(meta[qc_pass, , drop = FALSE])
stopifnot(identical(rownames(fractions), rownames(meta)))

detection_count <- colSums(fractions > 0)
detection_audit <- data.frame(
  CellType = cell_cols,
  Samples_detected = unname(detection_count[cell_cols]),
  Samples_tested = nrow(fractions),
  Detection_fraction = unname(detection_count[cell_cols]) / nrow(fractions),
  Mean_fraction = colMeans(fractions),
  Median_fraction = apply(fractions, 2, median),
  Maximum_fraction = apply(fractions, 2, max),
  Sparse_detection = unname(detection_count[cell_cols]) / nrow(fractions) < sparse_detection_fraction
)
write.csv(detection_audit, file.path(results_dir, "CIBERSORTx_cell_type_detection_audit.csv"), row.names = FALSE)

# Keep all cell types for the composition plots.
# Use only commonly detected cell types for statistical testing.
tested_cell_types <- detection_audit$CellType[!detection_audit$Sparse_detection]
excluded_cell_types <- detection_audit$CellType[detection_audit$Sparse_detection]

if (length(tested_cell_types) < 2) {
  stop("Fewer than two cell types passed the detection threshold.")
}

writeLines(
  c(
    paste("Detection threshold:", 100 * sparse_detection_fraction, "% of QC-passing samples"),
    paste("Cell types retained for testing:", length(tested_cell_types)),
    paste("Cell types excluded from testing:", length(excluded_cell_types)),
    if (length(excluded_cell_types)) excluded_cell_types else "None"
  ),
  file.path(results_dir, "CIBERSORTx_sparse_cell_types_excluded.txt")
)

# Transform the cell fractions used in the statistical analysis.
# Keep the original fractions for the composition plots.
fractions_model <- fractions[, tested_cell_types, drop = FALSE]
model_matrix <- t(asin(sqrt(fractions_model)))
colnames(model_matrix) <- rownames(meta)
constant_cell <- apply(model_matrix, 1, function(x) max(x) == min(x))
if (any(constant_cell)) {
  warning("Removing constant cell types: ", paste(rownames(model_matrix)[constant_cell], collapse = ", "))
  model_matrix <- model_matrix[!constant_cell, , drop = FALSE]
}

# Result and plotting functions
add_fraction_summary <- function(tab, reference_samples, comparison_samples) {
  reference_mean <- colMeans(fractions[reference_samples, , drop = FALSE])
  comparison_mean <- colMeans(fractions[comparison_samples, , drop = FALSE])
  tab$Mean_fraction_reference <- unname(reference_mean[tab$CellType])
  tab$Mean_fraction_comparison <- unname(comparison_mean[tab$CellType])
  tab$Raw_fraction_difference <- tab$Mean_fraction_comparison - tab$Mean_fraction_reference
  tab$Samples_detected <- unname(detection_count[tab$CellType])
  tab$Sparse_detection <- detection_audit$Sparse_detection[match(tab$CellType, detection_audit$CellType)]
  tab
}

save_results <- function(fit, coefficient, prefix, reference_samples, comparison_samples) {
  fit <- eBayes(fit, robust = TRUE)
  tab <- topTable(fit, coef = coefficient, number = Inf, sort.by = "P")
  tab$CellType <- rownames(tab)
  tab <- add_fraction_summary(tab, reference_samples, comparison_samples)
  tab <- tab[, c("CellType", "Mean_fraction_reference", "Mean_fraction_comparison", "Raw_fraction_difference", "Samples_detected", "Sparse_detection", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  tab <- tab[order(tab$adj.P.Val, tab$P.Value), ]
  sig <- tab[!is.na(tab$adj.P.Val) & tab$adj.P.Val < fdr_cutoff, , drop = FALSE]
  write.csv(tab, file.path(results_dir, paste0(prefix, "_full.csv")), row.names = FALSE)
  write.csv(sig, file.path(results_dir, paste0(prefix, "_significant.csv")), row.names = FALSE)
  list(full = tab, sig = sig)
}

fit_simple <- function(scores, metadata, design, coefficient, prefix, reference_samples, comparison_samples) {
  stopifnot(identical(colnames(scores), rownames(metadata)))
  fit <- lmFit(scores, design)
  save_results(fit, coefficient, prefix, reference_samples, comparison_samples)
}

make_result_plot <- function(result, title_text, out_file, top_n = 12) {
  plot_data <- result[!is.na(result$adj.P.Val) & result$adj.P.Val < fdr_cutoff, ]
  if (!nrow(plot_data)) return(invisible(NULL))

  plot_data <- head(plot_data[order(abs(plot_data$Raw_fraction_difference), decreasing = TRUE), ], top_n)
  plot_data$CellType <- factor(plot_data$CellType, levels = rev(plot_data$CellType))
  plot_data$Direction <- ifelse(plot_data$Raw_fraction_difference > 0, "Higher", "Lower")

  p <- ggplot(plot_data, aes(Raw_fraction_difference, CellType, fill = Direction)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.4) +
    scale_fill_manual(values = c(Higher = "#D55E00", Lower = "#0072B2")) +
    labs(title = title_text, subtitle = "BH FDR < 0.05", x = "Difference in estimated fraction", y = NULL, fill = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))

  ggsave(out_file, p, width = 9, height = 6.5, dpi = 300, bg = "white")
}

# Overall BPD effect adjusted for time and repeated samples
design_overall <- model.matrix(~ Time_point + Group, data = meta)
dupcor <- duplicateCorrelation(model_matrix, design_overall, block = meta$Infant_ID)
fit_overall <- lmFit(model_matrix, design_overall, block = meta$Infant_ID, correlation = dupcor$consensus.correlation)
overall <- save_results(fit_overall, "GroupBPD", "CELL_Overall_BPD_vs_nonBPD", rownames(meta)[meta$Group == "nonBPD"], rownames(meta)[meta$Group == "BPD"])
make_result_plot(overall$full, "Overall BPD vs nonBPD", file.path(figures_dir, "CELL_Overall_BPD_vs_nonBPD_barplot.png"))
writeLines(paste("Consensus repeated-measure correlation:", dupcor$consensus.correlation), file.path(results_dir, "CELL_Overall_repeated_measure_correlation.txt"))

# BPD versus nonBPD at each time point
timepoint_results <- list()
for (time in levels(meta$Time_point)) {
  metadata <- droplevels(meta[meta$Time_point == time, , drop = FALSE])
  scores <- model_matrix[, rownames(metadata), drop = FALSE]
  prefix <- paste0("CELL_", time, "_BPD_vs_nonBPD")
  timepoint_results[[time]] <- fit_simple(scores, metadata, model.matrix(~ Group, data = metadata), "GroupBPD", prefix, rownames(metadata)[metadata$Group == "nonBPD"], rownames(metadata)[metadata$Group == "BPD"])
  make_result_plot(timepoint_results[[time]]$full, paste0(time, ": BPD vs nonBPD"), file.path(figures_dir, paste0(prefix, "_barplot.png")))
}

# Paired longitudinal comparisons
run_progression <- function(group_name, reference_time, comparison_time) {
  metadata <- meta[meta$Group == group_name & meta$Time_point %in% c(reference_time, comparison_time), , drop = FALSE]
  pair_table <- table(metadata$Infant_ID, metadata$Time_point)
  complete_ids <- rownames(pair_table)[pair_table[, reference_time] == 1 & pair_table[, comparison_time] == 1]
  metadata <- droplevels(metadata[metadata$Infant_ID %in% complete_ids, , drop = FALSE])
  metadata$Time_point <- factor(metadata$Time_point, levels = c(reference_time, comparison_time))
  if (length(complete_ids) < 3) stop("Fewer than three complete pairs for ", group_name, ": ", comparison_time, " vs ", reference_time)

  scores <- model_matrix[, rownames(metadata), drop = FALSE]
  prefix <- paste0("CELL_", group_name, "_paired_", comparison_time, "_vs_", reference_time)
  result <- fit_simple(scores, metadata, model.matrix(~ Infant_ID + Time_point, data = metadata), paste0("Time_point", comparison_time), prefix, rownames(metadata)[metadata$Time_point == reference_time], rownames(metadata)[metadata$Time_point == comparison_time])
  result$Pair_count <- length(complete_ids)
  make_result_plot(result$full, paste0(group_name, ": ", comparison_time, " vs ", reference_time), file.path(figures_dir, paste0(prefix, "_barplot.png")))
  result
}

progression_results <- list(
  BPD_Day14_vs_Cord = run_progression("BPD", "Cord", "Day14"),
  BPD_Day28_vs_Day14 = run_progression("BPD", "Day14", "Day28"),
  BPD_Day28_vs_Cord = run_progression("BPD", "Cord", "Day28"),
  nonBPD_Day14_vs_Cord = run_progression("nonBPD", "Cord", "Day14"),
  nonBPD_Day28_vs_Day14 = run_progression("nonBPD", "Day14", "Day28"),
  nonBPD_Day28_vs_Cord = run_progression("nonBPD", "Cord", "Day28")
)

# Cell composition plot
long_data <- reshape(data.frame(Sample = rownames(fractions), fractions, check.names = FALSE), varying = cell_cols, v.names = "Fraction", timevar = "CellType", times = cell_cols, direction = "long")
long_data$Group <- meta[long_data$Sample, "Group"]
long_data$Time_point <- meta[long_data$Sample, "Time_point"]
long_data$Panel <- interaction(long_data$Group, long_data$Time_point, sep = ": ", drop = TRUE)
long_data$Panel <- factor(long_data$Panel, levels = c("nonBPD: Cord", "nonBPD: Day14", "nonBPD: Day28", "BPD: Cord", "BPD: Day14", "BPD: Day28"))

composition_plot <- ggplot(long_data, aes(Sample, Fraction, fill = CellType)) +
  geom_col(width = 0.95) +
  facet_wrap(~ Panel, nrow = 2, scales = "free_x") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "CIBERSORTx-estimated immune-cell composition", x = NULL, y = "Relative fraction", fill = "Cell type") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing.x = unit(0.8, "lines"), plot.title = element_text(face = "bold"))
ggsave(file.path(figures_dir, "CELL_composition_all_samples.png"), composition_plot, width = 16, height = 8, dpi = 300, bg = "white")

# Summary and run information
all_results <- c(list(Overall_BPD_vs_nonBPD = overall), setNames(timepoint_results, paste0(names(timepoint_results), "_BPD_vs_nonBPD")), progression_results)
summary_df <- data.frame(
  Analysis = names(all_results),
  Samples_or_pairs = unname(c(nrow(meta), vapply(levels(meta$Time_point), function(time) sum(meta$Time_point == time), integer(1)), vapply(progression_results, function(x) x$Pair_count, integer(1)))),
  Tested_cell_types = vapply(all_results, function(x) nrow(x$full), integer(1)),
  Significant_cell_types = vapply(all_results, function(x) nrow(x$sig), integer(1))
)
write.csv(summary_df, file.path(results_dir, "CELL_analysis_summary.csv"), row.names = FALSE)

writeLines(c(
  paste("Run time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "CIBERSORTx mode: relative fractions",
  paste("CIBERSORTx sample QC: deconvolution P-value <", deconvolution_p_cutoff),
  "Fraction transformation: arcsine square root",
  "Differential testing: limma with robust empirical Bayes",
  paste("FDR cutoff:", fdr_cutoff),
  paste("Sparse cell types excluded from differential testing when detected in less than", 100 * sparse_detection_fraction, "% of QC-passing samples"),
  paste("limma version:", as.character(packageVersion("limma")))
), file.path(results_dir, "CELL_run_parameters.txt"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "CELL_sessionInfo.txt"))

cat("\nQC-passing samples:", nrow(meta), "of", nrow(res), "\n")
cat("Cell types tested:", nrow(model_matrix), "\n")
cat("Sparse cell types excluded from testing:", length(excluded_cell_types), "\n")
if (length(excluded_cell_types)) cat("Excluded:", paste(excluded_cell_types, collapse = ", "), "\n")
cat("Repeated-measure correlation:", dupcor$consensus.correlation, "\n")
cat("\nAnalysis summary:\n")
print(summary_df, row.names = FALSE)
cat("\nResults saved in:", results_dir, "\n")
cat("Figures saved in:", figures_dir, "\n")
cat("Cell-type analysis completed.\n")


