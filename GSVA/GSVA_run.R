rm(list = ls())
gc()

# Packages and file paths
suppressPackageStartupMessages({
  library(GSVA)
  library(limma)
  library(msigdbr)
  library(BiocParallel)
  library(ggplot2)
  library(pheatmap)
})

project_dir <- "/home/sagar007/RNAseq/GSVA"
input_dir <- file.path(project_dir, "inputs")
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

fdr_cutoff <- 0.05
minimum_set_size <- 10
maximum_set_size <- 500

positive_colour <- "#E64B35"
negative_colour <- "#3C5488"
group_colours <- c(nonBPD = "#4DBBD5", BPD = "#E64B35")
time_colours <- c(Cord = "#3C5488", Day14 = "#00A087", Day28 = "#F39B7F")

# Read normalized expression and metadata
expression_file <- file.path(input_dir, "normalized_expression_SYMBOL.csv")
metadata_file <- file.path(input_dir, "Sample_metadata.txt")

if (!file.exists(expression_file) || !file.exists(metadata_file)) stop("GSVA input files are missing.")

expression <- as.matrix(read.csv(expression_file, row.names = 1, check.names = FALSE))
storage.mode(expression) <- "double"
metadata <- read.delim(metadata_file, check.names = FALSE)

required_columns <- c("Sample", "Infant_ID", "Group", "Time_point")
if (length(setdiff(required_columns, names(metadata))) > 0) stop("Required metadata columns are missing.")
if (anyDuplicated(rownames(expression)) || anyDuplicated(metadata$Sample)) stop("Duplicate gene or sample IDs found.")
if (any(!is.finite(expression))) stop("Expression matrix contains invalid values.")

rownames(metadata) <- metadata$Sample
metadata$Infant_ID <- factor(metadata$Infant_ID)
metadata$Group <- factor(metadata$Group, levels = c("nonBPD", "BPD"))
metadata$Time_point <- factor(metadata$Time_point, levels = c("Cord", "Day14", "Day28"))

if (anyNA(metadata[, required_columns])) stop("Metadata contains missing or unexpected values.")
if (!setequal(colnames(expression), rownames(metadata))) stop("Expression and metadata sample IDs do not match.")

metadata <- metadata[colnames(expression), , drop = FALSE]

# Hallmark pathway definitions
hallmark_data <- msigdbr(db_species = "HS", species = "Homo sapiens", collection = "H")
hallmark_original <- split(hallmark_data$gene_symbol, hallmark_data$gs_name)
hallmark_original <- lapply(hallmark_original, unique)

overlap_size <- vapply(hallmark_original, function(genes) sum(genes %in% rownames(expression)), integer(1))
pathway_audit <- data.frame(Pathway = names(hallmark_original), Original_size = lengths(hallmark_original), Overlap_size = overlap_size, Overlap_fraction = overlap_size / lengths(hallmark_original))
pathway_audit$Tested <- pathway_audit$Overlap_size >= minimum_set_size & pathway_audit$Overlap_size <= maximum_set_size
write.csv(pathway_audit, file.path(results_dir, "GSVA_Hallmark_overlap_audit.csv"), row.names = FALSE)

hallmark_sets <- lapply(hallmark_original, intersect, y = rownames(expression))
hallmark_sets <- hallmark_sets[lengths(hallmark_sets) >= minimum_set_size & lengths(hallmark_sets) <= maximum_set_size]
if (length(hallmark_sets) == 0) stop("No Hallmark pathways remain after filtering.")

# Calculate GSVA scores
gsva_parameters <- gsvaParam(exprData = expression, geneSets = hallmark_sets, kcdf = "Gaussian", minSize = minimum_set_size, maxSize = maximum_set_size, tau = 1, maxDiff = TRUE, absRanking = FALSE)
gsva_scores <- gsva(gsva_parameters, verbose = TRUE, BPPARAM = SerialParam())

if (any(!is.finite(gsva_scores))) stop("GSVA returned invalid scores.")
write.csv(gsva_scores, file.path(results_dir, "GSVA_Hallmark_scores_all_samples.csv"))

gene_set_size <- lengths(hallmark_sets)[rownames(gsva_scores)]

# Result and plotting functions
clean_pathway <- function(pathway) {
  pathway <- sub("^HALLMARK_", "", pathway)
  tools::toTitleCase(tolower(gsub("_", " ", pathway)))
}

save_results <- function(table, prefix) {
  table$Pathway <- rownames(table)
  table$Gene_set_size <- unname(gene_set_size[table$Pathway])
  table <- table[, c("Pathway", "Gene_set_size", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  table <- table[order(table$adj.P.Val, table$P.Value), ]
  significant <- table[!is.na(table$adj.P.Val) & table$adj.P.Val < fdr_cutoff, , drop = FALSE]

  write.csv(table, file.path(results_dir, paste0(prefix, "_full.csv")), row.names = FALSE)
  write.csv(significant, file.path(results_dir, paste0(prefix, "_significant.csv")), row.names = FALSE)
  list(full = table, significant = significant)
}

finish_model <- function(fit, coefficient, prefix) {
  if (!coefficient %in% colnames(fit$coefficients)) stop("Coefficient not found: ", coefficient)
  fit <- eBayes(fit, robust = TRUE, trend = sqrt(gene_set_size))
  table <- topTable(fit, coef = coefficient, number = Inf, sort.by = "P")
  save_results(table, prefix)
}

make_barplot <- function(results, title, prefix, top_n = 15) {
  plot_data <- results[!is.na(results$adj.P.Val) & results$adj.P.Val < fdr_cutoff, , drop = FALSE]
  plot_data <- head(plot_data[order(abs(plot_data$logFC), decreasing = TRUE), ], top_n)
  if (nrow(plot_data) == 0) return(invisible(NULL))

  plot_data$Pathway_label <- factor(clean_pathway(plot_data$Pathway), levels = rev(clean_pathway(plot_data$Pathway)))
  plot_data$Direction <- ifelse(plot_data$logFC > 0, "Higher", "Lower")

  p <- ggplot(plot_data, aes(logFC, Pathway_label, fill = Direction)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, colour = "grey35") +
    scale_fill_manual(values = c(Higher = positive_colour, Lower = negative_colour)) +
    labs(title = title, subtitle = "BH-FDR < 0.05", x = "Mean GSVA-score difference", y = NULL, fill = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))

  ggsave(file.path(figures_dir, paste0(prefix, "_barplot.png")), p, width = 9, height = 7, dpi = 300, bg = "white")
}

make_heatmap <- function(results, scores, metadata, title, prefix, top_n = 20) {
  significant <- results[!is.na(results$adj.P.Val) & results$adj.P.Val < fdr_cutoff, , drop = FALSE]
  pathways <- head(significant$Pathway[order(abs(significant$logFC), decreasing = TRUE)], top_n)
  pathways <- intersect(pathways, rownames(scores))
  if (length(pathways) < 2) return(invisible(NULL))

  matrix <- scores[pathways, , drop = FALSE]
  rownames(matrix) <- clean_pathway(rownames(matrix))
  annotation <- metadata[, c("Group", "Time_point"), drop = FALSE]
  annotation_colours <- list(Group = group_colours, Time_point = time_colours)

  png(file.path(figures_dir, paste0(prefix, "_heatmap.png")), 2400, 1900, res = 300)
  pheatmap(matrix, scale = "row", annotation_col = annotation, annotation_colors = annotation_colours, cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = FALSE, border_color = NA, color = colorRampPalette(c(negative_colour, "white", positive_colour))(100), fontsize_row = 9, main = title)
  dev.off()
}

save_figures <- function(result, scores, metadata, title, prefix) {
  make_barplot(result$full, title, prefix)
  make_heatmap(result$full, scores, metadata, title, prefix)
}

# Overall BPD effect adjusted for time and repeated samples
overall_design <- model.matrix(~ Time_point + Group, metadata)
repeated_correlation <- duplicateCorrelation(gsva_scores, overall_design, block = metadata$Infant_ID)$consensus.correlation
overall_fit <- lmFit(gsva_scores, overall_design, block = metadata$Infant_ID, correlation = repeated_correlation)
overall <- finish_model(overall_fit, "GroupBPD", "GSVA_Overall_BPD_vs_nonBPD")

writeLines(paste("Consensus repeated-measure correlation:", repeated_correlation), file.path(results_dir, "GSVA_Overall_repeated_measure_correlation.txt"))
save_figures(overall, gsva_scores, metadata, "Overall BPD vs nonBPD", "GSVA_Overall_BPD_vs_nonBPD")

# BPD versus nonBPD at each time point
run_timepoint <- function(time_point) {
  keep <- metadata$Time_point == time_point
  meta <- droplevels(metadata[keep, , drop = FALSE])
  scores <- gsva_scores[, rownames(meta), drop = FALSE]
  design <- model.matrix(~ Group, meta)
  prefix <- paste0("GSVA_", time_point, "_BPD_vs_nonBPD")
  result <- finish_model(lmFit(scores, design), "GroupBPD", prefix)
  save_figures(result, scores, meta, paste(time_point, "BPD vs nonBPD"), prefix)
  result
}

cord <- run_timepoint("Cord")
day14 <- run_timepoint("Day14")
day28 <- run_timepoint("Day28")

# Paired longitudinal comparisons
run_paired <- function(group, reference, comparison) {
  keep <- metadata$Group == group & metadata$Time_point %in% c(reference, comparison)
  meta <- droplevels(metadata[keep, , drop = FALSE])
  meta$Time_point <- factor(meta$Time_point, levels = c(reference, comparison))
  if (any(table(meta$Infant_ID, meta$Time_point) != 1)) stop("Incomplete sample pairs found.")

  scores <- gsva_scores[, rownames(meta), drop = FALSE]
  design <- model.matrix(~ Infant_ID + Time_point, meta)
  coefficient <- paste0("Time_point", comparison)
  prefix <- paste0("GSVA_", group, "_paired_", comparison, "_vs_", reference)
  result <- finish_model(lmFit(scores, design), coefficient, prefix)
  save_figures(result, scores, meta, paste(group, comparison, "vs", reference), prefix)
  result
}

bpd_day14_cord <- run_paired("BPD", "Cord", "Day14")
bpd_day28_day14 <- run_paired("BPD", "Day14", "Day28")
bpd_day28_cord <- run_paired("BPD", "Cord", "Day28")
nonbpd_day14_cord <- run_paired("nonBPD", "Cord", "Day14")
nonbpd_day28_day14 <- run_paired("nonBPD", "Day14", "Day28")
nonbpd_day28_cord <- run_paired("nonBPD", "Cord", "Day28")

# Analysis summary
analysis_names <- c("Overall_BPD_vs_nonBPD", "Cord_BPD_vs_nonBPD", "Day14_BPD_vs_nonBPD", "Day28_BPD_vs_nonBPD", "BPD_Day14_vs_Cord", "BPD_Day28_vs_Day14", "BPD_Day28_vs_Cord", "nonBPD_Day14_vs_Cord", "nonBPD_Day28_vs_Day14", "nonBPD_Day28_vs_Cord")
analysis_results <- list(overall, cord, day14, day28, bpd_day14_cord, bpd_day28_day14, bpd_day28_cord, nonbpd_day14_cord, nonbpd_day28_day14, nonbpd_day28_cord)

summary_table <- data.frame(Analysis = analysis_names, Tested_pathways = vapply(analysis_results, function(x) nrow(x$full), integer(1)), Significant_pathways = vapply(analysis_results, function(x) nrow(x$significant), integer(1)))
write.csv(summary_table, file.path(results_dir, "GSVA_analysis_summary.csv"), row.names = FALSE)

writeLines(c("Input: DESeq2 variance-stabilized expression values", "GSVA kernel: Gaussian", paste("Gene-set size:", minimum_set_size, "to", maximum_set_size), paste("FDR cutoff:", fdr_cutoff), "Differential testing: limma with robust eBayes trend by sqrt(gene-set size)", paste("GSVA version:", packageVersion("GSVA")), paste("limma version:", packageVersion("limma")), paste("msigdbr version:", packageVersion("msigdbr"))), file.path(results_dir, "GSVA_run_parameters.txt"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "GSVA_sessionInfo.txt"))

cat("\nExpression matrix:", nrow(expression), "genes x", ncol(expression), "samples\n")
cat("Hallmark pathways tested:", nrow(gsva_scores), "\n")
cat("Repeated-measure correlation:", repeated_correlation, "\n")
cat("\nAnalysis summary:\n")
print(summary_table, row.names = FALSE)
cat("\nResults saved in:", results_dir, "\n")
cat("Figures saved in:", figures_dir, "\n")
cat("GSVA analysis completed.\n")

