rm(list = ls())
gc()

# Packages and file paths
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(apeglm)
  library(openxlsx)
})

project_dir <- "/home/sagar007/RNAseq/DEG"
counts_file <- "/home/sagar007/RNAseq/counts/gene_counts_reverse.tsv"
metadata_file <- file.path(project_dir, "Sample_metadata.txt")
results_dir <- file.path(project_dir, "results")
plot_dir <- file.path(project_dir, "plots")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

fdr_cutoff <- 0.05
lfc_cutoff <- 1
min_count <- 10
min_samples <- 8

group_colours <- c(nonBPD = "#4DBBD5", BPD = "#E64B35")
time_colours <- c(Cord = "#3C5488", Day14 = "#00A087", Day28 = "#F39B7F")
time_shapes <- c(Cord = 16, Day14 = 17, Day28 = 15)

if (!file.exists(counts_file)) stop("Count file not found: ", counts_file)
if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)

# Read counts and metadata
counts_data <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
gene_column <- intersect(c("Geneid", "GeneID", "gene_id"), names(counts_data))[1]
if (is.na(gene_column)) gene_column <- names(counts_data)[1]

rownames(counts_data) <- counts_data[[gene_column]]
annotation_columns <- intersect(c("Geneid", "GeneID", "gene_id", "Chr", "Start", "End", "Strand", "Length"), names(counts_data))
counts <- as.matrix(counts_data[, !names(counts_data) %in% annotation_columns, drop = FALSE])
colnames(counts) <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(colnames(counts)))
storage.mode(counts) <- "integer"

metadata <- read.delim(metadata_file, check.names = FALSE)
required_columns <- c("Sample", "Infant_ID", "Group", "Time_point")

if (length(setdiff(required_columns, names(metadata))) > 0) stop("Required metadata columns are missing.")
if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)) || anyDuplicated(metadata$Sample)) stop("Duplicate gene or sample IDs found.")

rownames(metadata) <- metadata$Sample
metadata$Infant_ID <- factor(metadata$Infant_ID)
metadata$Group <- factor(metadata$Group, levels = c("nonBPD", "BPD"))
metadata$Time_point <- factor(metadata$Time_point, levels = c("Cord", "Day14", "Day28"))

if (anyNA(metadata[, required_columns])) stop("Metadata contains missing or unexpected values.")
if (!setequal(colnames(counts), rownames(metadata))) stop("Count and metadata sample IDs do not match.")

metadata <- metadata[colnames(counts), , drop = FALSE]

cat("\nSample distribution:\n")
print(table(metadata$Group, metadata$Time_point))
cat("Infants:", nlevels(metadata$Infant_ID), "\n")
cat("Samples:", nrow(metadata), "\n")

# Filter low-count genes
keep_gene <- rowSums(counts >= min_count) >= min_samples
counts_filtered <- counts[keep_gene, , drop = FALSE]

cat("Genes before filtering:", nrow(counts), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")

write.table(data.frame(Geneid = rownames(counts_filtered), counts_filtered, check.names = FALSE), file.path(results_dir, "filtered_gene_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Result and plotting functions
add_symbols <- function(results) {
  results$Ensembl <- rownames(results)
  results$Ensembl_clean <- sub("\\..*$", "", results$Ensembl)
  results$Symbol <- mapIds(org.Hs.eg.db, keys = results$Ensembl_clean, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  results
}

significant_genes <- function(results) {
  results[!is.na(results$padj) & results$padj < fdr_cutoff & abs(results$log2FoldChange) >= lfc_cutoff, , drop = FALSE]
}

format_deseq_results <- function(raw_results, shrunk_results) {
  raw <- as.data.frame(raw_results)
  shrunk <- as.data.frame(shrunk_results)
  raw$log2FoldChange_raw <- raw$log2FoldChange
  raw$log2FoldChange <- shrunk$log2FoldChange
  raw <- add_symbols(raw)
  raw <- raw[, c("Ensembl", "Ensembl_clean", "Symbol", "baseMean", "log2FoldChange", "log2FoldChange_raw", "lfcSE", "stat", "pvalue", "padj")]
  raw[order(raw$padj, na.last = TRUE), , drop = FALSE]
}

save_result_tables <- function(results, prefix) {
  significant <- significant_genes(results)
  ranked <- results[!is.na(results$stat), c("Ensembl", "Ensembl_clean", "Symbol", "stat"), drop = FALSE]
  ranked <- ranked[order(ranked$stat, decreasing = TRUE), ]

  write.csv(results, file.path(results_dir, paste0(prefix, "_full.csv")), row.names = FALSE)
  write.csv(significant, file.path(results_dir, paste0(prefix, "_significant.csv")), row.names = FALSE)
  write.csv(ranked, file.path(results_dir, paste0(prefix, "_FGSEA_ranked_statistics.csv")), row.names = FALSE)
  significant
}

make_volcano <- function(results, title, prefix, lfc_label = "Shrunken log2 fold change") {
  plot_data <- results
  plot_data$Status <- "Not significant"
  plot_data$Status[!is.na(plot_data$padj) & plot_data$padj < fdr_cutoff & plot_data$log2FoldChange >= lfc_cutoff] <- "Up"
  plot_data$Status[!is.na(plot_data$padj) & plot_data$padj < fdr_cutoff & plot_data$log2FoldChange <= -lfc_cutoff] <- "Down"
  plot_data$minus_log10_FDR <- -log10(pmax(plot_data$padj, .Machine$double.xmin))

  labels <- plot_data[plot_data$Status != "Not significant" & !is.na(plot_data$Symbol), ]
  labels <- head(labels[order(labels$padj), ], 10)

  p <- ggplot(plot_data, aes(log2FoldChange, minus_log10_FDR, colour = Status)) +
    geom_point(size = 1.6, alpha = 0.70) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2, colour = "grey45") +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = 2, colour = "grey45") +
    geom_text(data = labels, aes(label = Symbol), colour = "black", size = 3, check_overlap = TRUE, vjust = -0.7) +
    scale_colour_manual(values = c(Down = "#3C5488", `Not significant` = "grey75", Up = "#E64B35")) +
    labs(title = title, subtitle = paste0("BH-FDR < 0.05 and |", lfc_label, "| >= 1"), x = lfc_label, y = expression(-log[10](FDR)), colour = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))

  ggsave(file.path(plot_dir, paste0("Volcano_", prefix, ".png")), p, width = 8, height = 6.5, dpi = 300, bg = "white")
}

make_ma_plot <- function(results, title, prefix, lfc_label = "Shrunken log2 fold change") {
  plot_data <- results
  plot_data$Significant <- !is.na(plot_data$padj) & plot_data$padj < fdr_cutoff & abs(plot_data$log2FoldChange) >= lfc_cutoff

  p <- ggplot(plot_data, aes(baseMean, log2FoldChange, colour = Significant)) +
    geom_point(size = 1.3, alpha = 0.65) +
    geom_hline(yintercept = 0, colour = "grey35") +
    scale_x_log10() +
    scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#E64B35")) +
    labs(title = title, x = "Mean normalized expression", y = lfc_label, colour = "Significant") +
    theme_classic(base_size = 12) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))

  ggsave(file.path(plot_dir, paste0("MAplot_", prefix, ".png")), p, width = 8, height = 6, dpi = 300, bg = "white")
}

make_pca <- function(vsd, metadata, title, prefix, colour_by, shape_by = NULL) {
  pca <- prcomp(t(assay(vsd)))
  percent_var <- round(100 * summary(pca)$importance[2, 1:2], 1)

  plot_data <- data.frame(Sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  plot_data <- cbind(plot_data, metadata[plot_data$Sample, c("Infant_ID", "Group", "Time_point"), drop = FALSE])
  plot_data$Plot_colour <- if (colour_by == "Group") unname(group_colours[plot_data$Group]) else unname(time_colours[plot_data$Time_point])
  plot_data$Plot_shape <- if (is.null(shape_by)) 16 else unname(time_shapes[plot_data$Time_point])

  p <- ggplot(plot_data, aes(PC1, PC2, colour = .data[[colour_by]]))
  if (is.null(shape_by)) p <- p + geom_point(size = 3.8, alpha = 0.90)
  if (!is.null(shape_by)) p <- p + geom_point(aes(shape = .data[[shape_by]]), size = 3.8, alpha = 0.90)
  if (colour_by == "Group") p <- p + scale_colour_manual(values = group_colours)
  if (colour_by == "Time_point") p <- p + scale_colour_manual(values = time_colours)
  if (!is.null(shape_by)) p <- p + scale_shape_manual(values = time_shapes)
  p <- p + labs(title = title, x = paste0("PC1 (", percent_var[1], "%)"), y = paste0("PC2 (", percent_var[2], "%)"), colour = colour_by, shape = shape_by)
  p <- p + theme_classic(base_size = 12) + theme(plot.title = element_text(face = "bold"), legend.position = "right")

  ggsave(file.path(plot_dir, paste0("PCA_", prefix, ".png")), p, width = 8, height = 6, dpi = 300, bg = "white")
  plot_data
}

make_heatmap <- function(vsd, significant, metadata, title, prefix) {
  genes <- head(significant$Ensembl, 50)
  genes <- intersect(genes, rownames(vsd))
  if (length(genes) < 2) return(invisible(NULL))

  matrix <- assay(vsd)[genes, , drop = FALSE]
  matrix <- t(scale(t(matrix)))
  matrix <- matrix[complete.cases(matrix), , drop = FALSE]
  annotation <- metadata[, c("Group", "Time_point"), drop = FALSE]
  annotation_colours <- list(Group = group_colours, Time_point = time_colours)

  png(file.path(plot_dir, paste0("Heatmap_Top50_", prefix, ".png")), 2200, 2200, res = 300)
  pheatmap(matrix, annotation_col = annotation, annotation_colors = annotation_colours, show_rownames = FALSE, show_colnames = FALSE, border_color = NA, color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100), main = title)
  dev.off()
}

# DESeq2 analysis for time-specific and paired comparisons
run_deseq <- function(counts, metadata, design, coefficient, title, prefix, pca_colour, pca_shape = NULL) {
  cat("\nRunning:", title, "\n")

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = design)
  dds <- DESeq(dds)
  if (!coefficient %in% resultsNames(dds)) stop("Coefficient not found: ", coefficient)

  raw <- results(dds, name = coefficient, alpha = fdr_cutoff)
  shrunk <- lfcShrink(dds, coef = coefficient, type = "apeglm")
  results <- format_deseq_results(raw, shrunk)
  significant <- save_result_tables(results, prefix)
  vsd <- vst(dds, blind = FALSE)

  make_volcano(results, title, prefix)
  make_ma_plot(results, title, prefix)
  pca_data <- make_pca(vsd, metadata, title, prefix, pca_colour, pca_shape)
  make_heatmap(vsd, significant, metadata, title, prefix)

  cat("Significant DEGs:", nrow(significant), "\n")
  list(dds = dds, vsd = vsd, results = results, significant = significant, pca = pca_data)
}

# Overall repeated-measures analysis
run_overall <- function(counts, metadata) {
  title <- "Overall BPD vs nonBPD"
  prefix <- "DEG_Overall_BPD_vs_nonBPD"
  design <- model.matrix(~ Time_point + Group, metadata)

  y <- calcNormFactors(DGEList(counts))
  voom_data <- voom(y, design, plot = FALSE)
  correlation <- duplicateCorrelation(voom_data, design, block = metadata$Infant_ID)$consensus.correlation
  voom_data <- voom(y, design, plot = FALSE, block = metadata$Infant_ID, correlation = correlation)
  correlation <- duplicateCorrelation(voom_data, design, block = metadata$Infant_ID)$consensus.correlation
  fit <- eBayes(lmFit(voom_data, design, block = metadata$Infant_ID, correlation = correlation))

  table <- topTable(fit, coef = "GroupBPD", number = Inf, sort.by = "P")
  mean_cpm <- rowMeans(cpm(y, normalized.lib.sizes = TRUE))
  table$baseMean <- mean_cpm[rownames(table)]
  table$log2FoldChange <- table$logFC
  table$log2FoldChange_raw <- table$logFC
  table$lfcSE <- ifelse(table$t == 0, NA_real_, abs(table$logFC / table$t))
  table$stat <- table$t
  table$pvalue <- table$P.Value
  table$padj <- table$adj.P.Val
  table <- add_symbols(table)
  table <- table[, c("Ensembl", "Ensembl_clean", "Symbol", "baseMean", "log2FoldChange", "log2FoldChange_raw", "lfcSE", "stat", "pvalue", "padj")]
  table <- table[order(table$padj), ]

  significant <- save_result_tables(table, prefix)

  dds <- estimateSizeFactors(DESeqDataSetFromMatrix(counts, metadata, ~ Time_point + Group))
  vsd <- vst(dds, blind = FALSE)

  make_volcano(table, title, prefix, "Log2 fold change")
  make_ma_plot(table, title, prefix, "Log2 fold change")
  pca_data <- make_pca(vsd, metadata, title, prefix, "Group", "Time_point")
  make_heatmap(vsd, significant, metadata, title, prefix)

  writeLines(paste("Repeated-measure correlation =", correlation), file.path(results_dir, "Overall_repeated_measure_correlation.txt"))
  cat("\nRunning:", title, "\n")
  cat("Repeated-measure correlation:", correlation, "\n")
  cat("Significant DEGs:", nrow(significant), "\n")

  list(results = table, significant = significant, pca = pca_data)
}

# Run the ten comparisons
overall <- run_overall(counts_filtered, metadata)

run_timepoint <- function(time_point) {
  keep <- metadata$Time_point == time_point
  meta <- droplevels(metadata[keep, , drop = FALSE])
  run_deseq(counts_filtered[, keep, drop = FALSE], meta, ~ Group, "Group_BPD_vs_nonBPD", paste(time_point, "BPD vs nonBPD"), paste0("DEG_", time_point, "_BPD_vs_nonBPD"), "Group")
}

cord <- run_timepoint("Cord")
day14 <- run_timepoint("Day14")
day28 <- run_timepoint("Day28")

run_paired <- function(group, reference, comparison) {
  keep <- metadata$Group == group & metadata$Time_point %in% c(reference, comparison)
  meta <- droplevels(metadata[keep, , drop = FALSE])
  meta$Time_point <- factor(meta$Time_point, levels = c(reference, comparison))
  if (any(table(meta$Infant_ID, meta$Time_point) != 1)) stop("Incomplete sample pairs found.")

  coefficient <- paste0("Time_point_", comparison, "_vs_", reference)
  prefix <- paste0("DEG_", group, "_paired_", comparison, "_vs_", reference)
  title <- paste(group, comparison, "vs", reference)
  run_deseq(counts_filtered[, keep, drop = FALSE], meta, ~ Infant_ID + Time_point, coefficient, title, prefix, "Time_point")
}

bpd_day14_cord <- run_paired("BPD", "Cord", "Day14")
bpd_day28_day14 <- run_paired("BPD", "Day14", "Day28")
bpd_day28_cord <- run_paired("BPD", "Cord", "Day28")
nonbpd_day14_cord <- run_paired("nonBPD", "Cord", "Day14")
nonbpd_day28_day14 <- run_paired("nonBPD", "Day14", "Day28")
nonbpd_day28_cord <- run_paired("nonBPD", "Cord", "Day28")

# Summary and PCA sample key
analysis_names <- c("Overall_BPD_vs_nonBPD", "Cord_BPD_vs_nonBPD", "Day14_BPD_vs_nonBPD", "Day28_BPD_vs_nonBPD", "BPD_Day14_vs_Cord_paired", "BPD_Day28_vs_Day14_paired", "BPD_Day28_vs_Cord_paired", "nonBPD_Day14_vs_Cord_paired", "nonBPD_Day28_vs_Day14_paired", "nonBPD_Day28_vs_Cord_paired")
analysis_results <- list(overall, cord, day14, day28, bpd_day14_cord, bpd_day28_day14, bpd_day28_cord, nonbpd_day14_cord, nonbpd_day28_day14, nonbpd_day28_cord)

summary_table <- data.frame(Analysis = analysis_names, Significant_DEGs = vapply(analysis_results, function(x) nrow(x$significant), integer(1)))
write.csv(summary_table, file.path(results_dir, "DEG_analysis_summary.csv"), row.names = FALSE)

shared_symbols <- Reduce(intersect, list(na.omit(cord$significant$Symbol), na.omit(day14$significant$Symbol), na.omit(day28$significant$Symbol)))
write.csv(data.frame(Shared_Gene_Symbol = shared_symbols), file.path(results_dir, "Shared_significant_DEGs_Cord_Day14_Day28.csv"), row.names = FALSE)

pca_sheets <- setNames(lapply(analysis_results, function(x) x$pca), c("Overall", "Cord", "Day14", "Day28", "BPD_D14_Cord", "BPD_D28_D14", "BPD_D28_Cord", "nonBPD_D14_Cord", "nonBPD_D28_D14", "nonBPD_D28_Cord"))
write.xlsx(pca_sheets, file.path(results_dir, "PCA_sample_coordinates.xlsx"), rowNames = FALSE, overwrite = TRUE)

writeLines(capture.output(sessionInfo()), file.path(results_dir, "R_sessionInfo.txt"))

cat("\nAnalysis summary:\n")
print(summary_table, row.names = FALSE)
cat("Shared significant genes across Cord, Day14 and Day28:", length(shared_symbols), "\n")
cat("Results saved in:", results_dir, "\n")
cat("Plots saved in:", plot_dir, "\n")
cat("DEG analysis completed.\n")

