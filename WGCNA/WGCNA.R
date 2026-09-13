rm(list = ls())
gc()

# Packages and file paths
suppressPackageStartupMessages({
  library(WGCNA)
  library(lmerTest)
  library(emmeans)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()
set.seed(54321)

project_dir <- "/home/sagar007/RNAseq/WGCNA"
expr_file <- "/home/sagar007/RNAseq/GSVA/inputs/normalized_expression_SYMBOL.csv"
metadata_file <- file.path(project_dir, "Sample_metadata.txt")

plot_dir <- file.path(project_dir, "plots")
results_dir <- file.path(project_dir, "results")
module_dir <- file.path(project_dir, "module_gene_lists")
intermediate_dir <- file.path(project_dir, "intermediate")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(module_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

# Analysis settings
network_type <- "signed"
cor_type <- "bicor"
max_p_outliers <- 0.10
powers <- c(1:10, seq(12, 30, 2))
scale_free_target <- 0.85
fallback_power <- 12
min_module_size <- 30
merge_cut_height <- 0.15
max_block_size <- 5000
hub_top_n <- 20
fdr_cutoff <- 0.05

if (!file.exists(expr_file)) stop("Expression file not found: ", expr_file)
if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)

# Read expression and metadata
expr <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))
storage.mode(expr) <- "double"
datExpr <- as.data.frame(t(expr), check.names = FALSE)

metadata <- read.delim(metadata_file, check.names = FALSE)
required_cols <- c("Sample", "Infant_ID", "Group", "Time_point")

if (length(setdiff(required_cols, names(metadata))) > 0) stop("Required metadata columns are missing.")
if (anyDuplicated(metadata$Sample)) stop("Metadata contains duplicate sample IDs.")
if (any(!is.finite(as.matrix(datExpr)))) stop("Expression matrix contains missing or non-finite values.")

rownames(metadata) <- metadata$Sample
metadata$Infant_ID <- factor(metadata$Infant_ID)
metadata$Group <- factor(metadata$Group, levels = c("nonBPD", "BPD"))
metadata$Time_point <- factor(metadata$Time_point, levels = c("Cord", "Day14", "Day28"))

if (anyNA(metadata[, required_cols])) stop("Metadata contains missing or unexpected values.")
if (!setequal(rownames(datExpr), rownames(metadata))) stop("Expression and metadata sample IDs do not match.")

metadata <- metadata[rownames(datExpr), , drop = FALSE]

cat("Expression matrix:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")
cat("\nSample distribution:\n")
print(table(metadata$Group, metadata$Time_point))
cat("Infants:", nlevels(metadata$Infant_ID), "\n")

# Gene and sample QC
gene_mad <- apply(datExpr, 2, mad, constant = 1, na.rm = TRUE)
zero_mad_genes <- names(gene_mad)[!is.finite(gene_mad) | gene_mad <= 0]

if (length(zero_mad_genes) > 0) {
  datExpr <- datExpr[, !colnames(datExpr) %in% zero_mad_genes, drop = FALSE]
  write.table(zero_mad_genes, file.path(intermediate_dir, "WGCNA_removed_zero_MAD_genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

gsg <- goodSamplesGenes(datExpr, verbose = 3)
removed_samples <- rownames(datExpr)[!gsg$goodSamples]
removed_genes <- colnames(datExpr)[!gsg$goodGenes]

if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  metadata <- metadata[rownames(datExpr), , drop = FALSE]
}

if (nrow(datExpr) < 20) stop("Fewer than 20 samples remain after QC.")
if (ncol(datExpr) < 100) stop("Fewer than 100 genes remain after QC.")

qc_summary <- data.frame(
  Metric = c("Input_samples", "Input_genes", "Removed_zero_MAD_genes", "Removed_bad_samples", "Removed_bad_genes", "Final_samples", "Final_genes", "Unique_infants"),
  Value = c(ncol(expr), nrow(expr), length(zero_mad_genes), length(removed_samples), length(removed_genes), nrow(datExpr), ncol(datExpr), nlevels(droplevels(metadata$Infant_ID)))
)

write.csv(qc_summary, file.path(intermediate_dir, "WGCNA_input_QC.csv"), row.names = FALSE)

# Sample clustering and connectivity
sample_tree <- stats::hclust(stats::dist(datExpr), method = "average")

png(file.path(plot_dir, "WGCNA_sample_clustering.png"), 2200, 1400, res = 180)
par(mar = c(4, 5, 4, 2))
plot(stats::as.dendrogram(sample_tree), main = "Sample clustering", ylab = "Euclidean distance", leaflab = "perpendicular", cex = 0.50)
dev.off()

sample_plot_order <- data.frame(
  Plot_order = seq_along(sample_tree$order),
  Sample = rownames(datExpr)[sample_tree$order],
  Infant_ID = metadata$Infant_ID[sample_tree$order],
  Group = metadata$Group[sample_tree$order],
  Time_point = metadata$Time_point[sample_tree$order]
)

write.csv(sample_plot_order, file.path(intermediate_dir, "WGCNA_sample_dendrogram_order_traits.csv"), row.names = FALSE)

sample_adj <- adjacency(t(datExpr), type = network_type, power = 1, corFnc = cor_type, corOptions = list(maxPOutliers = max_p_outliers))
diag(sample_adj) <- 0
sample_connectivity <- rowSums(sample_adj)
z_connectivity <- as.numeric(scale(sample_connectivity))

sample_qc <- data.frame(
  Sample = rownames(datExpr),
  Infant_ID = metadata$Infant_ID,
  Group = metadata$Group,
  Time_point = metadata$Time_point,
  Connectivity = sample_connectivity,
  Z_connectivity = z_connectivity,
  Flag_below_minus_2.5 = z_connectivity < -2.5
)

write.csv(sample_qc, file.path(intermediate_dir, "WGCNA_sample_connectivity_QC.csv"), row.names = FALSE)
cat("Low-connectivity samples flagged for review:", sum(sample_qc$Flag_below_minus_2.5), "\n")

# Choose the soft-threshold power
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = network_type, corFnc = cor_type, corOptions = list(use = "pairwise.complete.obs", maxPOutliers = max_p_outliers), verbose = 5)

fit_table <- sft$fitIndices
fit_table$Signed_R2 <- -sign(fit_table[, 3]) * fit_table[, 2]
eligible_power <- which(is.finite(fit_table$Signed_R2) & fit_table$Signed_R2 >= scale_free_target)

if (length(eligible_power) > 0) {
  softPower <- fit_table[eligible_power[1], 1]
  power_reason <- paste0("lowest tested power with signed R2 >= ", scale_free_target)
} else {
  softPower <- fallback_power
  power_reason <- "fallback power for a signed network with more than 40 samples"
  warning("No tested power reached the scale-free topology target; using power 12.")
}

write.csv(fit_table, file.path(intermediate_dir, "soft_threshold_results.csv"), row.names = FALSE)
writeLines(c(paste("Chosen softPower =", softPower), paste("Reason =", power_reason)), file.path(intermediate_dir, "chosen_soft_power.txt"))
cat("Chosen soft-threshold power:", softPower, "\n")

png(file.path(plot_dir, "WGCNA_soft_threshold.png"), 2200, 1100, res = 180)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 2))
plot(fit_table[, 1], fit_table$Signed_R2, type = "n", xlab = "Soft-threshold power", ylab = "Signed scale-free topology fit (R2)", main = "Scale independence")
text(fit_table[, 1], fit_table$Signed_R2, labels = fit_table[, 1], cex = 0.8)
abline(h = scale_free_target, col = "red", lty = 2)
abline(v = softPower, col = "blue", lty = 3)
plot(fit_table[, 1], fit_table[, 5], type = "n", xlab = "Soft-threshold power", ylab = "Mean connectivity", main = "Mean connectivity")
text(fit_table[, 1], fit_table[, 5], labels = fit_table[, 1], cex = 0.8)
abline(v = softPower, col = "blue", lty = 3)
dev.off()

# Build the signed network
net <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = network_type,
  TOMType = network_type,
  corType = cor_type,
  maxPOutliers = max_p_outliers,
  maxBlockSize = max_block_size,
  minModuleSize = min_module_size,
  mergeCutHeight = merge_cut_height,
  deepSplit = 2,
  reassignThreshold = 1e-6,
  pamRespectsDendro = TRUE,
  numericLabels = TRUE,
  saveTOMs = FALSE,
  randomSeed = 54321,
  verbose = 5
)

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
names(moduleColors) <- colnames(datExpr)

MEs <- moduleEigengenes(datExpr, colors = moduleColors, excludeGrey = TRUE)$eigengenes
MEs <- orderMEs(MEs)
rownames(MEs) <- rownames(datExpr)

gene_module <- data.frame(Gene = colnames(datExpr), ModuleColor = unname(moduleColors))
module_sizes <- sort(table(moduleColors), decreasing = TRUE)

write.csv(MEs, file.path(results_dir, "module_eigengenes.csv"))
write.csv(gene_module, file.path(results_dir, "gene_module_assignments.csv"), row.names = FALSE)
write.csv(data.frame(ModuleColor = names(module_sizes), Size = as.integer(module_sizes)), file.path(results_dir, "module_sizes.csv"), row.names = FALSE)

cat("\nModule sizes:\n")
print(module_sizes)

for (i in seq_along(net$dendrograms)) {
  block_genes <- net$blockGenes[[i]]
  png(file.path(plot_dir, paste0("WGCNA_gene_dendrogram_block_", i, ".png")), 2200, 1200, res = 180)
  plotDendroAndColors(net$dendrograms[[i]], moduleColors[block_genes], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = paste("Gene dendrogram - block", i))
  dev.off()
}

# Module membership and hub genes
geneModuleMembership <- as.data.frame(bicor(datExpr, MEs, use = "pairwise.complete.obs", maxPOutliers = max_p_outliers))
colnames(geneModuleMembership) <- paste0("kME_", colnames(MEs))

gene_kme <- cbind(gene_module, geneModuleMembership)
own_kme_column <- paste0("kME_ME", gene_kme$ModuleColor)

gene_kme$kME_own_module <- vapply(seq_len(nrow(gene_kme)), function(i) {
  column <- own_kme_column[i]
  if (gene_kme$ModuleColor[i] == "grey" || !column %in% names(gene_kme)) return(NA_real_)
  as.numeric(gene_kme[i, column])
}, numeric(1))

write.csv(gene_kme, file.path(results_dir, "gene_module_membership_kME.csv"), row.names = FALSE)

modules <- sort(setdiff(unique(moduleColors), "grey"))
hub_list <- list()

for (module in modules) {
  module_genes <- gene_kme[gene_kme$ModuleColor == module, , drop = FALSE]
  module_genes <- module_genes[order(module_genes$kME_own_module, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  module_genes$Hub_rank <- seq_len(nrow(module_genes))
  hub_list[[module]] <- head(module_genes, hub_top_n)

  write.table(module_genes$Gene, file.path(module_dir, paste0("module_", module, "_genes.txt")), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.csv(hub_list[[module]], file.path(results_dir, paste0("hub_genes_", module, ".csv")), row.names = FALSE)
}

hub_all <- do.call(rbind, hub_list)
rownames(hub_all) <- NULL
write.csv(hub_all, file.path(results_dir, "hub_genes_all_modules.csv"), row.names = FALSE)

# Longitudinal module models
model_data <- cbind(metadata, MEs)

extract_tests <- function(fit, module_name) {
  overall <- contrast(emmeans(fit, ~ Group, weights = "equal"), method = list(BPD_vs_nonBPD = c(-1, 1)))
  overall <- as.data.frame(summary(overall, infer = TRUE))
  overall$Comparison <- "Overall_BPD_vs_nonBPD"

  group_by_time <- contrast(emmeans(fit, ~ Group | Time_point), method = list(BPD_vs_nonBPD = c(-1, 1)))
  group_by_time <- as.data.frame(summary(group_by_time, infer = TRUE))
  group_by_time$Comparison <- paste0(group_by_time$Time_point, "_BPD_vs_nonBPD")

  time_by_group <- contrast(emmeans(fit, ~ Time_point | Group), method = list(Day14_vs_Cord = c(-1, 1, 0), Day28_vs_Day14 = c(0, -1, 1), Day28_vs_Cord = c(-1, 0, 1)))
  time_by_group <- as.data.frame(summary(time_by_group, infer = TRUE))
  time_by_group$Comparison <- paste0(time_by_group$Group, "_", time_by_group$contrast)

  keep <- c("Comparison", "estimate", "SE", "df", "t.ratio", "p.value")
  results <- rbind(overall[, keep], group_by_time[, keep], time_by_group[, keep])
  results$Module <- module_name
  results$Random_effect_singular <- lme4::isSingular(fit, tol = 1e-4)
  results[, c("Module", keep, "Random_effect_singular")]
}

module_results <- vector("list", ncol(MEs))
model_objects <- vector("list", ncol(MEs))
names(model_objects) <- colnames(MEs)

for (i in seq_len(ncol(MEs))) {
  module_name <- colnames(MEs)[i]
  model_formula <- as.formula(paste0("`", module_name, "` ~ Group * Time_point + (1 | Infant_ID)"))
  fit <- lmer(model_formula, data = model_data, REML = TRUE)
  model_objects[[module_name]] <- fit
  module_results[[i]] <- extract_tests(fit, module_name)
}

all_module_results <- do.call(rbind, module_results)
rownames(all_module_results) <- NULL
all_module_results$FDR <- ave(all_module_results$p.value, all_module_results$Comparison, FUN = function(p) p.adjust(p, method = "BH"))
all_module_results$Significant_FDR05 <- !is.na(all_module_results$FDR) & all_module_results$FDR < fdr_cutoff

comparison_order <- c(
  "Overall_BPD_vs_nonBPD", "Cord_BPD_vs_nonBPD", "Day14_BPD_vs_nonBPD", "Day28_BPD_vs_nonBPD",
  "BPD_Day14_vs_Cord", "BPD_Day28_vs_Day14", "BPD_Day28_vs_Cord",
  "nonBPD_Day14_vs_Cord", "nonBPD_Day28_vs_Day14", "nonBPD_Day28_vs_Cord"
)

all_module_results <- all_module_results[order(match(all_module_results$Comparison, comparison_order), all_module_results$FDR), ]
sig_results <- all_module_results[all_module_results$Significant_FDR05, , drop = FALSE]

write.csv(all_module_results, file.path(results_dir, "module_association_all.csv"), row.names = FALSE)
write.csv(sig_results, file.path(results_dir, "significant_module_associations_FDR05.csv"), row.names = FALSE)

for (comparison in comparison_order) {
  results <- all_module_results[all_module_results$Comparison == comparison, , drop = FALSE]
  write.csv(results, file.path(results_dir, paste0("module_association_", comparison, ".csv")), row.names = FALSE)
}

analysis_summary <- do.call(rbind, lapply(comparison_order, function(comparison) {
  results <- all_module_results[all_module_results$Comparison == comparison, , drop = FALSE]
  data.frame(Analysis = comparison, Tested_modules = nrow(results), Significant_modules = sum(results$Significant_FDR05), Singular_models = sum(results$Random_effect_singular))
}))

rownames(analysis_summary) <- NULL
write.csv(analysis_summary, file.path(results_dir, "WGCNA_analysis_summary.csv"), row.names = FALSE)

# Module association heatmap
effect_mat <- matrix(NA_real_, nrow = ncol(MEs), ncol = length(comparison_order), dimnames = list(colnames(MEs), comparison_order))
fdr_mat <- effect_mat

for (i in seq_len(nrow(all_module_results))) {
  module <- all_module_results$Module[i]
  comparison <- all_module_results$Comparison[i]
  effect_mat[module, comparison] <- all_module_results$estimate[i]
  fdr_mat[module, comparison] <- all_module_results$FDR[i]
}

text_matrix <- matrix("", nrow(effect_mat), ncol(effect_mat), dimnames = dimnames(effect_mat))

for (i in seq_len(nrow(effect_mat))) {
  for (j in seq_len(ncol(effect_mat))) {
    text_matrix[i, j] <- paste0(formatC(effect_mat[i, j], digits = 2, format = "f"), "\nFDR=", formatC(fdr_mat[i, j], digits = 2, format = "g"))
  }
}

heatmap_limit <- max(abs(effect_mat), na.rm = TRUE)

png(file.path(plot_dir, "WGCNA_module_association_heatmap.png"), width = 3000, height = max(1400, 160 + 90 * nrow(effect_mat)), res = 180)
par(mar = c(12, 8, 4, 3))
labeledHeatmap(Matrix = effect_mat, xLabels = colnames(effect_mat), yLabels = rownames(effect_mat), ySymbols = rownames(effect_mat), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = text_matrix, setStdMargins = FALSE, cex.text = 0.45, cex.lab.x = 0.65, cex.lab.y = 0.75, zlim = c(-heatmap_limit, heatmap_limit), main = "Longitudinal module-eigengene associations")
dev.off()

png(file.path(plot_dir, "WGCNA_module_eigengene_dendrogram.png"), 1800, 1200, res = 180)
plotEigengeneNetworks(MEs, "Module eigengene network", marDendro = c(0, 4, 2, 0), marHeatmap = c(3, 4, 2, 2), plotDendrograms = TRUE, xLabelsAngle = 90)
dev.off()

# Save the workspace and summary
save(datExpr, metadata, sft, softPower, power_reason, net, moduleLabels, moduleColors, MEs, gene_kme, model_objects, all_module_results, file = file.path(results_dir, "WGCNA_workspace.RData"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "WGCNA_sessionInfo.txt"))

cat("\nAnalysis summary:\n")
print(analysis_summary, row.names = FALSE)
cat("\nWGCNA completed successfully.\n")
cat("Network:", network_type, cor_type, "with soft power", softPower, "\n")
cat("Samples:", nrow(datExpr), "from", nlevels(metadata$Infant_ID), "infants\n")
cat("Genes:", ncol(datExpr), "\n")
cat("Assigned modules:", length(modules), "(grey is unassigned)\n")
cat("Results saved in:", results_dir, "\n")
cat("Plots saved in:", plot_dir, "\n")

