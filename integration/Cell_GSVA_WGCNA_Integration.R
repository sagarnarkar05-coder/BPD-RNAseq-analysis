rm(list = ls())
gc()

# Packages and file paths
suppressPackageStartupMessages({
  library(nlme)
  library(ggplot2)
  library(pheatmap)
})

base_dir <- "/home/sagar007/RNAseq"
out_dir <- file.path(base_dir, "Stats", "cell_gsva_wgcna_integration")
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

metadata_candidates <- c(
  file.path(base_dir, "GSVA", "inputs", "Sample_metadata.txt"),
  file.path(base_dir, "WGCNA", "Sample_metadata.txt")
)
metadata_found <- metadata_candidates[file.exists(metadata_candidates)]
metadata_file <- if (length(metadata_found)) metadata_found[1] else NA_character_
gsva_file <- file.path(base_dir, "GSVA", "results", "GSVA_Hallmark_scores_all_samples.csv")
cell_file <- file.path(base_dir, "CELL", "CIBERSORTx_Results.txt")
wgcna_file <- file.path(base_dir, "WGCNA", "results", "module_eigengenes.csv")

files <- c(metadata_file, gsva_file, cell_file, wgcna_file)
if (is.na(metadata_file) || any(!file.exists(files))) {
  stop("One or more input files are missing.")
}

fdr_cutoff <- 0.05
cibersort_p_cutoff <- 0.05
min_nonzero_fraction <- 0.20

# Helper functions
zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

logit_fraction <- function(x) {
  n <- sum(is.finite(x))
  qlogis((x * (n - 1) + 0.5) / n)
}

fit_pair <- function(y, x, metadata, outcome, predictor, analysis) {
  dat <- data.frame(
    y = zscore(as.numeric(y)),
    x = zscore(as.numeric(x)),
    Infant_ID = metadata$Infant_ID,
    Group = metadata$Group,
    Time_point = metadata$Time_point
  )
  dat <- dat[complete.cases(dat), ]

  if (nrow(dat) < 15 || length(unique(dat$Infant_ID)) < 6 ||
      sd(dat$y) == 0 || sd(dat$x) == 0) {
    return(data.frame(
      Analysis = analysis, Outcome = outcome, Predictor = predictor,
      N_samples = nrow(dat), N_infants = length(unique(dat$Infant_ID)),
      Estimate = NA, SE = NA, DF = NA, t_value = NA, P_value = NA,
      Random_intercept_SD = NA, Residual_SD = NA,
      Model_status = "Insufficient data"
    ))
  }

  fit <- tryCatch(
    lme(
      y ~ x + Group * Time_point,
      random = ~1 | Infant_ID,
      data = dat,
      method = "REML",
      control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.frame(
      Analysis = analysis, Outcome = outcome, Predictor = predictor,
      N_samples = nrow(dat), N_infants = length(unique(dat$Infant_ID)),
      Estimate = NA, SE = NA, DF = NA, t_value = NA, P_value = NA,
      Random_intercept_SD = NA, Residual_SD = NA,
      Model_status = paste("Failed:", fit$message)
    ))
  }

  tt <- summary(fit)$tTable["x", ]
  vc <- VarCorr(fit)

  data.frame(
    Analysis = analysis,
    Outcome = outcome,
    Predictor = predictor,
    N_samples = nrow(dat),
    N_infants = length(unique(dat$Infant_ID)),
    Estimate = tt["Value"],
    SE = tt["Std.Error"],
    DF = tt["DF"],
    t_value = tt["t-value"],
    P_value = tt["p-value"],
    Random_intercept_SD = as.numeric(vc[1, "StdDev"]),
    Residual_SD = as.numeric(vc[nrow(vc), "StdDev"]),
    Model_status = "Success"
  )
}

run_models <- function(outcomes, predictors, metadata, analysis) {
  ans <- vector("list", nrow(outcomes) * nrow(predictors))
  k <- 1

  for (i in seq_len(nrow(outcomes))) {
    for (j in seq_len(nrow(predictors))) {
      ans[[k]] <- fit_pair(
        outcomes[i, ], predictors[j, ], metadata,
        rownames(outcomes)[i], rownames(predictors)[j], analysis
      )
      k <- k + 1
    }
  }

  ans <- do.call(rbind, ans)
  ans$FDR <- p.adjust(ans$P_value, method = "BH")
  ans$Significant_FDR_0.05 <- !is.na(ans$FDR) & ans$FDR < fdr_cutoff
  ans[order(ans$FDR, ans$P_value), ]
}

save_results <- function(x, prefix) {
  write.csv(
    x,
    file.path(table_dir, paste0(prefix, "_mixed_model_all_results.csv")),
    row.names = FALSE
  )
  write.csv(
    x[x$Significant_FDR_0.05, ],
    file.path(table_dir, paste0(prefix, "_mixed_model_significant.csv")),
    row.names = FALSE
  )
}

plot_heatmap <- function(x, y, filename, title) {
  r <- cor(t(x), t(y), method = "spearman", use = "pairwise.complete.obs")
  write.csv(r, file.path(table_dir, paste0(filename, "_spearman_correlations.csv")))

  heat_colours <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)
  heat_breaks <- seq(-1, 1, length.out = 102)

  png(file.path(figure_dir, paste0(filename, "_spearman_heatmap.png")), width = max(1800, ncol(r) * 140), height = max(1400, nrow(r) * 105), res = 200)
  pheatmap(r, color = heat_colours, breaks = heat_breaks, border_color = NA, fontsize = 9, fontsize_row = 8, fontsize_col = 8, angle_col = 45, main = title)
  dev.off()
}

plot_top <- function(results, outcomes, predictors, metadata, prefix, n = 6) {
  top <- head(results[results$Model_status == "Success", ], n)

  for (i in seq_len(nrow(top))) {
    out <- top$Outcome[i]
    pred <- top$Predictor[i]
    d <- data.frame(
      Outcome = as.numeric(outcomes[out, rownames(metadata)]),
      Predictor = as.numeric(predictors[pred, rownames(metadata)]),
      Group = metadata$Group,
      Time_point = metadata$Time_point
    )

    p <- ggplot(d, aes(Predictor, Outcome, colour = Group, shape = Time_point)) +
      geom_point(size = 2.8, alpha = 0.85) +
      geom_smooth(aes(x = Predictor, y = Outcome), inherit.aes = FALSE, data = d, method = "lm", formula = y ~ x, se = TRUE, colour = "grey35", fill = "grey80", linewidth = 0.8) +
      scale_colour_manual(values = c(nonBPD = "#0072B2", BPD = "#D55E00")) +
      labs(title = paste(out, "vs", pred), subtitle = paste0("Mixed-model BH FDR = ", signif(top$FDR[i], 3)), x = pred, y = out, colour = "Group", shape = "Time point") +
      theme_classic(base_size = 12) +
      theme(legend.position = "top", plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = "grey30"))

    name <- gsub("[^A-Za-z0-9]+", "_", paste(prefix, out, "vs", pred))
    ggsave(file.path(figure_dir, paste0(name, ".png")), p, width = 8, height = 6, dpi = 300, bg = "white")
  }
}

# Read the input files
metadata <- read.delim(metadata_file, check.names = FALSE)
required <- c("Sample", "Infant_ID", "Group", "Time_point")
if (length(setdiff(required, names(metadata))) > 0 || anyDuplicated(metadata$Sample)) {
  stop("Metadata columns or sample IDs are invalid.")
}
rownames(metadata) <- metadata$Sample
metadata$Infant_ID <- factor(metadata$Infant_ID)
metadata$Group <- factor(metadata$Group, levels = c("nonBPD", "BPD"))
metadata$Time_point <- factor(metadata$Time_point, levels = c("Cord", "Day14", "Day28"))
if (anyNA(metadata$Group) || anyNA(metadata$Time_point)) stop("Invalid metadata values.")

gsva <- as.matrix(read.csv(gsva_file, row.names = 1, check.names = FALSE))
storage.mode(gsva) <- "double"

cell_data <- read.delim(cell_file, row.names = 1, check.names = FALSE)
p_col <- intersect(c("P-value", "P.value"), names(cell_data))[1]
if (is.na(p_col)) stop("CIBERSORTx P-value column is missing.")
cibersort_p <- setNames(as.numeric(cell_data[[p_col]]), rownames(cell_data))
metric_cols <- intersect(
  c("P-value", "P.value", "Correlation", "RMSE", "Absolute score"),
  names(cell_data)
)
cells <- t(as.matrix(cell_data[, setdiff(names(cell_data), metric_cols), drop = FALSE]))
storage.mode(cells) <- "double"

me <- as.matrix(read.csv(wgcna_file, row.names = 1, check.names = FALSE))
storage.mode(me) <- "double"
if (all(rownames(me) %in% rownames(metadata))) {
  modules <- t(me)
} else if (all(colnames(me) %in% rownames(metadata))) {
  modules <- me
} else {
  stop("Could not determine WGCNA matrix orientation.")
}

# Match samples and apply QC filters
common <- Reduce(intersect, list(
  rownames(metadata), colnames(gsva), colnames(cells), colnames(modules),
  names(cibersort_p)[cibersort_p < cibersort_p_cutoff]
))
common <- rownames(metadata)[rownames(metadata) %in% common]
if (length(common) < 15) stop("Too few common QC-passing samples.")

metadata <- metadata[common, ]
gsva <- gsva[, common, drop = FALSE]
cells <- cells[, common, drop = FALSE]
modules <- modules[, common, drop = FALSE]

if (any(!is.finite(gsva)) || any(!is.finite(cells)) || any(!is.finite(modules))) {
  stop("An input matrix contains missing or non-finite values.")
}
if (any(cells < 0 | cells > 1)) stop("Cell fractions are outside 0-1.")

nonzero <- rowSums(cells > 0)
unique_n <- apply(cells, 1, function(x) length(unique(x)))
eligible <- nonzero >= ceiling(min_nonzero_fraction * length(common)) & unique_n >= 5

cell_qc <- data.frame(Cell_type = rownames(cells), Nonzero_samples = nonzero, Unique_values = unique_n, Eligible_for_models = eligible)
write.csv(cell_qc, file.path(table_dir, "CIBERSORTx_cell_type_model_QC.csv"), row.names = FALSE)

cells_raw <- cells[eligible, , drop = FALSE]
cells_logit <- t(apply(cells_raw, 1, logit_fraction))
dimnames(cells_logit) <- dimnames(cells_raw)

sample_qc <- data.frame(Sample = rownames(cell_data), CIBERSORTx_P_value = cibersort_p[rownames(cell_data)], Pass_P_lt_0.05 = cibersort_p[rownames(cell_data)] < cibersort_p_cutoff, Used_in_integration = rownames(cell_data) %in% common)
write.csv(sample_qc, file.path(table_dir, "CIBERSORTx_sample_integration_QC.csv"), row.names = FALSE)

cat("Samples:", length(common), "\n")
cat("Cell types:", nrow(cells_logit), "\n")
cat("GSVA pathways:", nrow(gsva), "\n")
cat("WGCNA modules:", nrow(modules), "\n")

# Run the mixed-effects models
cat("\nRunning Cell-GSVA models...\n")
cell_gsva <- run_models(gsva, cells_logit, metadata, "Cell-GSVA")
save_results(cell_gsva, "Cell_GSVA")

cat("Running Cell-WGCNA models...\n")
cell_wgcna <- run_models(modules, cells_logit, metadata, "Cell-WGCNA")
save_results(cell_wgcna, "Cell_WGCNA")

cat("Running WGCNA-GSVA models...\n")
wgcna_gsva <- run_models(gsva, modules, metadata, "WGCNA-GSVA")
save_results(wgcna_gsva, "WGCNA_GSVA")

# Select features for the main figures
priority_cells <- intersect(c("Neutrophils", "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2", "NK cells activated", "T cells CD8", "B cells naive"), rownames(cells_raw))
priority_pathways <- grep("INFLAMMATORY_RESPONSE|TNFA_SIGNALING|MTORC1|OXIDATIVE_PHOSPHORYLATION|GLYCOLYSIS|INTERFERON|COMPLEMENT|HYPOXIA|IL6_JAK_STAT3", rownames(gsva), value = TRUE)
priority_modules <- intersect(c("MEbrown", "MElightcyan", "MEmidnightblue", "MElightyellow", "MEblack", "MEtan", "MEgreen", "MEdarkturquoise"), rownames(modules))

# Correlation heatmaps
plot_heatmap(cells_raw[priority_cells, ], gsva[priority_pathways, ], "Priority_Cell_GSVA", "Cell fractions and GSVA pathways")
plot_heatmap(cells_raw[priority_cells, ], modules[priority_modules, ], "Priority_Cell_WGCNA", "Cell fractions and WGCNA modules")
plot_heatmap(modules[priority_modules, ], gsva[priority_pathways, ], "Priority_WGCNA_GSVA", "WGCNA modules and GSVA pathways")

# Plots of the strongest associations
plot_top(cell_gsva, gsva, cells_raw, metadata, "Cell_GSVA")
plot_top(cell_wgcna, modules, cells_raw, metadata, "Cell_WGCNA")
plot_top(wgcna_gsva, gsva, modules, metadata, "WGCNA_GSVA")

# Save the analysis summary
summary_table <- data.frame(
  Analysis = c("Cell-GSVA", "Cell-WGCNA", "WGCNA-GSVA"),
  Total_tests = c(nrow(cell_gsva), nrow(cell_wgcna), nrow(wgcna_gsva)),
  Successful_models = c(
    sum(cell_gsva$Model_status == "Success"),
    sum(cell_wgcna$Model_status == "Success"),
    sum(wgcna_gsva$Model_status == "Success")
  ),
  Significant_FDR_0.05 = c(
    sum(cell_gsva$Significant_FDR_0.05),
    sum(cell_wgcna$Significant_FDR_0.05),
    sum(wgcna_gsva$Significant_FDR_0.05)
  )
)

write.csv(summary_table, file.path(table_dir, "Integration_analysis_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "Integration_sessionInfo.txt"))

cat("\nIntegration completed.\n")
print(summary_table, row.names = FALSE)
cat("Results:", out_dir, "\n")

