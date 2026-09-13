rm(list = ls())
gc()

# Packages and file paths
suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(pheatmap)
  library(BiocParallel)
})

set.seed(12345)
options(timeout = 7200)

project_dir <- "/home/sagar007/RNAseq/FGSEA"
deg_dir <- "/home/sagar007/RNAseq/DEG/results"
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")
cache_file <- file.path(project_dir, "msigdb_cache", "human_msigdb_pathway_sets.rds")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "Enrichment"), recursive = TRUE, showWarnings = FALSE)

# Analysis settings
min_pathway_size <- 15
max_pathway_size <- 500
fdr_cutoff <- 0.05
figure_nes_cutoff <- 1.5
nperm_simple <- 10000
bpparam <- SerialParam(progressbar = FALSE)

positive_colour <- "#E64B35"
negative_colour <- "#3C5488"
neutral_colour <- "#F7F7F7"

comparisons <- data.frame(
  ID = c("Overall_BPD_vs_nonBPD", "Cord_BPD_vs_nonBPD", "Day14_BPD_vs_nonBPD", "Day28_BPD_vs_nonBPD", "BPD_Day14_vs_Cord", "BPD_Day28_vs_Day14", "BPD_Day28_vs_Cord", "nonBPD_Day14_vs_Cord", "nonBPD_Day28_vs_Day14", "nonBPD_Day28_vs_Cord"),
  Label = c("Overall BPD vs nonBPD", "Cord BPD vs nonBPD", "Day14 BPD vs nonBPD", "Day28 BPD vs nonBPD", "BPD Day14 vs Cord", "BPD Day28 vs Day14", "BPD Day28 vs Cord", "nonBPD Day14 vs Cord", "nonBPD Day28 vs Day14", "nonBPD Day28 vs Cord"),
  Rank_file = c("DEG_Overall_BPD_vs_nonBPD_FGSEA_ranked_statistics.csv", "DEG_Cord_BPD_vs_nonBPD_FGSEA_ranked_statistics.csv", "DEG_Day14_BPD_vs_nonBPD_FGSEA_ranked_statistics.csv", "DEG_Day28_BPD_vs_nonBPD_FGSEA_ranked_statistics.csv", "DEG_BPD_paired_Day14_vs_Cord_FGSEA_ranked_statistics.csv", "DEG_BPD_paired_Day28_vs_Day14_FGSEA_ranked_statistics.csv", "DEG_BPD_paired_Day28_vs_Cord_FGSEA_ranked_statistics.csv", "DEG_nonBPD_paired_Day14_vs_Cord_FGSEA_ranked_statistics.csv", "DEG_nonBPD_paired_Day28_vs_Day14_FGSEA_ranked_statistics.csv", "DEG_nonBPD_paired_Day28_vs_Cord_FGSEA_ranked_statistics.csv")
)

if (any(!file.exists(file.path(deg_dir, comparisons$Rank_file)))) stop("One or more DEG ranking files are missing.")

# Load MSigDB pathway collections
make_pathways <- function(collection, subcollection = NULL) {
  arguments <- list(db_species = "HS", species = "Homo sapiens", collection = collection)
  if (!is.null(subcollection)) arguments$subcollection <- subcollection
  genes <- do.call(msigdbr, arguments)
  genes <- genes[!is.na(genes$gs_name) & !is.na(genes$gene_symbol), c("gs_name", "gene_symbol")]
  genes <- unique(genes)
  split(genes$gene_symbol, genes$gs_name)
}

if (file.exists(cache_file)) {
  pathway_sets <- readRDS(cache_file)
  cat("Using cached MSigDB pathway sets.\n")
} else {
  cat("Downloading MSigDB pathway sets...\n")
  pathway_sets <- list(
    Hallmark = make_pathways("H"),
    Reactome = make_pathways("C2", "CP:REACTOME"),
    GO_BP = make_pathways("C5", "GO:BP"),
    KEGG = make_pathways("C2", "CP:KEGG_LEGACY")
  )
  saveRDS(pathway_sets, cache_file)
}

database_names <- names(pathway_sets)

for (database in database_names) {
  dir.create(file.path(results_dir, database), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(figures_dir, database), recursive = TRUE, showWarnings = FALSE)
}

pathway_summary <- data.frame(Database = database_names, Gene_sets = lengths(pathway_sets))
write.csv(pathway_summary, file.path(results_dir, "MSigDB_pathway_summary.csv"), row.names = FALSE)

# Ranking preparation
read_ranking <- function(file, comparison_id) {
  ranking <- read.csv(file, check.names = FALSE)
  if (!all(c("Symbol", "stat") %in% names(ranking))) stop("Ranking file is missing Symbol or stat.")

  input_rows <- nrow(ranking)
  missing_symbol <- sum(is.na(ranking$Symbol) | trimws(ranking$Symbol) == "")
  invalid_stat <- sum(!is.finite(ranking$stat))

  ranking <- ranking %>%
    transmute(Symbol = trimws(as.character(Symbol)), stat = as.numeric(stat)) %>%
    filter(!is.na(Symbol), Symbol != "", is.finite(stat))

  duplicate_symbols <- sum(duplicated(ranking$Symbol))

  ranking <- ranking %>%
    group_by(Symbol) %>%
    arrange(desc(abs(stat)), desc(stat), .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup()

  stats <- setNames(ranking$stat, ranking$Symbol)
  stats <- sort(stats, decreasing = TRUE)

  if (length(stats) < 1000) stop(comparison_id, ": fewer than 1,000 ranked genes remain.")
  if (all(stats >= 0) || all(stats <= 0)) stop(comparison_id, ": the ranking is one-sided.")

  qc <- data.frame(
    Comparison_ID = comparison_id,
    Input_rows = input_rows,
    Missing_symbol_rows = missing_symbol,
    Invalid_stat_rows = invalid_stat,
    Duplicate_symbol_rows_removed = duplicate_symbols,
    Ranked_unique_symbols = length(stats),
    Positive_statistics = sum(stats > 0),
    Negative_statistics = sum(stats < 0),
    Zero_statistics = sum(stats == 0),
    Tie_fraction = 1 - length(unique(stats)) / length(stats),
    Min_stat = min(stats),
    Max_stat = max(stats)
  )

  list(stats = stats, qc = qc)
}

# Run FGSEA for one comparison and database
run_fgsea <- function(stats, pathways, comparison_id, comparison_label, database) {
  overlap <- vapply(pathways, function(genes) sum(genes %in% names(stats)), integer(1))
  overlap_table <- data.frame(Comparison_ID = comparison_id, Comparison = comparison_label, Database = database, Pathway = names(pathways), Original_size = lengths(pathways), Overlap_size = overlap)
  overlap_table$Overlap_fraction <- overlap_table$Overlap_size / overlap_table$Original_size

  results <- fgseaMultilevel(pathways = pathways, stats = stats, minSize = min_pathway_size, maxSize = max_pathway_size, eps = 0, nPermSimple = nperm_simple, scoreType = "std", BPPARAM = bpparam)
  results <- as.data.frame(results)
  results$leadingEdge <- vapply(results$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  results$Comparison_ID <- comparison_id
  results$Comparison <- comparison_label
  results$Database <- database
  results <- results[, c("Comparison_ID", "Comparison", "Database", "pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge")]
  results <- results[order(results$padj, -abs(results$NES)), ]

  significant <- results[!is.na(results$padj) & results$padj < fdr_cutoff, , drop = FALSE]
  main_pathways <- significant[0, ]

  if (nrow(significant) > 0) {
    collapsed <- collapsePathways(as.data.table(significant), pathways, stats, pval.threshold = fdr_cutoff)
    main_pathways <- significant[significant$pathway %in% collapsed$mainPathways, , drop = FALSE]
  }

  prefix <- paste0("FGSEA_", comparison_id, "_", database)
  write.csv(results, file.path(results_dir, database, paste0(prefix, "_full.csv")), row.names = FALSE)
  write.csv(significant, file.path(results_dir, database, paste0(prefix, "_significant.csv")), row.names = FALSE)
  write.csv(main_pathways, file.path(results_dir, database, paste0(prefix, "_main_pathways.csv")), row.names = FALSE)

  list(full = results, significant = significant, main = main_pathways, overlap = overlap_table)
}

# Run all comparisons
all_results <- list()
all_rankings <- list()
ranking_qc <- list()
overlap_tables <- list()
summary_tables <- list()

for (i in seq_len(nrow(comparisons))) {
  comparison <- comparisons[i, ]
  ranking <- read_ranking(file.path(deg_dir, comparison$Rank_file), comparison$ID)
  all_rankings[[comparison$ID]] <- ranking$stats
  ranking_qc[[i]] <- ranking$qc

  cat("\n", comparison$Label, " - ", length(ranking$stats), " ranked genes\n", sep = "")

  for (database in database_names) {
    cat("  ", database, "\n", sep = "")
    analysis <- run_fgsea(ranking$stats, pathway_sets[[database]], comparison$ID, comparison$Label, database)
    key <- paste(comparison$ID, database, sep = "__")
    all_results[[key]] <- analysis$full
    overlap_tables[[key]] <- analysis$overlap
    summary_tables[[key]] <- data.frame(Comparison = comparison$Label, Database = database, Ranked_genes = length(ranking$stats), Tested_pathways = nrow(analysis$full), Significant_pathways = nrow(analysis$significant), Main_nonredundant_pathways = nrow(analysis$main))
  }
}

all_results <- bind_rows(all_results)
analysis_summary <- bind_rows(summary_tables)
ranking_qc <- bind_rows(ranking_qc)
pathway_overlap <- bind_rows(overlap_tables)

write.csv(all_results, file.path(results_dir, "FGSEA_all_results.csv"), row.names = FALSE)
write.csv(analysis_summary, file.path(results_dir, "FGSEA_analysis_summary.csv"), row.names = FALSE)
write.csv(ranking_qc, file.path(results_dir, "FGSEA_ranking_QC.csv"), row.names = FALSE)
write.csv(pathway_overlap, file.path(results_dir, "FGSEA_pathway_overlap_audit.csv"), row.names = FALSE)

# Prepare figure data
clean_pathway_name <- function(pathway) {
  pathway %>%
    str_remove("^HALLMARK_") %>%
    str_remove("^REACTOME_") %>%
    str_remove("^KEGG_") %>%
    str_remove("^GOBP_") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

figure_data <- all_results %>%
  filter(!is.na(padj), padj < fdr_cutoff, !is.na(NES), abs(NES) >= figure_nes_cutoff) %>%
  mutate(Pathway = clean_pathway_name(pathway), minus_log10_FDR = -log10(pmax(padj, .Machine$double.xmin)))

write.csv(figure_data, file.path(results_dir, "FGSEA_results_used_for_figures.csv"), row.names = FALSE)

plot_theme <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = "grey30"), axis.text = element_text(colour = "black"))

# Per-comparison NES bar plots
make_barplot <- function(data, title, file, top_n = 8) {
  positive <- data %>% filter(NES > 0) %>% slice_max(NES, n = top_n, with_ties = FALSE)
  negative <- data %>% filter(NES < 0) %>% slice_min(NES, n = top_n, with_ties = FALSE)
  plot_data <- bind_rows(positive, negative) %>% distinct(pathway, .keep_all = TRUE) %>% mutate(Pathway = fct_reorder(Pathway, NES), Direction = ifelse(NES > 0, "Positive", "Negative"))
  if (nrow(plot_data) == 0) return(invisible(NULL))

  p <- ggplot(plot_data, aes(NES, Pathway, fill = Direction)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, colour = "grey35") +
    scale_fill_manual(values = c(Negative = negative_colour, Positive = positive_colour)) +
    labs(title = title, subtitle = "FDR < 0.05 and |NES| >= 1.5", x = "Normalized enrichment score", y = NULL, fill = NULL) +
    plot_theme +
    theme(legend.position = "top")

  ggsave(file, p, width = 9, height = 7, dpi = 300, bg = "white")
}

# Cross-comparison dot plots
make_dotplot <- function(data, title, file, top_n = 20) {
  top_pathways <- data %>% group_by(Pathway) %>% summarise(Best_NES = max(abs(NES)), .groups = "drop") %>% slice_max(Best_NES, n = top_n, with_ties = FALSE)
  plot_data <- data %>% filter(Pathway %in% top_pathways$Pathway) %>% mutate(Pathway = factor(Pathway, levels = rev(top_pathways$Pathway)))
  if (nrow(plot_data) == 0) return(invisible(NULL))

  p <- ggplot(plot_data, aes(Comparison, Pathway, size = minus_log10_FDR, colour = NES)) +
    geom_point(alpha = 0.90) +
    scale_colour_gradient2(low = negative_colour, mid = neutral_colour, high = positive_colour, midpoint = 0) +
    scale_size_continuous(range = c(2.5, 8)) +
    labs(title = title, subtitle = "FDR < 0.05 and |NES| >= 1.5", x = NULL, y = NULL, size = expression(-log[10](FDR)), colour = "NES") +
    plot_theme +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  ggsave(file, p, width = 11, height = 8, dpi = 300, bg = "white")
}

# Cross-comparison NES heatmaps
make_heatmap <- function(data, title, file, top_n = 20) {
  top_pathways <- data %>% group_by(Pathway) %>% summarise(Significant_comparisons = n_distinct(Comparison), Best_NES = max(abs(NES)), .groups = "drop") %>% arrange(desc(Significant_comparisons), desc(Best_NES)) %>% slice_head(n = top_n)
  matrix_data <- data %>% filter(Pathway %in% top_pathways$Pathway) %>% group_by(Pathway, Comparison) %>% slice_max(abs(NES), n = 1, with_ties = FALSE) %>% ungroup() %>% select(Pathway, Comparison, NES) %>% pivot_wider(names_from = Comparison, values_from = NES, values_fill = 0)
  if (nrow(matrix_data) == 0) return(invisible(NULL))

  matrix <- as.data.frame(matrix_data)
  rownames(matrix) <- matrix$Pathway
  matrix$Pathway <- NULL
  matrix <- as.matrix(matrix)

  png(file, width = 2600, height = 2200, res = 300)
  pheatmap(matrix, color = colorRampPalette(c(negative_colour, "white", positive_colour))(100), cluster_rows = nrow(matrix) > 1, cluster_cols = FALSE, border_color = NA, na_col = "grey95", fontsize_row = 9, fontsize_col = 9, angle_col = 45, main = title)
  dev.off()
}

# Generate the main figures
comparison_groups <- list(
  Disease_timepoints = c("Cord BPD vs nonBPD", "Day14 BPD vs nonBPD", "Day28 BPD vs nonBPD"),
  BPD_progression = c("BPD Day14 vs Cord", "BPD Day28 vs Day14", "BPD Day28 vs Cord"),
  nonBPD_progression = c("nonBPD Day14 vs Cord", "nonBPD Day28 vs Day14", "nonBPD Day28 vs Cord")
)

for (database in database_names) {
  database_data <- figure_data %>% filter(Database == database)

  for (i in seq_len(nrow(comparisons))) {
    comparison <- comparisons[i, ]
    plot_data <- database_data %>% filter(Comparison_ID == comparison$ID)
    make_barplot(plot_data, paste(database, "-", comparison$Label), file.path(figures_dir, database, paste0(comparison$ID, "_barplot.png")))
  }

  for (group_name in names(comparison_groups)) {
    plot_data <- database_data %>% filter(Comparison %in% comparison_groups[[group_name]])
    title <- paste(database, "-", str_replace_all(group_name, "_", " "))
    make_dotplot(plot_data, title, file.path(figures_dir, database, paste0(database, "_", group_name, "_dotplot.png")))
    make_heatmap(plot_data, title, file.path(figures_dir, database, paste0(database, "_", group_name, "_heatmap.png")))
  }
}

# Enrichment curves for the strongest Hallmark pathways
hallmark_results <- all_results %>% filter(Database == "Hallmark", !is.na(padj), padj < fdr_cutoff)

for (i in seq_len(nrow(comparisons))) {
  comparison <- comparisons[i, ]
  top_pathways <- hallmark_results %>% filter(Comparison_ID == comparison$ID) %>% arrange(padj, desc(abs(NES))) %>% slice_head(n = 5)
  stats <- all_rankings[[comparison$ID]]

  for (pathway in top_pathways$pathway) {
    p <- plotEnrichment(pathway_sets$Hallmark[[pathway]], stats) +
      labs(title = clean_pathway_name(pathway), subtitle = comparison$Label, x = "Rank in ordered gene list", y = "Enrichment score") +
      plot_theme

    file_name <- paste0(comparison$ID, "__", str_replace_all(pathway, "[^A-Za-z0-9_]+", "_"), ".png")
    ggsave(file.path(figures_dir, "Enrichment", file_name), p, width = 8, height = 5.5, dpi = 300, bg = "white")
  }
}

# Save run information
saveRDS(pathway_sets, file.path(results_dir, "MSigDB_pathway_sets_used.rds"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "R_sessionInfo.txt"))
writeLines(c(paste("Ranking: signed DEG model statistic from all tested genes"), paste("FDR cutoff:", fdr_cutoff), paste("Figure NES cutoff:", figure_nes_cutoff), paste("Pathway size:", min_pathway_size, "to", max_pathway_size), paste("fgseaMultilevel eps: 0"), paste("nPermSimple:", nperm_simple), paste("msigdbr version:", packageVersion("msigdbr")), paste("fgsea version:", packageVersion("fgsea"))), file.path(results_dir, "FGSEA_run_parameters.txt"))

cat("\nAnalysis summary:\n")
print(analysis_summary)
cat("\nResults saved in:", results_dir, "\n")
cat("Figures saved in:", figures_dir, "\n")
cat("FGSEA analysis completed.\n")

