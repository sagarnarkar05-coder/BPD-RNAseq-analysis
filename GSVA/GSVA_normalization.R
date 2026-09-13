rm(list = ls())
gc()

# Packages and file paths
suppressPackageStartupMessages({
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

project_dir <- "/home/sagar007/RNAseq/GSVA"
counts_file <- "/home/sagar007/RNAseq/counts/gene_counts_reverse.tsv"
metadata_file <- file.path(project_dir, "Sample_metadata.txt")
output_dir <- file.path(project_dir, "inputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

minimum_count <- 10
minimum_samples <- 8

if (!file.exists(counts_file)) stop("Count file not found: ", counts_file)
if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)

# Read counts and metadata
counts_data <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
gene_column <- intersect(c("Geneid", "GeneID", "gene_id"), names(counts_data))[1]
if (is.na(gene_column)) stop("Gene ID column was not found.")

annotation_columns <- intersect(c("Geneid", "GeneID", "gene_id", "Chr", "Start", "End", "Strand", "Length"), names(counts_data))
counts <- as.matrix(counts_data[, !names(counts_data) %in% annotation_columns, drop = FALSE])
rownames(counts) <- counts_data[[gene_column]]
colnames(counts) <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(colnames(counts)))
if (any(abs(counts - round(counts)) > .Machine$double.eps^0.5)) stop("Raw counts must be integers.")
storage.mode(counts) <- "integer"

metadata <- read.delim(metadata_file, check.names = FALSE)
required_columns <- c("Sample", "Infant_ID", "Group", "Time_point")

if (length(setdiff(required_columns, names(metadata))) > 0) stop("Required metadata columns are missing.")
if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)) || anyDuplicated(metadata$Sample)) stop("Duplicate gene or sample IDs found.")
if (any(!is.finite(counts)) || any(counts < 0)) stop("Counts contain invalid values.")

rownames(metadata) <- metadata$Sample
metadata$Infant_ID <- factor(metadata$Infant_ID)
metadata$Group <- factor(metadata$Group, levels = c("nonBPD", "BPD"))
metadata$Time_point <- factor(metadata$Time_point, levels = c("Cord", "Day14", "Day28"))

if (anyNA(metadata[, required_columns])) stop("Metadata contains missing or unexpected values.")
if (!setequal(colnames(counts), rownames(metadata))) stop("Count and metadata sample IDs do not match.")

metadata <- metadata[colnames(counts), , drop = FALSE]

# Filter and variance-stabilize counts
keep_gene <- rowSums(counts >= minimum_count) >= minimum_samples
counts_filtered <- counts[keep_gene, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = metadata, design = ~ Time_point + Group)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
vsd <- vst(dds, blind = FALSE)
expression_ensembl <- assay(vsd)

# Convert Ensembl IDs to gene symbols
ensembl <- rownames(expression_ensembl)
ensembl_clean <- sub("\\..*$", "", ensembl)
symbols <- mapIds(org.Hs.eg.db, keys = ensembl_clean, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")

mapping <- data.frame(Ensembl = ensembl, Ensembl_clean = ensembl_clean, Symbol = unname(symbols))
mapped <- !is.na(mapping$Symbol) & mapping$Symbol != ""
expression_mapped <- expression_ensembl[mapped, , drop = FALSE]
mapped_symbols <- mapping$Symbol[mapped]

expression_symbol <- rowsum(expression_mapped, group = mapped_symbols, reorder = FALSE)
symbol_frequency <- table(mapped_symbols)
expression_symbol <- expression_symbol / as.numeric(symbol_frequency[rownames(expression_symbol)])

constant_symbol <- apply(expression_symbol, 1, function(x) max(x) == min(x))
expression_symbol <- expression_symbol[!constant_symbol, , drop = FALSE]

if (any(!is.finite(expression_symbol)) || anyDuplicated(rownames(expression_symbol))) stop("The final expression matrix is invalid.")

# Save normalized expression and QC
write.csv(expression_ensembl, file.path(output_dir, "normalized_expression_ENSEMBL.csv"))
write.csv(expression_symbol, file.path(output_dir, "normalized_expression_SYMBOL.csv"))
write.csv(mapping, file.path(output_dir, "GSVA_gene_mapping.csv"), row.names = FALSE)
write.table(metadata, file.path(output_dir, "Sample_metadata.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

qc_summary <- data.frame(
  Metric = c("Samples", "Genes_before_filtering", "Genes_after_count_filter", "Mapped_Ensembl_rows", "Unique_symbols_before_constant_filter", "Constant_symbols_removed", "Final_unique_symbols", "Minimum_count", "Minimum_samples"),
  Value = c(ncol(counts), nrow(counts), nrow(counts_filtered), sum(mapped), length(unique(mapped_symbols)), sum(constant_symbol), nrow(expression_symbol), minimum_count, minimum_samples)
)

write.csv(qc_summary, file.path(output_dir, "GSVA_normalization_QC.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "GSVA_normalization_sessionInfo.txt"))

cat("\nSamples:", ncol(counts), "\n")
print(table(metadata$Group, metadata$Time_point))
cat("Genes before filtering:", nrow(counts), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")
cat("Final matrix:", nrow(expression_symbol), "genes x", ncol(expression_symbol), "samples\n")
cat("Results saved in:", output_dir, "\n")
cat("GSVA normalization completed.\n")

