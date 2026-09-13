suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# File paths
project_dir <- "/home/sagar007/RNAseq/CELL"
counts_file <- "/home/sagar007/RNAseq/counts/gene_counts_reverse.tsv"
gtf_file <- "/home/sagar007/RNAseq/reference/gencode.v44.annotation.gtf"
metadata_file <- file.path(project_dir, "Sample_metadata.txt")
mixture_file <- file.path(project_dir, "CIBERSORTx_mixture_TPM.txt")
mapping_file <- file.path(project_dir, "CIBERSORTx_gene_mapping.csv")
qc_file <- file.path(project_dir, "CIBERSORTx_input_QC.csv")
length_cache <- file.path(project_dir, "gencode_v44_gene_exonic_lengths.rds")

for (file in c(counts_file, gtf_file, metadata_file)) {
  if (!file.exists(file)) stop("File not found: ", file)
}

# Read counts and metadata
counts_df <- read.delim(counts_file, check.names = FALSE)
gene_col <- intersect(c("GeneID", "Geneid", "gene_id"), names(counts_df))
if (length(gene_col) != 1) stop("Could not identify exactly one gene-ID column.")

counts <- as.matrix(counts_df[, setdiff(names(counts_df), gene_col), drop = FALSE])
storage.mode(counts) <- "numeric"
rownames(counts) <- counts_df[[gene_col]]

if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) stop("Duplicate gene or sample IDs found.")
if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) stop("Count matrix contains invalid values.")

meta <- read.delim(metadata_file, check.names = FALSE)
required <- c("Sample", "Infant_ID", "Time_point", "Group")
if (length(setdiff(required, names(meta)))) stop("Required metadata columns are missing.")
if (anyDuplicated(meta$Sample)) stop("Duplicate sample IDs found in metadata.")

rownames(meta) <- meta$Sample
if (!setequal(colnames(counts), rownames(meta))) stop("Counts and metadata samples do not match.")
meta <- meta[colnames(counts), , drop = FALSE]

# Calculate non-overlapping exonic length for each gene
if (file.exists(length_cache)) {
  cat("Loading cached gene lengths:\n ", length_cache, "\n")
  gene_length_bp <- readRDS(length_cache)
} else {
  cat("Importing exons from:\n ", gtf_file, "\n")
  exons <- import(gtf_file, format = "gtf", feature.type = "exon")
  if (!"gene_id" %in% names(mcols(exons))) stop("The GTF has no gene_id attribute.")

  exons <- exons[!is.na(exons$gene_id) & nzchar(exons$gene_id)]
  exons$gene_id_clean <- sub("\\..*$", "", exons$gene_id)
  exons_by_gene <- split(exons, exons$gene_id_clean)
  gene_length_bp <- vapply(exons_by_gene, function(x) sum(width(reduce(x))), numeric(1))
  gene_length_bp <- gene_length_bp[is.finite(gene_length_bp) & gene_length_bp > 0]
  saveRDS(gene_length_bp, length_cache)
  rm(exons, exons_by_gene)
  gc()
}

# Convert counts to TPM
ensembl <- rownames(counts)
ensembl_clean <- sub("\\..*$", "", ensembl)
length_bp <- unname(gene_length_bp[ensembl_clean])
length_found <- is.finite(length_bp) & length_bp > 0
if (mean(length_found) < 0.90) stop("Gene lengths found for less than 90% of count rows.")

counts_length <- counts[length_found, , drop = FALSE]
rpk <- sweep(counts_length, 1, length_bp[length_found] / 1000, "/")
rpk_sum <- colSums(rpk)
if (any(!is.finite(rpk_sum)) || any(rpk_sum <= 0)) stop("Invalid RPK library total.")
tpm <- sweep(rpk, 2, rpk_sum / 1e6, "/")

# Map Ensembl IDs to gene symbols
symbol <- mapIds(org.Hs.eg.db, keys = ensembl_clean[length_found], keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
mapping <- data.frame(Ensembl = ensembl[length_found], Ensembl_clean = ensembl_clean[length_found], Exonic_length_bp = length_bp[length_found], Symbol = unname(symbol))

symbol_found <- !is.na(mapping$Symbol) & nzchar(mapping$Symbol)
tpm_symbol <- rowsum(tpm[symbol_found, , drop = FALSE], group = mapping$Symbol[symbol_found], reorder = FALSE)
tpm_symbol <- tpm_symbol[rowSums(tpm_symbol) > 0, , drop = FALSE]

if (anyNA(tpm_symbol) || any(!is.finite(tpm_symbol)) || any(tpm_symbol < 0)) stop("Final TPM matrix contains invalid values.")
if (anyDuplicated(rownames(tpm_symbol))) stop("Duplicate gene symbols remain.")

# Save CIBERSORTx mixture and QC files
mixture <- data.frame(GeneSymbol = rownames(tpm_symbol), tpm_symbol, check.names = FALSE)
write.table(mixture, mixture_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(mapping, mapping_file, row.names = FALSE)

input_qc <- data.frame(
  Metric = c("Samples", "Count_matrix_genes", "Genes_with_GTF_length", "Genes_without_GTF_length", "Genes_with_symbol", "Genes_without_symbol", "Final_unique_symbols", "Minimum_TPM", "Maximum_TPM"),
  Value = c(ncol(counts), nrow(counts), sum(length_found), sum(!length_found), sum(symbol_found), sum(!symbol_found), nrow(tpm_symbol), min(tpm_symbol), max(tpm_symbol))
)
write.csv(input_qc, qc_file, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(project_dir, "CIBERSORTx_input_sessionInfo.txt"))

cat("\nSamples:", ncol(tpm_symbol), "\n")
cat("Genes in count matrix:", nrow(counts), "\n")
cat("Genes with GTF-derived length:", sum(length_found), "\n")
cat("Mapped genes before symbol collapsing:", sum(symbol_found), "\n")
cat("Final unique symbols:", nrow(tpm_symbol), "\n")
cat("Mixture column sums after symbol filtering:\n")
print(summary(colSums(tpm_symbol)))
cat("\nSaved:", mixture_file, "\n")
cat("CIBERSORTx input preparation completed.\n")

