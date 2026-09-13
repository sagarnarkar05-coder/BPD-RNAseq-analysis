files <- list.files(
    path = ".",
    pattern = "_ReadsPerGene\\.out\\.tab$",
    full.names = TRUE
)

if (length(files) == 0) {
    stop("No STAR ReadsPerGene.out.tab files were found.")
}

cat("Count files found:", length(files), "\n")

if (length(files) != 66) {
    warning("Expected 66 files, but found ", length(files))
}

# STAR column 4 = reverse-stranded counts
count_column <- 4

count_list <- lapply(files, function(file) {

    x <- read.delim(
        file,
        header = FALSE,
        stringsAsFactors = FALSE
    )

    if (ncol(x) < 4) {
        stop("File has fewer than four columns: ", file)
    }

    # Remove STAR summary rows such as N_unmapped and N_multimapping
    x <- x[!grepl("^N_", x[[1]]), c(1, count_column)]

    sample <- sub(
        "_ReadsPerGene\\.out\\.tab$",
        "",
        basename(file)
    )

    colnames(x) <- c("GeneID", sample)

    x
})

gene_ids <- count_list[[1]]$GeneID

for (i in seq_along(count_list)) {
    if (!identical(count_list[[i]]$GeneID, gene_ids)) {
        stop(
            "Gene order is inconsistent in file: ",
            basename(files[i])
        )
    }
}

count_matrix <- data.frame(
    GeneID = gene_ids,
    check.names = FALSE
)

for (x in count_list) {
    count_matrix[[colnames(x)[2]]] <- x[[2]]
}

if (anyDuplicated(count_matrix$GeneID)) {
    stop("Duplicated Gene IDs detected.")
}

if (any(is.na(count_matrix))) {
    stop("Missing values detected in the count matrix.")
}

write.table(
    count_matrix,
    file = "STAR_reverse_stranded_count_matrix.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "Created STAR_reverse_stranded_count_matrix.tsv\n",
    "Genes:", nrow(count_matrix), "\n",
    "Samples:", ncol(count_matrix) - 1, "\n"
)
