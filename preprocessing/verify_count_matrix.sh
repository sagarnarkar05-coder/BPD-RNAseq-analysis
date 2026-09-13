#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$HOME/RNAseq"
COUNT_FILE="$PROJECT_DIR/counts/gene_counts_reverse.tsv"

if [[ ! -s "$COUNT_FILE" ]]; then
    echo "ERROR: Count matrix not found or empty:" >&2
    echo "$COUNT_FILE" >&2
    exit 1
fi

sample_count=$(
    head -1 "$COUNT_FILE" |
        tr '\t' '\n' |
        tail -n +2 |
        wc -l
)

duplicate_samples=$(
    head -1 "$COUNT_FILE" |
        tr '\t' '\n' |
        tail -n +2 |
        sort |
        uniq -d |
        wc -l
)

line_count=$(wc -l < "$COUNT_FILE")
gene_count=$((line_count - 1))

echo "Count file: $COUNT_FILE"
echo "Genes: $gene_count"
echo "Samples: $sample_count"
echo "Duplicated sample names: $duplicate_samples"

echo
echo "First two rows:"
head -2 "$COUNT_FILE" | cut -f1-6 | column -t -s $'\t'

if [[ "$sample_count" -eq 66 &&
      "$duplicate_samples" -eq 0 ]]; then
    echo
    echo "PASS: Count matrix structure is valid."
else
    echo
    echo "WARNING: Count matrix validation failed."
    exit 1
fi
