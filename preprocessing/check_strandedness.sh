#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$HOME/RNAseq"
ALIGN_DIR="$PROJECT_DIR/alignment"
OUT_FILE="$ALIGN_DIR/strandedness_check.tsv"

printf "Sample\tUnstranded\tForward\tReverse\tForward_fraction\tReverse_fraction\n" \
    > "$OUT_FILE"

shopt -s nullglob
files=("$ALIGN_DIR"/*_ReadsPerGene.out.tab)

if (( ${#files[@]} == 0 )); then
    echo "ERROR: No STAR ReadsPerGene files found." >&2
    exit 1
fi

for file in "${files[@]}"; do
    sample=$(basename "$file" _ReadsPerGene.out.tab)

    awk -v sample="$sample" '
        $1 !~ /^N_/ {
            unstranded += $2
            forward += $3
            reverse += $4
        }
        END {
            total_stranded = forward + reverse

            if (total_stranded > 0) {
                forward_fraction = forward / total_stranded
                reverse_fraction = reverse / total_stranded
            } else {
                forward_fraction = 0
                reverse_fraction = 0
            }

            printf "%s\t%.0f\t%.0f\t%.0f\t%.4f\t%.4f\n",
                   sample,
                   unstranded,
                   forward,
                   reverse,
                   forward_fraction,
                   reverse_fraction
        }
    ' "$file" >> "$OUT_FILE"
done

echo "Samples checked: ${#files[@]}"
echo "Output: $OUT_FILE"

echo
echo "First 10 samples:"
column -t -s $'\t' "$OUT_FILE" | head -11

echo
awk -F'\t' '
    NR > 1 {
        forward_total += $3
        reverse_total += $4
    }
    END {
        if (reverse_total > forward_total * 2) {
            print "Overall interpretation: REVERSE-STRANDED"
            print "Use STAR ReadsPerGene column 4."
        } else if (forward_total > reverse_total * 2) {
            print "Overall interpretation: FORWARD-STRANDED"
            print "Use STAR ReadsPerGene column 3."
        } else {
            print "Overall interpretation: UNSTRANDED OR AMBIGUOUS"
            print "Use STAR ReadsPerGene column 2 only after confirming library preparation."
        }
    }
' "$OUT_FILE"
