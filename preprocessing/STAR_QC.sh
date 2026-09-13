#!/bin/bash

printf "Sample\tInput_reads\tUnique_pct\tMulti_pct\n" \
> STAR_alignment_QC.tsv

for log in *_Log.final.out
do

sample=$(basename "$log" _Log.final.out)

input=$(awk -F'|' '/Number of input reads/ {gsub(/[[:space:]]/,"",$2);print $2}' "$log")

unique=$(awk -F'|' '/Uniquely mapped reads %/ {gsub(/[[:space:]%]/,"",$2);print $2}' "$log")

multi=$(awk -F'|' '/% of reads mapped to multiple loci/ {gsub(/[[:space:]%]/,"",$2);print $2}' "$log")

printf "%s\t%s\t%s\t%s\n" \
"$sample" "$input" "$unique" "$multi" \
>> STAR_alignment_QC.tsv

done

