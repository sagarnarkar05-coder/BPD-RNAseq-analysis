
#!/bin/bash

THREADS=12

for fq in ~/RNAseq/fastq/*_1.fastq.gz
do

sample=$(basename "$fq" _1.fastq.gz)

STAR \
--runThreadN $THREADS \
--genomeDir ~/RNAseq/reference \
--readFilesIn "$fq" \
--readFilesCommand zcat \
--outFileNamePrefix ~/RNAseq/alignment/${sample}_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts

done

