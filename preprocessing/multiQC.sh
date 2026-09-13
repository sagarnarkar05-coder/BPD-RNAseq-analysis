#!/bin/bash

mkdir -p ~/RNAseq/multiqc_fastqc

multiqc ~/RNAseq/fastqc \
-o ~/RNAseq/multiqc_fastqc \
--force

