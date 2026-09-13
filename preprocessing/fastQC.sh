#!/bin/bash

mkdir -p ~/RNAseq/fastqc

fastqc ~/RNAseq/fastq/*.fastq.gz \
-o ~/RNAseq/fastqc \
-t 12

