#!/bin/bash

samtools quickcheck -v \
*_Aligned.sortedByCoord.out.bam \
> bam_quickcheck_errors.txt

