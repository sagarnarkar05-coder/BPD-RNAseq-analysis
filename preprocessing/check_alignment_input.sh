#!/bin/bash

for log in ~/RNAseq/alignment/*_Log.out
do

sample=$(basename "$log" _Log.out)

grep readFilesIn "$log"

done

