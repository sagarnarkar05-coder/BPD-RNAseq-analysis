# BPD longitudinal RNA-seq analysis

This repository contains the scripts I used for the longitudinal RNA-seq reanalysis of preterm infants with and without bronchopulmonary dysplasia (BPD).

The RNA-seq data are publicly available from NCBI GEO under accession **GSE220135**. For this study I focused on the male infants who had samples available at all three time points: cord blood, Day 14 and Day 28.

The final dataset used here contains **66 samples from 22 infants**:

- 8 infants with BPD
- 14 infants without BPD
- 3 samples per infant

The main purpose of this reanalysis was to look at how the blood transcriptome changes over time in BPD and non-BPD preterm infants, rather than only comparing the two groups at a single time point.

## What is included

The workflow covers:

- FASTQ quality control
- STAR alignment and gene counting
- differential expression analysis
- FGSEA pathway enrichment
- GSVA pathway activity analysis
- CIBERSORTx immune-cell deconvolution
- WGCNA
- integration of cell fractions, GSVA scores and WGCNA modules
- scripts used to generate selected manuscript figures

## Repository layout

```text
BPD-RNAseq-analysis/
│
├── README.md
├── LICENSE
│
├── metadata/
│   └── Sample_metadata.txt
│
├── preprocessing/
│   ├── fastQC.sh
│   ├── multiQC.sh
│   ├── STAR_alignment.sh
│   ├── check_alignment_input.sh
│   ├── STAR_QC.sh
│   ├── BAM_integrity.sh
│   ├── check_strandedness.sh
│   ├── build_star_count_matrix.R
│   └── verify_count_matrix.sh
│
├── DEG/
│   ├── DEG.R
│   └── DEG_counts.R
│
├── FGSEA/
│   └── FGSEA.R
│
├── GSVA/
│   ├── GSVA_normalization.R
│   └── GSVA_run.R
│
├── CELL/
│   ├── CIBERSORTx_input.R
│   ├── CELL_analysis.R
│   └── CELL_trajectory_plots.R
│
├── WGCNA/
│   └── WGCNA.R
│
├── integration/
│   └── Cell_GSVA_WGCNA_Integration.R
│
└── reproducibility/
    └── software_versions.txt
```

## Order of analysis

The scripts were used in this order:

1. `fastQC.sh`
2. `multiQC.sh`
3. `STAR_alignment.sh`
4. `check_alignment_input.sh`
5. `STAR_QC.sh`
6. `BAM_integrity.sh`
7. `check_strandedness.sh`
8. `build_star_count_matrix.R`
9. `verify_count_matrix.sh`
10. `DEG.R`
11. `DEG_counts.R`
12. `FGSEA.R`
13. `GSVA_normalization.R`
14. `GSVA_run.R`
15. `CIBERSORTx_input.R`
16. CIBERSORTx Fractions
17. `CELL_analysis.R`
18. `CELL_trajectory_plots.R`
19. `WGCNA.R`
20. `Cell_GSVA_WGCNA_Integration.R`

## RNA-seq preprocessing

I used the following reference files and software for preprocessing:

- Genome assembly: **GRCh38**
- Annotation: **GENCODE v44**
- STAR: **2.7.11b**

FastQC was used to check the raw FASTQ files and MultiQC was used to summarize the QC reports.

STAR was used for alignment and gene counting with:

```text
--quantMode GeneCounts
```

The libraries were identified as **reverse-stranded**, so STAR column 4 from `ReadsPerGene.out.tab` was used to build the count matrix.

The final count matrix contained all 66 samples used in the downstream analyses.

## Differential expression

`DEG.R` was used for the gene-level analysis.

The following comparisons were included:

**BPD versus non-BPD**
- Overall
- Cord blood
- Day 14
- Day 28

**BPD progression**
- Day 14 vs Cord
- Day 28 vs Day 14
- Day 28 vs Cord

**non-BPD progression**
- Day 14 vs Cord
- Day 28 vs Day 14
- Day 28 vs Cord

For the longitudinal comparisons, repeated samples from the same infant were analysed as paired measurements.

Main settings:

- minimum count = 10
- gene retained if this count was reached in at least 8 samples
- BH-FDR < 0.05
- |log2 fold change| >= 1
- apeglm shrinkage for log2 fold-change estimates

The script also saves ranked gene statistics for FGSEA.

`DEG_counts.R` was used to summarize the number of significant upregulated and downregulated genes across the different comparisons and to generate the DEG count figure used in the manuscript.

## FGSEA

`FGSEA.R` uses the ranked statistics produced by the DEG analysis.

I tested the following pathway collections:

- Hallmark
- Reactome
- Gene Ontology Biological Process
- KEGG Legacy

Main settings:

- minimum pathway size = 15
- maximum pathway size = 500
- BH-FDR < 0.05
- `fgseaMultilevel`
- `nPermSimple = 10000`

The same BPD vs non-BPD and longitudinal comparisons were used here.

## GSVA

GSVA was run in two steps.

`GSVA_normalization.R` prepares the expression matrix by filtering low-count genes, applying DESeq2 variance-stabilizing transformation and converting Ensembl IDs to gene symbols.

`GSVA_run.R` then calculates sample-level Hallmark pathway scores and performs the statistical comparisons.

Main settings:

- Hallmark gene sets
- Gaussian GSVA kernel
- minimum gene-set size = 10
- maximum gene-set size = 500
- BH-FDR < 0.05
- limma for statistical testing
- infant ID used to account for repeated measurements where required

The normalized expression matrix generated here was also used as input for WGCNA.

## CIBERSORTx analysis

`CIBERSORTx_input.R` converts the count matrix to TPM and prepares the input file for CIBERSORTx.

Gene lengths were calculated from the GENCODE v44 annotation and Ensembl IDs were converted to gene symbols before the final mixture file was created.

CIBERSORTx Fractions was run using:

- relative mode
- LM22 signature matrix
- 1000 permutations
- B-mode batch correction
- quantile normalization disabled

The resulting CIBERSORTx output was then analysed with `CELL_analysis.R`.

The downstream cell analysis includes:

- CIBERSORTx sample QC
- cell-type detection checks
- exclusion of sparsely detected cell types from statistical testing
- overall BPD vs non-BPD analysis
- time point-specific comparisons
- paired longitudinal comparisons
- BH-FDR correction

All 66 samples passed the CIBERSORTx deconvolution QC threshold in the reported analysis.

`CELL_trajectory_plots.R` was used to generate the trajectory plots for neutrophils, naive CD4 T cells and naive B cells.

## WGCNA

`WGCNA.R` was used to identify co-expression modules across the longitudinal dataset.

Main settings:

- signed network
- biweight midcorrelation (`bicor`)
- scale-free topology target = 0.85
- minimum module size = 30
- merge cut height = 0.15
- maximum block size = 5000
- BH-FDR < 0.05

The script includes sample and gene QC, soft-threshold selection, module detection, module eigengene analysis and hub-gene ranking.

The module eigengenes were saved for the final integration analysis.

## Integration analysis

`Cell_GSVA_WGCNA_Integration.R` brings together:

- CIBERSORTx immune-cell fractions
- GSVA Hallmark pathway scores
- WGCNA module eigengenes

I tested three types of associations:

1. cell type vs GSVA pathway
2. cell type vs WGCNA module
3. WGCNA module vs GSVA pathway

The mixed-effects models include BPD group, time point and their interaction, with infant ID included as a random intercept.

BH-FDR < 0.05 was used for the final multiple-testing correction.

## Metadata

`Sample_metadata.txt` contains:

- sample accession
- infant ID
- time point
- BPD group
- gestational age
- birth weight

The sample IDs were checked against the relevant count/expression matrices before each analysis.

## Software

The main software used in the pipeline includes:

- FastQC 0.12.1
- MultiQC 1.18
- STAR 2.7.11b
- samtools 1.23
- R 4.6.1

Important R/Bioconductor packages include DESeq2, limma, edgeR, apeglm, fgsea, msigdbr, GSVA, WGCNA, nlme, lmerTest, emmeans, ggplot2 and pheatmap.

The exact preprocessing versions are recorded in:

```text
reproducibility/software_versions.txt
```

Several R scripts also save `sessionInfo()` files during the analysis.

## Notes for rerunning the code

The scripts in this repository are the scripts used for the reported analysis. Some of them still contain the original local paths, for example:

```text
/home/sagar007/RNAseq/
```

Anyone rerunning the analysis on another computer will need to replace these paths with their own project directory.

Large files such as FASTQ files, BAM files, STAR genome indexes and other large intermediate files are not included in this repository.

## Data availability

The RNA-seq data used in this reanalysis are publicly available from NCBI GEO:

**GSE220135**

## Code availability

The scripts used for the analyses reported in the manuscript are provided in this repository.

A versioned release of the repository will be archived in Zenodo for long-term access.

**Zenodo DOI:** TO BE ADDED

## Citation

The manuscript citation will be added after publication.

## License

See the `LICENSE` file.
