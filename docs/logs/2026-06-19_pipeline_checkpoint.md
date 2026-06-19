# Project 2B checkpoint — 2026-06-19

## Project
Breast-Cancer-RNAseq-DE-Analysis

## Biological comparison
ER+ primary breast tumors vs TNBC primary breast tumors

## Samples tested today
- SRR1313090 = ER+ pilot sample
- SRR1313132 = TNBC pilot sample

## Completed today
- Fixed DNS / SRA tools issue.
- Confirmed rnaseq conda environment works.
- Downloaded and converted SRR1313090 and SRR1313132.
- Ran FastQC on both paired-end samples.
- Generated MultiQC report.
- Built and tested Salmon transcriptome quantification using GRCh38 Ensembl cDNA index.
- Successfully quantified:
  - results/salmon/SRR1313090/quant.sf
  - results/salmon/SRR1313132/quant.sf
- Built a storage-safe batch script:
  - detects existing FASTQ
  - detects .fastq.gz
  - runs FastQC
  - runs Salmon
  - compresses FASTQ after successful processing
  - removes temporary SRA folders
- Compressed pilot FASTQs from ~45G to ~12G.
- Confirmed project has ~345G free disk space after cleanup.

## Current status
2-sample pilot pipeline is successful.
Do not run full 20-sample batch until script is switched from srr_test_list.txt to clean full SRR list.
