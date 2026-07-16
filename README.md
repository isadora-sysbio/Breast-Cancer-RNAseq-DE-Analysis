# Breast Cancer RNA-seq Differential Expression Analysis

## Overview

This project explores transcriptomic differences between breast cancer molecular subtypes using publicly available RNA-seq data.

The objective is to develop a reproducible bioinformatics workflow for identifying differentially expressed genes and interpreting their potential biological relevance in breast cancer.

This project represents an application of computational genomics approaches to translational cancer research.

---

## Research Question

Which genes and biological pathways show differential expression between:

- Estrogen receptor-positive (ER+) breast cancer
- Triple-negative breast cancer (TNBC)

Understanding molecular differences between breast cancer subtypes may provide insights into disease biology and potential biomarkers.

---

## Biological Background

Breast cancer is a heterogeneous disease composed of molecularly distinct subtypes with differences in:

- receptor status
- gene expression profiles
- therapeutic response
- prognosis

ER-positive tumors frequently depend on estrogen receptor signaling, while TNBC represents a clinically aggressive subtype lacking ER, PR, and HER2 expression.

---

## Dataset

Publicly available human breast cancer RNA-seq data obtained from GEO.

Study:

- GSE58135

Cohort design:

- 10 ER-positive breast cancer primary tumor samples
- 10 triple-negative breast cancer primary tumor samples

Biological replicates were selected to enable differential expression analysis.

---

## Workflow

### 1. Cohort Selection

- Identification of relevant GEO samples
- Biological inclusion/exclusion criteria
- Metadata organization

### 2. RNA-seq Processing

Raw sequencing data processing:

FASTQ files

↓

Quality control

↓

Transcript quantification with Salmon

↓

Gene-level expression matrix

### 3. Differential Expression Analysis

Statistical analysis using:

- R
- DESeq2
- tximport

Outputs:

- Differentially expressed genes
- Log2 fold changes
- Adjusted p-values
- PCA plot
- Volcano plot
- Heatmap of top differentially expressed genes

### 4. Biological Interpretation

Functional interpretation using:

- Direction-specific gene lists
- GO Biological Process enrichment
- Redundancy-reduced GO term summaries

---

## Differential Expression Direction

The DESeq2 comparison was structured as:

ER_pos vs TNBC

Therefore:

- Positive log2FoldChange = higher expression in ER_pos tumors
- Negative log2FoldChange = higher expression in TNBC tumors

---

## Differential Expression Results

Using an FDR cutoff of padj < 0.05:

- Total significant DE genes: 3,674
- Higher in ER_pos: 1,778
- Higher in TNBC: 1,896

Significant protein-coding DE genes:

- Total: 3,280
- Higher in ER_pos: 1,515
- Higher in TNBC: 1,765

---

## GO Biological Process Enrichment

After differential expression analysis, significant protein-coding DE genes were separated by direction:

- ER_pos-up genes: 1,515
- TNBC-up genes: 1,765

GO Biological Process enrichment was performed separately for each direction-specific gene set using `clusterProfiler`, with the DESeq2-tested protein-coding genes used as the enrichment universe.

### Main findings

- ER_pos-up genes were enriched for biological processes related to hormone regulation, hormone transport/secretion, signal release, and cilium/microtubule-associated organization.
- TNBC-up genes were enriched for mitotic division, nuclear division, chromosome segregation, sister chromatid segregation, organelle fission, and spindle/checkpoint-associated processes.

GO redundancy reduction was performed using `clusterProfiler::simplify()` with a semantic similarity cutoff of 0.7.

These results suggest that ER_pos-up genes reflect hormone/luminal-associated transcriptional programs, while TNBC-up genes show a strong proliferative/cell-cycle-associated transcriptional signature.

These enrichment results indicate overrepresentation of GO terms among subtype-associated gene sets and do not prove causal pathway activation.

---

## Tools and Technologies

### Programming

- R
- Bash/Linux

### Bioinformatics

- Salmon
- tximport
- DESeq2
- clusterProfiler
- org.Hs.eg.db
- enrichplot
- ggplot2
- GEO/SRA resources
- Bioconductor ecosystem

### Data Management

- Git
- GitHub
- Reproducible project organization

---

## Project Goals

This project aims to build practical experience in:

- RNA-seq analysis
- computational genomics
- cancer biology
- transcriptomic interpretation
- pathway enrichment analysis
- reproducible research workflows

The long-term objective is to apply computational approaches to translational medicine and human genetics research.

---

## Project Status

Pipeline stages completed:

- Biological question definition
- GEO cohort selection
- Metadata organization
- Computational environment setup
- Salmon quantification
- tximport gene-level import
- DESeq2 differential expression analysis
- PCA visualization
- Volcano plot
- Top DE gene heatmap
- Direction-specific protein-coding gene lists
- GO Biological Process enrichment
- GO redundancy reduction

Next planned steps:

- Add Reactome/pathway-level enrichment
- Improve figure organization and captions
- Expand biological interpretation
- Prepare manuscript-style project summary
