# Breast Cancer RNA-seq Differential Expression Analysis

## Overview

This project explores transcriptomic differences between breast cancer molecular subtypes using publicly available RNA-seq data.

The objective is to develop a reproducible bioinformatics workflow for identifying differentially expressed genes and interpreting their potential biological relevance in breast cancer.

This project represents an application of computational genomics approaches to translational cancer research.

---

# Research Question

Which genes and biological pathways show differential expression between:

- Estrogen receptor-positive (ER+) breast cancer
- Triple-negative breast cancer (TNBC)

Understanding molecular differences between breast cancer subtypes may provide insights into disease biology and potential biomarkers.

---

# Biological Background

Breast cancer is a heterogeneous disease composed of molecularly distinct subtypes with differences in:

- receptor status
- gene expression profiles
- therapeutic response
- prognosis

ER-positive tumors frequently depend on estrogen receptor signaling, while TNBC represents a clinically aggressive subtype lacking ER, PR, and HER2 expression.

---

# Dataset

Publicly available human breast cancer RNA-seq data obtained from GEO.

Study:

GSE58135

Cohort design:

- 10 ER-positive breast cancer samples
- 10 Triple-negative breast cancer samples

Biological replicates were selected to enable differential expression analysis.

---

# Workflow

## 1. Cohort Selection

- Identification of relevant GEO samples
- Biological inclusion/exclusion criteria
- Metadata organization

## 2. RNA-seq Processing

Raw sequencing data processing:

FASTQ files

↓

Quality Control

↓

Transcript quantification

↓

Gene-level expression matrix

## 3. Differential Expression Analysis

Statistical analysis using:

- R
- DESeq2

Outputs:

- Differentially expressed genes
- Log2 fold changes
- Adjusted p-values
- Volcano plots
- Heatmaps

## 4. Biological Interpretation

Future analyses:

- Pathway enrichment
- Functional annotation
- Biological interpretation of significant genes

---

# Tools and Technologies

## Programming

- R
- Bash/Linux

## Bioinformatics

- Salmon
- DESeq2
- GEO/SRA resources
- Bioconductor ecosystem

## Data Management

- Git
- GitHub
- Reproducible project organization

---

# Project Goals

This project aims to build practical experience in:

- RNA-seq analysis
- computational genomics
- cancer biology
- transcriptomic interpretation
- reproducible research workflows

The long-term objective is to apply computational approaches to translational medicine and human genetics research.

---

# Project Status

🚧 Currently in development

Pipeline stages completed:

✅ Biological question definition  
✅ GEO cohort selection  
✅ Metadata organization  
✅ Computational environment setup  

Upcoming steps:

- RNA-seq quantification
- Differential expression analysis
- Biological interpretation
