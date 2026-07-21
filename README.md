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

## Key Figures

### PCA plot

**File:** `figures/PCA_ER_pos_vs_TNBC.png`

The PCA plot summarizes overall sample-level expression structure after variance-stabilizing transformation. It shows whether ER_pos and TNBC samples separate based on broad transcriptomic patterns. PCA is used for exploratory visualization and quality assessment, not for identifying individual significant genes.

### Volcano plot

**File:** `figures/Volcano_ER_pos_vs_TNBC.png`

The volcano plot displays gene-level differential expression results from DESeq2. The x-axis represents log2 fold change, showing direction and magnitude of expression difference. The y-axis represents statistical significance. In this project, positive log2FoldChange indicates higher expression in ER_pos tumors, while negative log2FoldChange indicates higher expression in TNBC tumors.

### Top DE gene heatmap

**File:** `figures/Heatmap_top50_DE_genes_ER_pos_vs_TNBC.png`

The heatmap shows scaled expression values for the top 50 differentially expressed genes ranked by adjusted p-value. Rows represent genes and columns represent samples. The heatmap provides visual support that the strongest differentially expressed genes separate ER_pos and TNBC tumors into subtype-associated expression patterns.

### GO Biological Process dotplots

**Files:**

- `figures/GO_BP_dotplot_ER_pos_up_genes.png`
- `figures/GO_BP_dotplot_TNBC_up_genes.png`
- `figures/GO_BP_dotplot_ER_pos_up_genes_simplified.png`
- `figures/GO_BP_dotplot_TNBC_up_genes_simplified.png`

The GO dotplots summarize biological processes enriched among direction-specific gene sets. Dot position represents GeneRatio, dot size represents the number of genes in each GO term, and color represents adjusted p-value. Simplified dotplots reduce redundant GO terms and are used for clearer biological interpretation.

---

## Key Findings

1. ER_pos and TNBC tumors show distinct transcriptomic profiles based on PCA, differential expression analysis, heatmap visualization, and pathway enrichment.

2. DESeq2 identified 3,674 significant differentially expressed genes using an FDR cutoff of `padj < 0.05`.

3. Among significant protein-coding genes, 1,515 were higher in ER_pos tumors and 1,765 were higher in TNBC tumors.

4. ER_pos-up genes were enriched for hormone regulation, hormone transport/secretion, signal release, and cilium/microtubule-associated biological processes.

5. TNBC-up genes were enriched for mitotic division, nuclear division, chromosome segregation, sister chromatid segregation, organelle fission, and spindle/checkpoint-associated biological processes.

6. These results suggest that ER_pos tumors show hormone/luminal-associated transcriptional programs, while TNBC tumors show a stronger proliferative and cell-cycle-associated transcriptional signature.

7. These findings are based on association from public RNA-seq data and pathway overrepresentation analysis. They do not prove causal pathway activation or direct treatment response.

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

---

## Reproducibility and Resume Instructions

This project was run in a Conda environment named `rnaseq`.

To resume the project locally:

```bash
conda activate rnaseq
cd ~/Breast-Cancer-RNAseq-DE-Analysis
```

Main saved analysis objects are stored in:

```text
results/r_objects/
```

Main result tables are stored in:

```text
results/
results/pathway_inputs/
results/pathway_outputs/
```

Main figures are stored in:

```text
figures/
```

Key R objects include:

- `dds_after_DESeq2.rds`
- `res_annotated_ER_pos_vs_TNBC.rds`
- `significant_protein_coding_DE_genes_FDR005_ER_pos_vs_TNBC.rds`
- `ego_ER_pos_BP.rds`
- `ego_TNBC_BP.rds`
- `ego_ER_pos_BP_simplified.rds`
- `ego_TNBC_BP_simplified.rds`

The Conda environment files are included for reproducibility:

- `rnaseq_environment.yml`
- `rnaseq_conda_explicit.txt`

Large raw sequencing files and reference/index files are not intended to be tracked directly in GitHub.

---

## Final Project Status

Core analysis complete.

This repository now contains a reproducible public RNA-seq analysis comparing ER-positive and triple-negative breast cancer tumors using Salmon quantification, tximport, DESeq2 differential expression analysis, visualization, and GO Biological Process enrichment.

The project is considered portfolio-complete as an RNA-seq differential expression and pathway interpretation project.

### Completed outputs

- Curated cohort design: 10 ER_pos primary tumor samples and 10 TNBC primary tumor samples
- Salmon transcript quantification
- Gene-level import with tximport
- DESeq2 differential expression analysis
- PCA plot
- Volcano plot
- Top differentially expressed gene heatmap
- Significant DE gene tables
- Protein-coding filtered DE gene tables
- Direction-specific ER_pos-up and TNBC-up gene lists
- GO Biological Process enrichment
- Simplified GO Biological Process dotplots
- README documentation with biological interpretation and reproducibility notes

### Limitations

- This is a public-data reanalysis project and does not include experimental validation.
- The analysis identifies transcriptomic associations, not causal mechanisms.
- GO enrichment suggests overrepresented biological processes but does not prove pathway activation.
- Clinical treatment response cannot be inferred directly from these RNA-seq results.
- Future extensions could include Reactome enrichment, GSEA, validation in an independent cohort, or comparison with clinical metadata if available.

### Portfolio interpretation

This project demonstrates practical skills in RNA-seq analysis, differential expression, pathway enrichment, biological interpretation, Linux/R workflow management, and reproducible GitHub documentation.
