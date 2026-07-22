# Project 2B reconstructed analysis scripts

These scripts were reconstructed on 2026-07-22 from the surviving Project 2B evidence package.

## Confidence levels

- `run_salmon_batch.sh`: exact surviving script.
- `04_GO_BP_enrichment.R`: direct cleanup of the saved 2026-07-16 R history.
- `05_GO_BP_simplification_and_plots.R`: direct cleanup of the saved GO histories.
- `01_tximport_and_DESeq2.R`: evidence-based reconstruction from the DESeq2 session information, tximport annotation note, saved object names, output tables, and documented workflow.
- `02_DESeq2_QC_and_visualization.R`: evidence-based reconstruction from preserved figure and R-object names.
- `03_DE_gene_annotation_and_export.R`: evidence-based reconstruction from preserved filenames and historical result counts.
- `06_Reactome_enrichment.R`: inferred reconstruction. No dedicated Reactome command history survived, so historical outputs must be retained and compared.

These files do not replace the original `.Rhistory`, `.RData`, `.rds`, logs, figures, or CSV files. Keep those permanently as historical research evidence.

## Safety behavior

The scripts refuse to overwrite an existing output by default. This protects the historical Project 2B results.

Test them in a separate copy of the project. Only after the backup has been verified, overwrite protection can be disabled for a deliberate rerun:

```bash
export PROJECT2B_ALLOW_OVERWRITE=YES
```

Close the terminal or run `unset PROJECT2B_ALLOW_OVERWRITE` to restore protection.

## Installation into the project

From the Project 2B root:

```bash
mkdir -p scripts/reconstructed
cp /path/to/reconstructed/scripts/*.R scripts/reconstructed/
cp /path/to/reconstructed/scripts/run_salmon_batch.sh scripts/reconstructed/
```

The scripts call `scripts/00_project_setup.R`. Therefore, the recommended installation is to copy the `.R` files directly into the existing `scripts/` directory, after backing up that directory:

```bash
cp scripts scripts_before_reconstruction_2026-07-22 -a
cp /path/to/reconstructed/scripts/*.R scripts/
```

Do not replace the existing `run_salmon_batch.sh`; the packaged copy is included only as evidence.

## Run order

Activate the original environment and enter the project:

```bash
conda activate rnaseq
cd ~/Breast-Cancer-RNAseq-DE-Analysis
```

Check R syntax without running the analyses:

```bash
for script in scripts/0*.R; do
  Rscript -e "parse(file='$script')" || exit 1
done
```

Then, in a clean project test copy, run:

```bash
Rscript scripts/01_tximport_and_DESeq2.R
Rscript scripts/02_DESeq2_QC_and_visualization.R
Rscript scripts/03_DE_gene_annotation_and_export.R
Rscript scripts/04_GO_BP_enrichment.R
Rscript scripts/05_GO_BP_simplification_and_plots.R
Rscript scripts/06_Reactome_enrichment.R
```

`00_project_setup.R` is sourced automatically and is not run separately.

## Required preserved inputs

- `metadata/sample_sheet.csv`
- all 20 `results/salmon/<SRR>/quant.sf` directories
- `reference/Homo_sapiens.GRCh38.112.gtf`
- the original conda/R environment or compatible Bioconductor packages

## Validation targets

Compare reconstructed outputs with the preserved historical evidence rather than deleting or replacing it.

Historical targets recorded in the project:

- 20 samples: 10 ER_pos and 10 TNBC.
- Positive log2 fold change means higher expression in ER_pos.
- 3,674 genes with `padj < 0.05`.
- 1,778 higher in ER_pos and 1,896 higher in TNBC.
- 3,280 significant protein-coding genes.
- 1,515 protein-coding higher in ER_pos and 1,765 higher in TNBC.
- GO BP simplification reduced ER_pos terms from 67 to 31 and TNBC terms from 147 to 47.
- The tximport note recorded 13,448 unmatched transcript IDs, approximately 7%.

Small differences can occur if package versions, annotation contents, factor levels, filtering, or database releases differ. A mismatch is a signal to compare the reconstructed steps against the preserved `.rds` objects and CSV files—not a reason to overwrite the original evidence.
