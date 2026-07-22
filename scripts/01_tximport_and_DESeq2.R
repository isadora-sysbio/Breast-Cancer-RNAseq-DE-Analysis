# Project 2B — Salmon import and DESeq2 model
# Evidence-based reconstruction from saved session information, project notes,
# object names, and output filenames. Run after Salmon quantification is complete.

source("scripts/00_project_setup.R")
require_packages(c("tximport", "DESeq2", "rtracklayer"))

sample_design <- load_sample_design()
quant_files <- file.path("results", "salmon", sample_design$SRR, "quant.sf")
names(quant_files) <- rownames(sample_design)

missing_quant <- quant_files[!file.exists(quant_files)]
if (length(missing_quant) > 0) {
  stop(
    "Missing Salmon quant.sf file(s):\n",
    paste(missing_quant, collapse = "\n")
  )
}

# Build transcript-to-gene mapping from the exact Ensembl GRCh38.112 GTF.
gtf_path <- "reference/Homo_sapiens.GRCh38.112.gtf"
assert_file(gtf_path)
gtf <- rtracklayer::import(gtf_path)
transcript_rows <- gtf[S4Vectors::mcols(gtf)$type == "transcript"]

tx2gene <- unique(data.frame(
  TXNAME = as.character(S4Vectors::mcols(transcript_rows)$transcript_id),
  GENEID = as.character(S4Vectors::mcols(transcript_rows)$gene_id),
  stringsAsFactors = FALSE
))
tx2gene <- tx2gene[
  !is.na(tx2gene$TXNAME) & tx2gene$TXNAME != "" &
    !is.na(tx2gene$GENEID) & tx2gene$GENEID != "",
]

if (nrow(tx2gene) == 0L) {
  stop("The GTF import produced an empty transcript-to-gene mapping.")
}

# Preserve the known annotation compatibility issue as a reproducible QC table.
first_quant <- read.delim(quant_files[[1]], stringsAsFactors = FALSE)
if (!"Name" %in% names(first_quant)) {
  stop("The Salmon quant.sf file does not contain the expected Name column.")
}
missing_transcripts <- setdiff(unique(first_quant$Name), tx2gene$TXNAME)
write_csv_safe(
  data.frame(transcript_id = missing_transcripts),
  "results/qc/tximport_missing_transcripts.csv"
)
message(
  "Transcripts absent from tx2gene: ", length(missing_transcripts),
  " (historical note reported 13,448; compare after the run)."
)

# Import Salmon transcript estimates and summarize them to gene level.
txi <- tximport::tximport(
  files = quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "no"
)
save_rds_safe(txi, "results/r_objects/txi_salmon_gene_level.rds")

col_data <- sample_design[, "Group", drop = FALSE]
dds <- DESeq2::DESeqDataSetFromTximport(
  txi = txi,
  colData = col_data,
  design = ~ Group
)

# Preserve the historical workflow order: size factors first, then DESeq2.
dds <- DESeq2::estimateSizeFactors(dds)
size_factor_table <- data.frame(
  sample_id = names(DESeq2::sizeFactors(dds)),
  size_factor = unname(DESeq2::sizeFactors(dds))
)
write_csv_safe(size_factor_table, "results/deseq2/DESeq2_size_factors.csv")

dds <- DESeq2::DESeq(dds)
save_rds_safe(dds, "results/r_objects/dds_after_DESeq2.rds")

# Positive log2 fold change means higher expression in ER_pos than TNBC.
res <- DESeq2::results(
  dds,
  contrast = c("Group", "ER_pos", "TNBC"),
  alpha = 0.05
)
res <- res[order(res$padj, na.last = TRUE),]
save_rds_safe(res, "results/r_objects/res_ER_pos_vs_TNBC.rds")

res_table <- as.data.frame(res)
res_table$gene_id <- rownames(res_table)
res_table <- res_table[, c("gene_id", setdiff(names(res_table), "gene_id"))]
write_csv_safe(
  res_table,
  "results/deseq2/ER_pos_vs_TNBC_DESeq2_results.csv"
)

raw_counts <- as.data.frame(DESeq2::counts(dds, normalized = FALSE))
raw_counts$gene_id <- rownames(raw_counts)
raw_counts <- raw_counts[, c("gene_id", setdiff(names(raw_counts), "gene_id"))]
write_csv_safe(raw_counts, "counts/gene_counts_matrix.csv")

write_lines_safe(
  capture.output(sessionInfo()),
  "docs/R_sessionInfo_DESeq2_reconstructed.txt"
)

message("DESeq2 model completed.")
message("Contrast: ER_pos versus TNBC; positive log2FoldChange = higher in ER_pos.")
