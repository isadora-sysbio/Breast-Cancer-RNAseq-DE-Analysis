# Project 2B — gene annotation, FDR filtering, and pathway input exports
# Evidence-based reconstruction from preserved filenames and reported counts.

source("scripts/00_project_setup.R")
require_packages(c("rtracklayer"))

res_path <- "results/r_objects/res_ER_pos_vs_TNBC.rds"
gtf_path <- "reference/Homo_sapiens.GRCh38.112.gtf"
assert_file(res_path)
assert_file(gtf_path)

res <- readRDS(res_path)
res_annotated <- as.data.frame(res)
res_annotated$gene_id <- rownames(res_annotated)

gtf <- rtracklayer::import(gtf_path)
gene_rows <- gtf[S4Vectors::mcols(gtf)$type == "gene"]
gene_metadata <- S4Vectors::mcols(gene_rows)

biotype_column <- if ("gene_biotype" %in% names(gene_metadata)) {
  "gene_biotype"
} else if ("gene_type" %in% names(gene_metadata)) {
  "gene_type"
} else {
  stop("The GTF has neither a gene_biotype nor gene_type column.")
}

gene_annotation <- unique(data.frame(
  gene_id = as.character(gene_metadata$gene_id),
  gene_name = as.character(gene_metadata$gene_name),
  gene_biotype = as.character(gene_metadata[[biotype_column]]),
  stringsAsFactors = FALSE
))
gene_annotation <- gene_annotation[!duplicated(gene_annotation$gene_id),]

annotation_index <- match(res_annotated$gene_id, gene_annotation$gene_id)
res_annotated$gene_name <- gene_annotation$gene_name[annotation_index]
res_annotated$gene_biotype <- gene_annotation$gene_biotype[annotation_index]

res_annotated$direction <- "Not_significant"
res_annotated$direction[
  !is.na(res_annotated$padj) & res_annotated$padj < 0.05 &
    res_annotated$log2FoldChange > 0
] <- "Higher_in_ER_pos"
res_annotated$direction[
  !is.na(res_annotated$padj) & res_annotated$padj < 0.05 &
    res_annotated$log2FoldChange < 0
] <- "Higher_in_TNBC"

preferred_columns <- c(
  "gene_id", "gene_name", "gene_biotype", "baseMean",
  "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "direction"
)
res_annotated <- res_annotated[, preferred_columns]
res_annotated <- res_annotated[order(res_annotated$padj, na.last = TRUE),]

save_rds_safe(
  res_annotated,
  "results/r_objects/res_annotated_ER_pos_vs_TNBC.rds"
)
write_csv_safe(
  res_annotated,
  "results/deseq2/ER_pos_vs_TNBC_DESeq2_results_annotated.csv"
)

significant <- res_annotated[
  !is.na(res_annotated$padj) & res_annotated$padj < 0.05,
]
significant <- significant[order(significant$padj, -abs(significant$log2FoldChange)),]

protein_coding <- significant[
  significant$gene_biotype == "protein_coding" &
    !is.na(significant$gene_name) & significant$gene_name != "",
]

top_ranked <- significant[order(significant$padj, -abs(significant$log2FoldChange)),]

save_rds_safe(
  significant,
  "results/r_objects/significant_DE_genes_FDR005_ER_pos_vs_TNBC.rds"
)
save_rds_safe(
  protein_coding,
  "results/r_objects/significant_protein_coding_DE_genes_FDR005_ER_pos_vs_TNBC.rds"
)
save_rds_safe(
  top_ranked,
  "results/r_objects/top_DE_genes_ranked_ER_pos_vs_TNBC.rds"
)

write_csv_safe(
  significant,
  "results/significant_DE_genes_FDR005_ER_pos_vs_TNBC.csv"
)
write_csv_safe(
  protein_coding,
  "results/significant_protein_coding_DE_genes_FDR005_ER_pos_vs_TNBC.csv"
)
write_csv_safe(
  top_ranked,
  "results/top_DE_genes_ranked_ER_pos_vs_TNBC.csv"
)

ER_pos_up <- protein_coding[protein_coding$direction == "Higher_in_ER_pos",]
TNBC_up <- protein_coding[protein_coding$direction == "Higher_in_TNBC",]

write_csv_safe(
  ER_pos_up,
  "results/pathway_inputs/ER_pos_up_DE_genes_FDR005_protein_coding.csv"
)
write_csv_safe(
  TNBC_up,
  "results/pathway_inputs/TNBC_up_DE_genes_FDR005_protein_coding.csv"
)
write_lines_safe(
  unique(ER_pos_up$gene_name),
  "results/pathway_inputs/ER_pos_up_gene_symbols_FDR005_protein_coding.txt"
)
write_lines_safe(
  unique(TNBC_up$gene_name),
  "results/pathway_inputs/TNBC_up_gene_symbols_FDR005_protein_coding.txt"
)

# Reconstructed marker panel: verify against the already-preserved historical CSV.
marker_symbols <- c(
  "ESR1", "PGR", "GATA3", "FOXA1", "TFF3", "CA12", "NAT1",
  "FOXC1", "KRT5", "KRT14", "KRT17", "EGFR", "VGLL1", "ELF5"
)
known_marker_check <- res_annotated[
  !is.na(res_annotated$gene_name) & res_annotated$gene_name %in% marker_symbols,
]
known_marker_check <- known_marker_check[
  match(marker_symbols, known_marker_check$gene_name, nomatch = 0L),
]
write_csv_safe(
  known_marker_check,
  "results/deseq2/known_marker_gene_check.csv"
)

observed_counts <- c(
  significant_total = nrow(significant),
  significant_ER_pos = sum(significant$direction == "Higher_in_ER_pos"),
  significant_TNBC = sum(significant$direction == "Higher_in_TNBC"),
  protein_coding_total = nrow(protein_coding),
  protein_coding_ER_pos = nrow(ER_pos_up),
  protein_coding_TNBC = nrow(TNBC_up)
)
expected_counts <- c(
  significant_total = 3674,
  significant_ER_pos = 1778,
  significant_TNBC = 1896,
  protein_coding_total = 3280,
  protein_coding_ER_pos = 1515,
  protein_coding_TNBC = 1765
)

count_check <- data.frame(
  category = names(expected_counts),
  expected_historical = unname(expected_counts),
  observed_reconstructed = unname(observed_counts),
  matches = unname(expected_counts == observed_counts)
)
write_csv_safe(count_check, "results/deseq2/reconstruction_count_check.csv")

if (!all(count_check$matches)) {
  warning(
    "One or more reconstructed DE counts differ from the preserved historical counts. ",
    "Do not overwrite the historical results; inspect reconstruction_count_check.csv."
  )
}

print(count_check)
message("Annotation and DE export stage completed.")
