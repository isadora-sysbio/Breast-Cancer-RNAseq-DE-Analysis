# Project 2B — GO Biological Process over-representation analysis
# Direct reconstruction from R_history_GO_BP_enrichment_2026-07-16.Rhistory.

source("scripts/00_project_setup.R")
require_packages(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2"))

sig_path <- "results/r_objects/significant_protein_coding_DE_genes_FDR005_ER_pos_vs_TNBC.rds"
annotated_path <- "results/r_objects/res_annotated_ER_pos_vs_TNBC.rds"
assert_file(sig_path)
assert_file(annotated_path)

sig_protein_coding_DE_genes <- readRDS(sig_path)
res_annotated <- readRDS(annotated_path)

ER_pos_up_genes <- sig_protein_coding_DE_genes[
  sig_protein_coding_DE_genes$direction == "Higher_in_ER_pos",
]
TNBC_up_genes <- sig_protein_coding_DE_genes[
  sig_protein_coding_DE_genes$direction == "Higher_in_TNBC",
]

ER_pos_up_symbols <- unique(ER_pos_up_genes$gene_name)
TNBC_up_symbols <- unique(TNBC_up_genes$gene_name)

ER_pos_entrez <- clusterProfiler::bitr(
  ER_pos_up_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db::org.Hs.eg.db
)
TNBC_entrez <- clusterProfiler::bitr(
  TNBC_up_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db::org.Hs.eg.db
)

tested_protein_coding_genes <- res_annotated[
  !is.na(res_annotated$padj) &
    res_annotated$gene_biotype == "protein_coding" &
    !is.na(res_annotated$gene_name) &
    res_annotated$gene_name != "",
]
universe_symbols <- unique(tested_protein_coding_genes$gene_name)
universe_entrez <- clusterProfiler::bitr(
  universe_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db::org.Hs.eg.db
)

ER_pos_entrez_ids <- unique(ER_pos_entrez$ENTREZID)
TNBC_entrez_ids <- unique(TNBC_entrez$ENTREZID)
universe_entrez_ids <- unique(universe_entrez$ENTREZID)

ego_ER_pos_BP <- clusterProfiler::enrichGO(
  gene = ER_pos_entrez_ids,
  universe = universe_entrez_ids,
  OrgDb = org.Hs.eg.db::org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_TNBC_BP <- clusterProfiler::enrichGO(
  gene = TNBC_entrez_ids,
  universe = universe_entrez_ids,
  OrgDb = org.Hs.eg.db::org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_ER_pos_BP_df <- as.data.frame(ego_ER_pos_BP)
ego_TNBC_BP_df <- as.data.frame(ego_TNBC_BP)

write_csv_safe(
  ego_ER_pos_BP_df,
  "results/pathway_outputs/GO_BP_enrichment_ER_pos_up_genes.csv"
)
write_csv_safe(
  ego_TNBC_BP_df,
  "results/pathway_outputs/GO_BP_enrichment_TNBC_up_genes.csv"
)
save_rds_safe(ego_ER_pos_BP, "results/r_objects/ego_ER_pos_BP.rds")
save_rds_safe(ego_TNBC_BP, "results/r_objects/ego_TNBC_BP.rds")

p_ER_pos_BP <- enrichplot::dotplot(ego_ER_pos_BP, showCategory = 15) +
  ggplot2::ggtitle("GO Biological Process enrichment: ER_pos-up genes")
p_TNBC_BP <- enrichplot::dotplot(ego_TNBC_BP, showCategory = 15) +
  ggplot2::ggtitle("GO Biological Process enrichment: TNBC-up genes")

save_plot_safe(
  p_ER_pos_BP,
  "figures/GO_BP_dotplot_ER_pos_up_genes.png"
)
save_plot_safe(
  p_TNBC_BP,
  "figures/GO_BP_dotplot_TNBC_up_genes.png"
)

# Historical combined top-10 comparison table.
top_ER_pos_GO_BP <- head(ego_ER_pos_BP_df, 10)
top_TNBC_GO_BP <- head(ego_TNBC_BP_df, 10)
top_ER_pos_GO_BP$Gene_set <- "ER_pos_up"
top_TNBC_GO_BP$Gene_set <- "TNBC_up"
combined_top_GO_BP <- rbind(top_ER_pos_GO_BP, top_TNBC_GO_BP)
combined_top_GO_BP <- combined_top_GO_BP[, c(
  "Gene_set", "ID", "Description", "GeneRatio", "BgRatio",
  "p.adjust", "qvalue", "Count", "geneID"
)]
write_csv_safe(
  combined_top_GO_BP,
  "results/pathway_outputs/top10_GO_BP_enrichment_ER_pos_up_vs_TNBC_up.csv"
)

summary_lines <- c(
  "# Project 2B pathway enrichment summary",
  "",
  "Comparison: ER_pos breast tumors vs TNBC breast tumors",
  "",
  "Differential-expression direction:",
  "- Positive log2FoldChange = higher expression in ER_pos tumors",
  "- Negative log2FoldChange = higher expression in TNBC tumors",
  "",
  paste0("Significant DE genes (padj < 0.05): ", nrow(sig_protein_coding_DE_genes)),
  paste0("- Protein-coding higher in ER_pos: ", nrow(ER_pos_up_genes)),
  paste0("- Protein-coding higher in TNBC: ", nrow(TNBC_up_genes)),
  "",
  "GO Biological Process enrichment was performed using a DESeq2-tested protein-coding gene universe.",
  "Benjamini-Hochberg adjustment was used.",
  "",
  "Caution: over-representation does not prove causal pathway activation."
)
write_lines_safe(
  summary_lines,
  "results/pathway_outputs/project2B_pathway_enrichment_summary.md"
)

write_lines_safe(
  capture.output(sessionInfo()),
  "logs/sessionInfo_Project2B_GO_BP_reconstructed.txt"
)

message("GO BP enrichment completed.")
