# Project 2B — Reactome pathway over-representation analysis
# Evidence-based reconstruction from preserved output/object names and plots.
# No dedicated Reactome R history survived, so compare regenerated outputs
# against the preserved historical CSV and RDS files before replacing them.

source("scripts/00_project_setup.R")
require_packages(c(
  "ReactomePA", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2"
))

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

tested_protein_coding_genes <- res_annotated[
  !is.na(res_annotated$padj) &
    res_annotated$gene_biotype == "protein_coding" &
    !is.na(res_annotated$gene_name) & res_annotated$gene_name != "",
]

ER_pos_entrez <- clusterProfiler::bitr(
  unique(ER_pos_up_genes$gene_name),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db::org.Hs.eg.db
)
TNBC_entrez <- clusterProfiler::bitr(
  unique(TNBC_up_genes$gene_name),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db::org.Hs.eg.db
)
universe_entrez <- clusterProfiler::bitr(
  unique(tested_protein_coding_genes$gene_name),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db::org.Hs.eg.db
)

reactome_ER_pos_up <- ReactomePA::enrichPathway(
  gene = unique(ER_pos_entrez$ENTREZID),
  universe = unique(universe_entrez$ENTREZID),
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  readable = TRUE
)
reactome_TNBC_up <- ReactomePA::enrichPathway(
  gene = unique(TNBC_entrez$ENTREZID),
  universe = unique(universe_entrez$ENTREZID),
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  readable = TRUE
)

write_csv_safe(
  as.data.frame(reactome_ER_pos_up),
  "results/reactome_outputs/Reactome_enrichment_ER_pos_up_genes.csv"
)
write_csv_safe(
  as.data.frame(reactome_TNBC_up),
  "results/reactome_outputs/Reactome_enrichment_TNBC_up_genes.csv"
)
save_rds_safe(
  reactome_ER_pos_up,
  "results/r_objects/reactome_ER_pos_up.rds"
)
save_rds_safe(
  reactome_TNBC_up,
  "results/r_objects/reactome_TNBC_up.rds"
)

plot_Reactome_ER_pos <- enrichplot::dotplot(
  reactome_ER_pos_up,
  showCategory = 15
) + ggplot2::ggtitle("Reactome enrichment: ER_pos-up genes")
plot_Reactome_TNBC <- enrichplot::dotplot(
  reactome_TNBC_up,
  showCategory = 15
) + ggplot2::ggtitle("Reactome enrichment: TNBC-up genes")

save_plot_safe(
  plot_Reactome_ER_pos,
  "figures/Reactome_dotplot_ER_pos_up_genes.png"
)
save_plot_safe(
  plot_Reactome_TNBC,
  "figures/Reactome_dotplot_TNBC_up_genes.png"
)
save_rds_safe(
  plot_Reactome_ER_pos,
  "results/r_objects/plot_Reactome_dotplot_ER_pos_up_genes.rds"
)
save_rds_safe(
  plot_Reactome_TNBC,
  "results/r_objects/plot_Reactome_dotplot_TNBC_up_genes.rds"
)

write_lines_safe(
  capture.output(sessionInfo()),
  "logs/sessionInfo_Project2B_Reactome_reconstructed.txt"
)

message("Reactome enrichment reconstruction completed.")
message("Compare row counts and top pathways with the preserved historical outputs.")
