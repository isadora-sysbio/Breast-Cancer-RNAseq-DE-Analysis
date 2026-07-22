# Project 2B — GO BP redundancy reduction and simplified plots
# Direct reconstruction from the two preserved GO history files.

source("scripts/00_project_setup.R")
require_packages(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2"))

ER_path <- "results/r_objects/ego_ER_pos_BP.rds"
TNBC_path <- "results/r_objects/ego_TNBC_BP.rds"
assert_file(ER_path)
assert_file(TNBC_path)

ego_ER_pos_BP <- readRDS(ER_path)
ego_TNBC_BP <- readRDS(TNBC_path)

ego_ER_pos_BP_simplified <- clusterProfiler::simplify(
  ego_ER_pos_BP,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)
ego_TNBC_BP_simplified <- clusterProfiler::simplify(
  ego_TNBC_BP,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

ego_ER_pos_BP_simplified_df <- as.data.frame(ego_ER_pos_BP_simplified)
ego_TNBC_BP_simplified_df <- as.data.frame(ego_TNBC_BP_simplified)

write_csv_safe(
  ego_ER_pos_BP_simplified_df,
  "results/pathway_outputs/GO_BP_enrichment_ER_pos_up_genes_simplified.csv"
)
write_csv_safe(
  ego_TNBC_BP_simplified_df,
  "results/pathway_outputs/GO_BP_enrichment_TNBC_up_genes_simplified.csv"
)
save_rds_safe(
  ego_ER_pos_BP_simplified,
  "results/r_objects/ego_ER_pos_BP_simplified.rds"
)
save_rds_safe(
  ego_TNBC_BP_simplified,
  "results/r_objects/ego_TNBC_BP_simplified.rds"
)

p_ER_pos_BP_simplified <- enrichplot::dotplot(
  ego_ER_pos_BP_simplified,
  showCategory = 15
) + ggplot2::ggtitle("Simplified GO BP enrichment: ER_pos-up genes")

p_TNBC_BP_simplified <- enrichplot::dotplot(
  ego_TNBC_BP_simplified,
  showCategory = 15
) + ggplot2::ggtitle("Simplified GO BP enrichment: TNBC-up genes")

save_plot_safe(
  p_ER_pos_BP_simplified,
  "figures/GO_BP_dotplot_ER_pos_up_genes_simplified.png"
)
save_plot_safe(
  p_TNBC_BP_simplified,
  "figures/GO_BP_dotplot_TNBC_up_genes_simplified.png"
)
save_rds_safe(
  p_ER_pos_BP_simplified,
  "results/r_objects/plot_GO_BP_dotplot_ER_pos_up_genes_simplified.rds"
)
save_rds_safe(
  p_TNBC_BP_simplified,
  "results/r_objects/plot_GO_BP_dotplot_TNBC_up_genes_simplified.rds"
)

reduction_table <- data.frame(
  gene_set = c("ER_pos_up", "TNBC_up"),
  original_terms = c(
    nrow(as.data.frame(ego_ER_pos_BP)),
    nrow(as.data.frame(ego_TNBC_BP))
  ),
  simplified_terms = c(
    nrow(ego_ER_pos_BP_simplified_df),
    nrow(ego_TNBC_BP_simplified_df)
  )
)
write_csv_safe(
  reduction_table,
  "results/pathway_outputs/GO_BP_redundancy_reduction_counts.csv"
)

message("GO BP simplification completed with cutoff = 0.7.")
print(reduction_table)
