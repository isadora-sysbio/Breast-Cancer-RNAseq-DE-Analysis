# Project 2B — VST, PCA, heatmap, and volcano plot
# Evidence-based reconstruction from preserved object and figure names.

source("scripts/00_project_setup.R")
require_packages(c("DESeq2", "ggplot2", "pheatmap"))

dds_path <- "results/r_objects/dds_after_DESeq2.rds"
res_path <- "results/r_objects/res_ER_pos_vs_TNBC.rds"
assert_file(dds_path)
assert_file(res_path)

dds <- readRDS(dds_path)
res <- readRDS(res_path)

vsd <- DESeq2::vst(dds, blind = FALSE)
save_rds_safe(vsd, "results/r_objects/vsd_for_PCA.rds")

# PCA
pca_data <- DESeq2::plotPCA(vsd, intgroup = "Group", returnData = TRUE)
percent_variance <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot2::ggplot(
  pca_data,
  ggplot2::aes(x = PC1, y = PC2, color = Group, shape = Group)
) +
  ggplot2::geom_point(size = 3.5) +
  ggplot2::xlab(paste0("PC1: ", percent_variance[[1]], "% variance")) +
  ggplot2::ylab(paste0("PC2: ", percent_variance[[2]], "% variance")) +
  ggplot2::ggtitle("PCA: ER_pos vs TNBC") +
  ggplot2::theme_bw()

save_plot_safe(
  pca_plot,
  "figures/PCA_ER_pos_vs_TNBC.png",
  width = 8,
  height = 6
)
save_rds_safe(pca_plot, "results/r_objects/pca_plot.rds")

# Top-50 heatmap: prioritize adjusted P value, then effect-size magnitude.
res_table <- as.data.frame(res)
res_table$gene_id <- rownames(res_table)
eligible <- res_table[!is.na(res_table$padj),]
eligible <- eligible[order(eligible$padj, -abs(eligible$log2FoldChange)),]
top50_gene_ids <- head(eligible$gene_id, 50L)

vst_matrix <- SummarizedExperiment::assay(vsd)[top50_gene_ids, , drop = FALSE]
row_z_scores <- t(scale(t(vst_matrix)))
row_z_scores[is.na(row_z_scores)] <- 0

column_annotation <- data.frame(
  Group = SummarizedExperiment::colData(vsd)$Group,
  row.names = colnames(vsd)
)

heatmap_path <- "figures/Heatmap_top50_DE_genes_ER_pos_vs_TNBC.png"
prepare_output(heatmap_path)
heatmap_object <- pheatmap::pheatmap(
  row_z_scores,
  annotation_col = column_annotation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Top 50 differential genes: ER_pos vs TNBC",
  filename = heatmap_path,
  width = 11,
  height = 10
)

# Volcano plot
volcano_data <- res_table
volcano_data$direction <- "Not_significant"
volcano_data$direction[
  !is.na(volcano_data$padj) & volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange > 0
] <- "Higher_in_ER_pos"
volcano_data$direction[
  !is.na(volcano_data$padj) & volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange < 0
] <- "Higher_in_TNBC"
volcano_data$minus_log10_padj <- -log10(volcano_data$padj)
volcano_data$minus_log10_padj[!is.finite(volcano_data$minus_log10_padj)] <- NA_real_

volcano_plot <- ggplot2::ggplot(
  volcano_data,
  ggplot2::aes(
    x = log2FoldChange,
    y = minus_log10_padj,
    color = direction
  )
) +
  ggplot2::geom_point(alpha = 0.55, size = 1.2, na.rm = TRUE) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggplot2::labs(
    title = "Volcano plot: ER_pos vs TNBC",
    x = "log2 fold change (positive = higher in ER_pos)",
    y = "-log10 adjusted P value",
    color = "Direction"
  ) +
  ggplot2::theme_bw()

save_plot_safe(
  volcano_plot,
  "figures/Volcano_ER_pos_vs_TNBC.png",
  width = 9,
  height = 7
)
save_rds_safe(volcano_plot, "results/r_objects/volcano_plot.rds")

checkpoint_path <- "results/r_objects/project2B_DESeq2_heatmap_checkpoint.RData"
prepare_output(checkpoint_path)
save(
  dds,
  vsd,
  res,
  pca_plot,
  volcano_plot,
  top50_gene_ids,
  row_z_scores,
  heatmap_object,
  file = checkpoint_path
)

message("PCA, top-50 heatmap, and volcano plot completed.")
