#!/usr/bin/env Rscript
# Sample-level QC for the primary RNA-seq analysis: sequencing/count metrics,
# sample-sample correlation, hierarchical clustering, and PCA.
#
# Works off the clean count matrix produced by count_matrix.R. Counts are
# variance-stabilised (DESeq2 vst, or a log2(CPM+1) fallback when there are too
# few genes for vst to fit) before correlation/PCA so that library size and the
# mean-variance trend do not drive the structure.
#
# Usage: Rscript sample_qc.R <gene_counts_csv> <metadata_csv> <out_dir>

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(ggrepel)
  library(glue)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript sample_qc.R <gene_counts_csv> <metadata_csv> <out_dir>")
}
counts_csv <- args[1]
metadata_path <- args[2]
out_dir <- args[3]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

counts_df <- read.csv(counts_csv, check.names = FALSE, stringsAsFactors = FALSE)
counts <- as.matrix(counts_df[, -1, drop = FALSE])
rownames(counts) <- counts_df[[1]]
mode(counts) <- "numeric"

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)

# Count matrix columns carry sample_id where available, otherwise the FASTQ name
key_col <- if ("sample_id" %in% colnames(metadata) &&
                all(colnames(counts) %in% metadata$sample_id)) "sample_id" else "sample"
coldata <- metadata[match(colnames(counts), metadata[[key_col]]), , drop = FALSE]
rownames(coldata) <- colnames(counts)

if (!"condition" %in% colnames(coldata)) coldata$condition <- NA_character_
coldata$condition <- ifelse(is.na(coldata$condition) | !nzchar(coldata$condition),
                            "unspecified", coldata$condition)
coldata$condition <- factor(coldata$condition, levels = unique(coldata$condition))

## ---- Sequencing / count metrics ---------------------------------------------
metrics <- data.frame(
  sample = colnames(counts),
  condition = as.character(coldata$condition),
  library_size = colSums(counts),
  genes_detected = colSums(counts > 0),
  genes_detected_min10 = colSums(counts >= 10),
  median_nonzero_count = apply(counts, 2, function(x) median(x[x > 0])),
  check.names = FALSE
)
write.csv(metrics, file.path(out_dir, "sample_metrics.csv"), row.names = FALSE)

## ---- Variance-stabilised expression -----------------------------------------
keep <- rowSums(counts >= 10) >= max(2, floor(0.1 * ncol(counts)))
counts_filt <- counts[keep, , drop = FALSE]
message(glue("{nrow(counts_filt)} of {nrow(counts)} genes pass the expression filter"))

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_filt),
  colData = coldata,
  design = ~1
)

transform_label <- "variance-stabilising transformation (DESeq2 vst)"
mat <- tryCatch(
  assay(vst(dds, blind = TRUE)),
  error = function(e) {
    warning(glue("vst failed ({e$message}); falling back to log2(CPM + 1)"))
    transform_label <<- "log2(CPM + 1)"
    cpm <- sweep(counts_filt, 2, colSums(counts_filt), "/") * 1e6
    log2(cpm + 1)
  }
)
writeLines(transform_label, file.path(out_dir, "transformation.txt"))

## ---- Sample-sample correlation ----------------------------------------------
cor_spearman <- cor(mat, method = "spearman")
cor_pearson <- cor(mat, method = "pearson")
write.csv(cor_spearman, file.path(out_dir, "sample_correlation_spearman.csv"))
write.csv(cor_pearson, file.path(out_dir, "sample_correlation_pearson.csv"))

annotation_col <- data.frame(condition = coldata$condition)
rownames(annotation_col) <- colnames(mat)

draw_correlation_heatmap <- function(cor_mat, title) {
  pheatmap(
    cor_mat,
    annotation_col = annotation_col,
    annotation_row = annotation_col,
    display_numbers = ncol(cor_mat) <= 16,
    number_format = "%.2f",
    fontsize_number = 6,
    fontsize_row = 8,
    fontsize_col = 8,
    main = title,
    silent = TRUE
  )
}

save_plot <- function(plot_obj, base_name, width = 10, height = 8) {
  pdf(file.path(out_dir, glue("{base_name}.pdf")), width = width, height = height)
  print(plot_obj)
  dev.off()
  png(file.path(out_dir, glue("{base_name}.png")), width = width, height = height,
      units = "in", res = 200)
  print(plot_obj)
  dev.off()
}

save_plot(
  draw_correlation_heatmap(cor_spearman, "Sample-sample Spearman correlation"),
  "sample_correlation_spearman_heatmap"
)
save_plot(
  draw_correlation_heatmap(cor_pearson, "Sample-sample Pearson correlation"),
  "sample_correlation_pearson_heatmap"
)

## ---- Hierarchical clustering ------------------------------------------------
sample_dist <- dist(t(mat))
write.csv(as.matrix(sample_dist), file.path(out_dir, "sample_distance_euclidean.csv"))

hc <- hclust(sample_dist, method = "complete")
saveRDS(hc, file.path(out_dir, "sample_clustering.rds"))

plot_dendrogram <- function() {
  plot(hc,
       main = "Sample clustering (Euclidean distance, complete linkage)",
       xlab = "", sub = "", cex = 0.8)
}
pdf(file.path(out_dir, "sample_clustering_dendrogram.pdf"), width = 10, height = 7)
plot_dendrogram()
dev.off()
png(file.path(out_dir, "sample_clustering_dendrogram.png"), width = 10, height = 7,
    units = "in", res = 200)
plot_dendrogram()
dev.off()

save_plot(
  pheatmap(
    as.matrix(sample_dist),
    clustering_distance_rows = sample_dist,
    clustering_distance_cols = sample_dist,
    annotation_col = annotation_col,
    fontsize_row = 8,
    fontsize_col = 8,
    main = "Sample-sample Euclidean distance",
    silent = TRUE
  ),
  "sample_distance_heatmap"
)

## ---- PCA ---------------------------------------------------------------------
# Restrict to the most variable genes so PCA reflects biological structure
# rather than the long tail of near-constant genes.
n_top <- min(500, nrow(mat))
top_var <- head(order(MatrixGenerics::rowVars(mat), decreasing = TRUE), n_top)
pca <- prcomp(t(mat[top_var, , drop = FALSE]), center = TRUE, scale. = FALSE)

percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- data.frame(
  sample = rownames(pca$x),
  # Keep the metadata ordering of the groups so the legend reads in
  # experimental order rather than alphabetically
  condition = coldata$condition,
  pca$x[, seq_len(min(5, ncol(pca$x))), drop = FALSE],
  check.names = FALSE
)
write.csv(pca_df, file.path(out_dir, "pca_coordinates.csv"), row.names = FALSE)
write.csv(
  data.frame(
    component = paste0("PC", seq_along(percent_var)),
    percent_variance = percent_var
  ),
  file.path(out_dir, "pca_variance_explained.csv"),
  row.names = FALSE
)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sample), size = 2.5, max.overlaps = 30, show.legend = FALSE) +
  labs(
    title = glue("PCA - top {n_top} variable genes"),
    subtitle = transform_label,
    x = glue("PC1 ({percent_var[1]}% variance)"),
    y = glue("PC2 ({percent_var[2]}% variance)"),
    colour = "Condition"
  ) +
  theme_bw()
save_plot(pca_plot, "pca_plot", width = 10, height = 7)

scree_df <- data.frame(
  component = factor(paste0("PC", seq_along(percent_var)),
                     levels = paste0("PC", seq_along(percent_var))),
  percent_variance = percent_var
)
scree_plot <- ggplot(head(scree_df, 10), aes(x = component, y = percent_variance)) +
  geom_col(fill = "#1f4788") +
  labs(title = "PCA scree plot", x = NULL, y = "Variance explained (%)") +
  theme_bw()
save_plot(scree_plot, "pca_scree_plot", width = 8, height = 5)

message("Sample QC complete.")
