#!/usr/bin/env Rscript
# limma-voom differential expression analysis
# Usage: Rscript limma_voom_analysis.R counts.txt metadata.csv output_dir comparisons_config.yaml
#
# Parallel to deseq2_analysis.R: reads the same featureCounts matrix, the same
# metadata.csv and the same comparisons config, so the two methods can be
# compared contrast-for-contrast. Samples are grouped by NPC line-ID prefix into
# line_groups, each `comparisons` entry contrasts two line_groups, and each is
# run under every combination of `run_variants` (collapse_replicates x
# include_male_samples).
#
# Counts/metadata loading and the annotation helpers are deliberately kept
# self-contained (mirroring deseq2_analysis.R rather than sharing a module), so
# that editing one method's script cannot break the other's pipeline rule.

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(tidyverse)
  library(glue)
  library(AnnotationDbi)
  library(yaml)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)

default_counts <- "output/feature_count/all_samples_counts.txt"
default_meta <- "metadata/metadata.csv"
default_out <- "output/limma_voom"
default_comparisons_config <- "proj_src/deseq2_comparisons.yaml"

counts_file <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_counts
meta_file <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_meta
out_dir <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_out
comparisons_config_path <- if (length(args) >= 4 && nzchar(args[4])) args[4] else default_comparisons_config

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

## ---- Load counts -----------------------------------------------------------

counts <- read.table(counts_file, header = TRUE, row.names = 1)
count_matrix <- counts[, 6:ncol(counts)] # first 5 remaining cols are featureCounts metadata

format_sample_id <- function(colname) {
  colname |>
    str_remove("^output\\.hisat2_alignment\\.") |>
    str_remove("_align_sorted_markdup\\.bam$") |>
    str_replace_all("\\.", "-")
}
colnames(count_matrix) <- sapply(colnames(count_matrix), format_sample_id)

## ---- Load metadata ---------------------------------------------------------

sample_map <- tibble::tibble(
  sample_col = colnames(count_matrix),
  index_pair = str_extract(colnames(count_matrix), "[ACGT]+-[ACGT]+")
)

# janitor::clean_names() renders "i5barcode_NovaSeqV1.5" differently across
# janitor versions (2.2.1 splits the camel case into i5barcode_nova_seq_v1_5),
# so resolve the i5 column by pattern rather than hardcoding one spelling.
meta_raw <- read.csv(meta_file) |> janitor::clean_names()
i5_col <- grep("^i5barcode", names(meta_raw), value = TRUE)
if (length(i5_col) != 1) {
  stop(glue(
    "Expected exactly one i5barcode column in {meta_file}, found: {paste(i5_col, collapse = ', ')}"
  ))
}

meta <- meta_raw |>
  dplyr::mutate(
    index_pair = paste0(i7barcode, "-", .data[[i5_col]]),
    line = str_extract(sample, "^NPC[0-9]+"),
    sex = str_trim(sex),
    condition = factor(condition)
  ) |>
  dplyr::inner_join(sample_map, by = "index_pair") |>
  tibble::column_to_rownames("sample_col")

if (nrow(meta) == 0) {
  stop(glue("No samples in {meta_file} matched count matrix columns in {counts_file}"))
}

count_matrix <- count_matrix[, rownames(meta)]

## ---- QC MDS/PCA across all samples -----------------------------------------
## voom works on log-CPM, so the all-sample ordination uses TMM-normalised
## log-CPM here rather than DESeq2's VST.

dge_all <- DGEList(counts = count_matrix, samples = meta)
dge_all <- calcNormFactors(dge_all, method = "TMM")
keep_all <- filterByExpr(dge_all, group = meta$condition)
dge_all <- dge_all[keep_all, , keep.lib.sizes = FALSE]
logcpm_all <- cpm(dge_all, log = TRUE, prior.count = 3)

pca_from_logcpm <- function(logcpm, coldata, colour_var, labels, title) {
  # Match DESeq2's plotPCA(): top 500 most variable genes, samples in rows.
  rv <- matrixStats::rowVars(as.matrix(logcpm))
  select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
  pca <- prcomp(t(logcpm[select, , drop = FALSE]))
  percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2))
  df <- tibble::tibble(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    group = as.character(coldata[[colour_var]]),
    sample_label = as.character(labels)
  )
  ggplot(df, aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = sample_label), size = 3, max.overlaps = 20) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    labs(title = title, color = colour_var) +
    theme_bw()
}

pca_vars_all <- c("condition", "line", "sex") |> set_names()
pca_plots_all <- map(pca_vars_all, function(v) {
  pca_from_logcpm(logcpm_all, meta, v, meta$sample_id, glue("PCA (log-CPM) - all samples - {v}"))
})
pdf(file.path("results", "pca_plots_all_limma_voom.pdf"), width = 12, height = 6)
print(pca_plots_all)
dev.off()

## ---- Helpers ---------------------------------------------------------------

safe_filename <- function(x) {
  x <- stringr::str_replace_all(x, "\\+", "plus")
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9_-]", "_")
  tolower(x)
}

annotate_gene_symbols <- function(df) {
  gene_ids_nover <- sub("\\.\\d+$", "", df$gene)
  species <- NULL
  keytype <- "ENSEMBL"
  if (any(grepl("^ENSMUSG", gene_ids_nover))) {
    species <- "mouse"
  } else if (any(grepl("^ENSG", gene_ids_nover))) {
    species <- "human"
  } else if (any(grepl("^FBgn", gene_ids_nover))) {
    species <- "drosophila"
    keytype <- "FLYBASE"
  }

  df$gene_symbol <- NA_character_
  if (!is.null(species)) {
    OrgDb <- NULL
    if (species == "mouse" && requireNamespace("org.Mm.eg.db", quietly = TRUE)) OrgDb <- get("org.Mm.eg.db", envir = asNamespace("org.Mm.eg.db"))
    if (species == "human" && requireNamespace("org.Hs.eg.db", quietly = TRUE)) OrgDb <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))
    if (species == "drosophila" && requireNamespace("org.Dm.eg.db", quietly = TRUE)) OrgDb <- get("org.Dm.eg.db", envir = asNamespace("org.Dm.eg.db"))
    if (!is.null(OrgDb)) {
      map_df <- AnnotationDbi::select(OrgDb, keys = unique(gene_ids_nover), keytype = keytype, columns = c("SYMBOL"))
      if (!is.null(map_df) && nrow(map_df) > 0) {
        names(map_df)[names(map_df) == keytype] <- "gene_nover"
        names(map_df)[names(map_df) == "SYMBOL"] <- "gene_symbol_mapped"
        df$gene_nover <- gene_ids_nover
        df <- dplyr::left_join(df, map_df, by = "gene_nover")
        df$gene_symbol <- df$gene_symbol_mapped
        df$gene_nover <- NULL
        df$gene_symbol_mapped <- NULL
      }
    }
  }
  other_cols <- setdiff(colnames(df), c("gene", "gene_symbol"))
  df[, c("gene", "gene_symbol", other_cols), drop = FALSE]
}

# DESeq2's collapseReplicates() sums counts across technical replicates. Do the
# same here so `collapse_replicates = TRUE` means the same thing in both scripts.
collapse_count_replicates <- function(mat, groupby) {
  groupby <- as.character(groupby)
  collapsed <- t(rowsum(t(mat), group = groupby, reorder = FALSE))
  collapsed
}

## ---- Read line-group comparisons config ------------------------------------

comparisons_config <- yaml::read_yaml(comparisons_config_path)
stopifnot(
  "comparisons_config must define line_groups" = !is.null(comparisons_config$line_groups),
  "comparisons_config must define comparisons" = !is.null(comparisons_config$comparisons)
)

line_groups <- comparisons_config$line_groups
comparisons <- comparisons_config$comparisons
run_variants <- comparisons_config$run_variants
collapse_options <- if (!is.null(run_variants$collapse_replicates)) unlist(run_variants$collapse_replicates) else FALSE
male_options <- if (!is.null(run_variants$include_male_samples)) unlist(run_variants$include_male_samples) else TRUE

## ---- Run each comparison x run_variant combination -------------------------

manifest_rows <- list()

for (cmp in comparisons) {
  cmp_name <- cmp$name
  group_a <- cmp$group_a
  group_b <- cmp$group_b
  lines_a <- line_groups[[group_a]]
  lines_b <- line_groups[[group_b]]

  if (is.null(lines_a) || is.null(lines_b)) {
    warning(glue("Skipping comparison {cmp_name}: unknown line group(s) {group_a}/{group_b}"))
    next
  }

  for (collapse_replicates in collapse_options) {
    for (include_male_samples in male_options) {
      variant_label <- safe_filename(glue("collapse_{collapse_replicates}_male_{include_male_samples}"))
      message(glue("Running {cmp_name} [{variant_label}]"))

      meta_sub <- meta |>
        dplyr::filter(line %in% c(lines_a, lines_b))
      if (!isTRUE(include_male_samples)) {
        meta_sub <- meta_sub |> dplyr::filter(!str_detect(sex, regex("^male$", ignore_case = TRUE)))
      }
      meta_sub <- meta_sub |>
        dplyr::mutate(line_group = factor(ifelse(line %in% lines_a, group_a, group_b), levels = c(group_a, group_b)))

      if (nrow(meta_sub) < 2 || any(table(meta_sub$line_group) == 0)) {
        warning(glue("Skipping {cmp_name} [{variant_label}]: one group has zero samples after filtering"))
        next
      }

      counts_sub <- count_matrix[, rownames(meta_sub), drop = FALSE]
      sample_labels <- as.character(meta_sub$sample_id)
      group_vec <- meta_sub$line_group

      if (isTRUE(collapse_replicates)) {
        collapse_group <- str_remove(meta_sub$replicates, "-[0-9]+$")
        counts_sub <- collapse_count_replicates(counts_sub, collapse_group)
        # One row of metadata per collapsed group; line_group is constant within
        # a replicate group, so taking the first row per group is well defined.
        keep_first <- !duplicated(collapse_group)
        group_vec <- factor(
          as.character(meta_sub$line_group[keep_first]),
          levels = c(group_a, group_b)
        )
        # Label the pooled column by its collapsed group name rather than by the
        # first technical replicate, matching the deseq2 script's relabelling.
        sample_labels <- colnames(counts_sub)
      }

      if (any(table(group_vec) < 1) || nlevels(droplevels(group_vec)) < 2) {
        warning(glue("Skipping {cmp_name} [{variant_label}]: fewer than two groups after collapsing"))
        next
      }

      fit_result <- tryCatch(
        {
          dge <- DGEList(counts = counts_sub, group = group_vec)
          keep <- filterByExpr(dge, group = group_vec)
          dge <- dge[keep, , keep.lib.sizes = FALSE]
          dge <- calcNormFactors(dge, method = "TMM")

          # No-intercept design + explicit contrast so that a positive logFC
          # means "up in group_a", matching DESeq2's
          # contrast = c("line_group", group_a, group_b).
          design <- model.matrix(~ 0 + group_vec)
          colnames(design) <- levels(group_vec)

          voom_path <- file.path(out_dir, cmp_name, variant_label, "voom_mean_variance.pdf")
          dir.create(dirname(voom_path), showWarnings = FALSE, recursive = TRUE)
          pdf(voom_path, width = 7, height = 6)
          v <- voom(dge, design, plot = TRUE)
          dev.off()

          contrast_expr <- glue("{group_a}-{group_b}")
          contrast_matrix <- makeContrasts(contrasts = contrast_expr, levels = design)
          fit <- lmFit(v, design)
          fit <- contrasts.fit(fit, contrast_matrix)
          fit <- eBayes(fit)
          list(fit = fit, v = v, dge = dge, design = design)
        },
        error = function(e) {
          warning(glue("limma-voom fit failed for {cmp_name} [{variant_label}]: {e$message}"))
          NULL
        }
      )
      if (is.null(fit_result)) next

      fit <- fit_result$fit
      v <- fit_result$v
      dge <- fit_result$dge

      res_df <- topTable(fit, coef = 1, number = Inf, sort.by = "P")
      res_df$gene <- rownames(res_df)
      res_df <- annotate_gene_symbols(res_df)

      combo_dir <- file.path(out_dir, cmp_name, variant_label)
      dir.create(combo_dir, showWarnings = FALSE, recursive = TRUE)

      # adj.P.Val < 0.1 mirrors the DESeq2 script's padj < 0.1 cut so the two
      # manifests' significance counts are on the same footing.
      n_sig <- sum(!is.na(res_df$adj.P.Val) & res_df$adj.P.Val < 0.1)
      base_name <- safe_filename(glue("deg_{group_a}_v_{group_b}_{variant_label}_{n_sig}"))
      results_path <- file.path(combo_dir, glue("{base_name}.csv"))
      fit_path <- file.path(combo_dir, glue("{base_name}_efit.rds"))
      logcpm_path <- file.path(combo_dir, glue("{base_name}_logcpm.csv"))
      pca_path <- file.path(combo_dir, glue("{base_name}_pca.pdf"))

      write.csv(res_df, file = results_path, row.names = FALSE)
      saveRDS(fit, file = fit_path)
      write.csv(cpm(dge, log = TRUE, prior.count = 3), logcpm_path)

      pca_plot <- tryCatch(
        pca_from_logcpm(
          v$E,
          data.frame(line_group = group_vec),
          "line_group",
          sample_labels,
          glue("PCA (log-CPM) - {cmp_name} - {variant_label}")
        ),
        error = function(e) {
          warning(glue("PCA failed for {cmp_name} [{variant_label}]: {e$message}"))
          NULL
        }
      )
      if (!is.null(pca_plot)) {
        pdf(pca_path, width = 8, height = 6)
        print(pca_plot)
        dev.off()
      }

      manifest_rows[[length(manifest_rows) + 1]] <- tibble::tibble(
        comparison = cmp_name,
        group_a = group_a,
        group_b = group_b,
        collapse_replicates = collapse_replicates,
        include_male_samples = include_male_samples,
        n_samples = ncol(dge),
        n_group_a = sum(as.character(group_vec) == group_a),
        n_group_b = sum(as.character(group_vec) == group_b),
        n_genes_after_filter = nrow(dge),
        n_sig_adj_p_0_1 = n_sig,
        results_csv = results_path,
        efit_rds = fit_path
      )
    }
  }
}

manifest <- if (length(manifest_rows) > 0) dplyr::bind_rows(manifest_rows) else tibble::tibble()
readr::write_csv(manifest, file.path(out_dir, "limma_voom_comparisons_manifest.csv"))

## ---- Primary-comparison copies for downstream rules ------------------------
## Mirrors the deseq2 script: the first comparison, run without replicate
## collapsing and with all samples, is copied to a fixed path.

primary_cmp <- comparisons[[1]]$name
primary_row <- manifest |>
  dplyr::filter(comparison == primary_cmp, !collapse_replicates, include_male_samples)
if (nrow(primary_row) > 0) {
  file.copy(primary_row$efit_rds[1], file.path(out_dir, "efit.rds"), overwrite = TRUE)
  file.copy(primary_row$results_csv[1], file.path(out_dir, "limma_voom_results.csv"), overwrite = TRUE)
} else {
  warning(glue("Primary combination not found in manifest for {primary_cmp}; efit.rds/limma_voom_results.csv not written"))
}
