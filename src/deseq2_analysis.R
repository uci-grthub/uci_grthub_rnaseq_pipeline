#!/usr/bin/env Rscript
# DESeq2 analysis script
# Usage: Rscript deseq2_analysis.R counts.txt metadata.csv output_dir comparisons_config.yaml
#
# Comparisons are driven by proj_src/deseq2_comparisons.yaml (see proj_src/notes.md):
# samples are selected into `groups` by condition/treatment metadata filters,
# each `comparisons` entry contrasts two groups (group_a vs group_b, group_b as
# the reference level), and each is run under every `run_variants` value
# (include_male_samples).

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(glue)
  library(AnnotationDbi)
  library(yaml)
  library(iSEEde)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)

default_counts <- "output/feature_count/all_samples_counts.txt"
default_meta <- "metadata/metadata.csv"
default_out <- "output/deseq2"
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
gene_lengths <- counts[, 5]

format_sample_id <- function(colname) {
  colname |>
    str_remove("^output\\.hisat2_alignment\\.") |>
    str_remove("_align_sorted_markdup\\.bam$") |>
    str_replace_all("\\.", "-")
}
colnames(count_matrix) <- sapply(colnames(count_matrix), format_sample_id)

## ---- Load metadata -----------------------------------------------------------

sample_map <- tibble::tibble(
  sample_col = colnames(count_matrix),
  index_pair = str_extract(colnames(count_matrix), "[ACGT]+-[ACGT]+")
)

meta <- read.csv(meta_file) |>
  janitor::clean_names() |>
  dplyr::mutate(
    index_pair = glue("{i7barcode}-{i5barcode_novaseqv1_5}") |> as.character(),
    line = str_extract(sample, "[0-9]+[A-Z][0-9]+"),  # e.g. "CNiMG: 26C5: ..." -> "26C5"
    sex = str_trim(sex),
    treatment = str_trim(treatment),
    condition = factor(str_trim(condition))
  ) |>
  dplyr::inner_join(sample_map, by = "index_pair") |>
  tibble::column_to_rownames("sample_col")

if (nrow(meta) == 0) {
  stop(glue("No samples in {meta_file} matched count matrix columns in {counts_file}"))
}

count_matrix <- count_matrix[, rownames(meta)]

count_matrix |>
  tibble::rownames_to_column("ensgene") |>
  readr::write_csv(file.path(out_dir, "all_sample_counts.csv"))

## ---- QC PCA across all samples -----------------------------------------------

dds_all <- DESeqDataSetFromMatrix(countData = count_matrix, colData = meta, design = ~condition)
dds_all <- DESeq(dds_all)
vsd_all <- vst(dds_all, blind = FALSE)
pca_vars_all <- c("condition", "line", "sex") |> set_names()
pca_plots_all <- map(pca_vars_all, function(v) {
  p <- plotPCA(vsd_all, intgroup = v, returnData = TRUE)
  percentVar <- round(100 * attr(p, "percentVar"))
  ggplot(p, aes(PC1, PC2, color = !!sym(v))) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = name), size = 3, max.overlaps = 20) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    labs(title = glue("PCA - all samples - {v}")) +
    theme_bw()
})
pdf(file.path("results", "pca_plots_all.pdf"), width = 12, height = 6)
print(pca_plots_all)
dev.off()

## ---- Helpers -----------------------------------------------------------

safe_filename <- function(x) {
  x <- stringr::str_replace_all(x, "\\+", "plus")
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9_-]", "_")
  tolower(x)
}

# Return a logical vector over meta rows matching every column filter in `spec`
# (a named list such as list(condition = "CN", treatment = "AB")).
group_filter <- function(meta, spec) {
  keep <- rep(TRUE, nrow(meta))
  for (col in names(spec)) {
    if (!col %in% colnames(meta)) {
      stop(glue("Group filter references unknown metadata column '{col}'"))
    }
    keep <- keep & (as.character(meta[[col]]) == as.character(spec[[col]]))
  }
  keep
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

## ---- Read condition/treatment comparisons config -----------------------------------------------------------

comparisons_config <- yaml::read_yaml(comparisons_config_path)
stopifnot(
  "comparisons_config must define groups" = !is.null(comparisons_config$groups),
  "comparisons_config must define comparisons" = !is.null(comparisons_config$comparisons)
)

groups <- comparisons_config$groups
comparisons <- comparisons_config$comparisons
run_variants <- comparisons_config$run_variants
male_options <- if (!is.null(run_variants$include_male_samples)) unlist(run_variants$include_male_samples) else TRUE

## ---- Run each comparison x run_variant combination -----------------------------------------------------------

manifest_rows <- list()

for (cmp in comparisons) {
  cmp_name <- cmp$name
  group_a <- cmp$group_a
  group_b <- cmp$group_b
  spec_a <- groups[[group_a]]
  spec_b <- groups[[group_b]]

  if (is.null(spec_a) || is.null(spec_b)) {
    warning(glue("Skipping comparison {cmp_name}: unknown group(s) {group_a}/{group_b}"))
    next
  }

  for (include_male_samples in male_options) {
    variant_label <- safe_filename(glue("male_{include_male_samples}"))
    message(glue("Running {cmp_name} [{variant_label}]"))

    meta_sub <- meta
    meta_sub$comparison_group <- NA_character_
    meta_sub$comparison_group[group_filter(meta_sub, spec_a)] <- group_a
    meta_sub$comparison_group[group_filter(meta_sub, spec_b)] <- group_b
    meta_sub <- meta_sub[!is.na(meta_sub$comparison_group), , drop = FALSE]

    if (!isTRUE(include_male_samples)) {
      meta_sub <- meta_sub[!str_detect(meta_sub$sex, regex("^male$", ignore_case = TRUE)), , drop = FALSE]
    }
    meta_sub$comparison_group <- factor(meta_sub$comparison_group, levels = c(group_a, group_b))

    if (nrow(meta_sub) < 2 || any(table(meta_sub$comparison_group) == 0)) {
      warning(glue("Skipping {cmp_name} [{variant_label}]: one group has zero samples after filtering"))
      next
    }

    counts_sub <- count_matrix[, rownames(meta_sub), drop = FALSE]
    dds <- DESeqDataSetFromMatrix(countData = counts_sub, colData = meta_sub, design = ~comparison_group)

    dds <- tryCatch(DESeq(dds), error = function(e) {
      warning(glue("DESeq failed for {cmp_name} [{variant_label}]: {e$message}"))
      NULL
    })
    if (is.null(dds)) next

    res <- tryCatch(
      results(dds, contrast = c("comparison_group", group_a, group_b)),
      error = function(e) {
        warning(glue("results() failed for {cmp_name} [{variant_label}]: {e$message}"))
        NULL
      }
    )
    if (is.null(res)) next

    combo_dir <- file.path(out_dir, cmp_name, variant_label)
    dir.create(combo_dir, showWarnings = FALSE, recursive = TRUE)

    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df <- annotate_gene_symbols(res_df)
    if ("padj" %in% colnames(res_df)) res_df <- res_df[order(res_df$padj, na.last = TRUE), , drop = FALSE]
    write.csv(res_df, file = file.path(combo_dir, "results.csv"), row.names = FALSE)
    saveRDS(dds, file = file.path(combo_dir, "dds.rds"))
    write.csv(counts(dds, normalized = TRUE), file.path(combo_dir, "normalized_counts.csv"))

    pca_plot <- tryCatch(
      {
        vsd <- vst(dds, blind = FALSE)
        p <- plotPCA(vsd, intgroup = "comparison_group", returnData = TRUE)
        percentVar <- round(100 * attr(p, "percentVar"))
        ggplot(p, aes(PC1, PC2, color = comparison_group)) +
          geom_point(size = 3) +
          geom_text_repel(aes(label = name), size = 3, max.overlaps = 20) +
          xlab(paste0("PC1: ", percentVar[1], "% variance")) +
          ylab(paste0("PC2: ", percentVar[2], "% variance")) +
          labs(title = glue("PCA - {cmp_name} - {variant_label}")) +
          theme_bw()
      },
      error = function(e) {
        warning(glue("PCA failed for {cmp_name} [{variant_label}]: {e$message}"))
        NULL
      }
    )
    if (!is.null(pca_plot)) {
      pdf(file.path(combo_dir, "pca.pdf"), width = 8, height = 6)
      print(pca_plot)
      dev.off()
    }

    manifest_rows[[length(manifest_rows) + 1]] <- tibble::tibble(
      comparison = cmp_name,
      group_a = group_a,
      group_b = group_b,
      include_male_samples = include_male_samples,
      n_samples = ncol(dds),
      n_group_a = sum(as.character(colData(dds)$comparison_group) == group_a),
      n_group_b = sum(as.character(colData(dds)$comparison_group) == group_b),
      results_csv = file.path(combo_dir, "results.csv"),
      dds_rds = file.path(combo_dir, "dds.rds")
    )
  }
}

manifest <- if (length(manifest_rows) > 0) dplyr::bind_rows(manifest_rows) else tibble::tibble()
readr::write_csv(manifest, file.path(out_dir, "deseq2_comparisons_manifest.csv"))

## ---- Legacy single dds.rds / deseq2_results.csv for downstream rules -----------------------------------------------------------
## Points at the first comparison, run without replicate collapsing and with all samples,
## since generate_isee_app / the shiny app still expect one fixed dds.rds.

primary_cmp <- comparisons[[1]]$name
primary_dir <- file.path(out_dir, primary_cmp, "male_true")
if (dir.exists(primary_dir)) {
  file.copy(file.path(primary_dir, "dds.rds"), file.path(out_dir, "dds.rds"), overwrite = TRUE)
  file.copy(file.path(primary_dir, "results.csv"), file.path(out_dir, "deseq2_results.csv"), overwrite = TRUE)
} else {
  warning(glue("Primary combination directory not found: {primary_dir}; dds.rds/deseq2_results.csv not written"))
}
