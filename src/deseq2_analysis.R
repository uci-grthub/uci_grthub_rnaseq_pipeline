#!/usr/bin/env Rscript
# DESeq2 analysis script
# Usage: Rscript deseq2_analysis.R <counts> metadata.csv output_dir comparisons_config.yaml
#
# <counts> may be any of three shapes, detected by file extension:
#   .txt   featureCounts output (HISAT2 branch)
#   .csv   plain gene x sample matrix from src/salmon_gene_counts.R
#   .rds   gene-level SummarizedExperiment from src/salmon_gene_counts_tximeta.R
#
# The .rds path is preferred for the Salmon branch: the dds inherits rowRanges
# and transcriptome provenance from tximeta, so results carry genomic
# coordinates without a second lookup.
#
# Comparisons are driven by config.deseq2_comparisons.yaml:
# `group_column` names the metadata column samples are grouped on (e.g. condition,
# or line for the NPC line-ID studies), `groups` maps a group name to the values of
# that column it covers, each `comparisons` entry contrasts two groups, and each is
# run under every combination of `run_variants`
# (collapse_replicates x include_male_samples).

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
default_comparisons_config <- "config.deseq2_comparisons.yaml"

counts_file <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_counts
meta_file <- if (length(args) >= 2 && nzchar(args[2])) args[2] else default_meta
out_dir <- if (length(args) >= 3 && nzchar(args[3])) args[3] else default_out
comparisons_config_path <- if (length(args) >= 4 && nzchar(args[4])) args[4] else default_comparisons_config

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

## ---- Load counts -----------------------------------------------------------

format_sample_id <- function(colname) {
  colname |>
    str_remove("^output\\.hisat2_alignment\\.") |>
    str_remove("_align_sorted_markdup\\.bam$") |>
    str_replace_all("\\.", "-")
}

gse <- NULL

if (grepl("\\.rds$", counts_file, ignore.case = TRUE)) {
  # tximeta branch: a gene-level SummarizedExperiment. Counts were summarised
  # with countsFromAbundance="lengthScaledTPM", so the length correction is
  # already inside them and DESeqDataSet() will not apply offsets on top --
  # it reports "using just counts from tximeta". That is what we want here,
  # because collapseReplicates() below sums assays, and summing per-sample
  # normalisation factors would be meaningless.
  message(glue("Reading tximeta SummarizedExperiment from {counts_file}"))
  gse <- readRDS(counts_file)
  if (!inherits(gse, "SummarizedExperiment")) {
    stop(glue("{counts_file} is not a SummarizedExperiment"))
  }
  count_matrix <- as.data.frame(round(SummarizedExperiment::assay(gse, "counts")))
  gene_lengths <- NULL

  meta_raw <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (all(c("sample", "sample_id") %in% names(meta_raw)) &&
      all(colnames(count_matrix) %in% meta_raw$sample_id)) {
    translated <- meta_raw$sample[match(colnames(count_matrix), meta_raw$sample_id)]
    colnames(count_matrix) <- translated
    colnames(gse) <- translated
  }
} else if (grepl("\\.csv$", counts_file, ignore.case = TRUE)) {
  # Salmon branch: a plain gene x sample matrix from src/salmon_gene_counts.R.
  # Already integer and length-corrected (countsFromAbundance="lengthScaledTPM"),
  # so it needs no offset here.
  message(glue("Reading CSV count matrix from {counts_file}"))
  count_matrix <- read.csv(counts_file, header = TRUE, row.names = 1, check.names = FALSE)
  gene_lengths <- NULL

  # That script labels columns with the metadata sample_id, which carries no
  # index pair. Translate back to the FASTQ name so the barcode join below
  # behaves identically for both branches.
  meta_raw <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (all(c("sample", "sample_id") %in% names(meta_raw)) &&
      all(colnames(count_matrix) %in% meta_raw$sample_id)) {
    colnames(count_matrix) <- meta_raw$sample[match(colnames(count_matrix), meta_raw$sample_id)]
  }
} else {
  # featureCounts branch
  counts <- read.table(counts_file, header = TRUE, row.names = 1)
  count_matrix <- counts[, 6:ncol(counts)] # first 5 remaining cols are featureCounts metadata
  gene_lengths <- counts[, 5]
  colnames(count_matrix) <- sapply(colnames(count_matrix), format_sample_id)
}

## ---- Load metadata -----------------------------------------------------------

sample_map <- tibble::tibble(
  sample_col = colnames(count_matrix),
  index_pair = str_extract(colnames(count_matrix), "[ACGT]+-[ACGT]+")
)

meta <- read.csv(meta_file) |>
  janitor::clean_names()

# janitor::clean_names() renders "i5barcode_NovaSeqV1.5" differently across
# versions (i5barcode_novaseqv1_5 vs i5barcode_nova_seq_v1_5), so resolve the
# barcode columns by prefix instead of hardcoding either spelling.
pick_col <- function(prefix) {
  hits <- grep(paste0("^", prefix), names(meta), value = TRUE)
  if (length(hits) == 0) stop(glue("Metadata {meta_file} has no column starting with '{prefix}'"))
  hits[1]
}
i7_col <- pick_col("i7barcode")
i5_col <- pick_col("i5barcode")

# sex is optional in the metadata sheet, and blank entries mean "unknown"
# rather than a distinct group. Normalise both to NA so the male filter and
# the QC PCA below can treat them the same way.
if (!"sex" %in% names(meta)) {
  meta$sex <- NA_character_
}

meta <- meta |>
  dplyr::mutate(
    index_pair = paste0(.data[[i7_col]], "-", .data[[i5_col]]),
    line = str_extract(sample, "^NPC[0-9]+"),
    sex = dplyr::na_if(str_trim(sex), ""),
    condition = factor(condition)
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

## ---- dds construction --------------------------------------------------------

# Build from the SummarizedExperiment when there is one, so rowRanges and the
# tximeta provenance travel with the dds; fall back to the plain matrix
# otherwise. Both produce the same counts.
make_dds <- function(cols, coldata, design) {
  if (!is.null(gse)) {
    sub <- gse[, cols, drop = FALSE]
    SummarizedExperiment::colData(sub) <- S4Vectors::DataFrame(coldata)
    DESeqDataSet(sub, design = design)
  } else {
    DESeqDataSetFromMatrix(
      countData = count_matrix[, cols, drop = FALSE],
      colData = coldata,
      design = design
    )
  }
}

## ---- QC PCA across all samples -----------------------------------------------

dds_all <- make_dds(rownames(meta), meta, ~condition)
dds_all <- DESeq(dds_all)
vsd_all <- vst(dds_all, blind = FALSE)
# Only colour by variables that are actually populated -- plotPCA errors on a
# column that is entirely NA, and an all-blank column carries no information.
pca_vars_all <- c("condition", "line", "sex") |>
  keep(function(v) v %in% names(colData(vsd_all)) && any(!is.na(colData(vsd_all)[[v]]))) |>
  set_names()
pca_plots_all <- map(pca_vars_all, function(v) {
  p <- plotPCA(vsd_all, intgroup = v, returnData = TRUE)
  p$sample_label <- as.character(colData(vsd_all)$sample_id)
  percentVar <- round(100 * attr(p, "percentVar"))
  ggplot(p, aes(PC1, PC2, color = !!sym(v))) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = sample_label), size = 3, max.overlaps = 20) +
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

## ---- Read comparison-group config -----------------------------------------------------------

comparisons_config <- yaml::read_yaml(comparisons_config_path)
# line_groups is the older key name from the NPC line studies, where the
# grouping column was always `line`; `groups` + `group_column` generalise it.
sample_groups <- comparisons_config$groups %||% comparisons_config$line_groups
group_column <- comparisons_config$group_column %||% "line"
stopifnot(
  "comparisons_config must define groups (or line_groups)" = !is.null(sample_groups),
  "comparisons_config must define comparisons" = !is.null(comparisons_config$comparisons)
)
if (!group_column %in% names(meta)) {
  stop(glue("group_column '{group_column}' is not a column in {meta_file}"))
}

comparisons <- comparisons_config$comparisons
run_variants <- comparisons_config$run_variants
collapse_options <- if (!is.null(run_variants$collapse_replicates)) unlist(run_variants$collapse_replicates) else FALSE
male_options <- if (!is.null(run_variants$include_male_samples)) unlist(run_variants$include_male_samples) else TRUE

## ---- Run each comparison x run_variant combination -----------------------------------------------------------

manifest_rows <- list()

for (cmp in comparisons) {
  cmp_name <- cmp$name
  group_a <- cmp$group_a
  group_b <- cmp$group_b
  levels_a <- sample_groups[[group_a]]
  levels_b <- sample_groups[[group_b]]

  if (is.null(levels_a) || is.null(levels_b)) {
    warning(glue("Skipping comparison {cmp_name}: unknown group(s) {group_a}/{group_b}"))
    next
  }

  for (collapse_replicates in collapse_options) {
    for (include_male_samples in male_options) {
      variant_label <- safe_filename(glue("collapse_{collapse_replicates}_male_{include_male_samples}"))
      message(glue("Running {cmp_name} [{variant_label}]"))

      meta_sub <- meta |>
        dplyr::filter(as.character(.data[[group_column]]) %in% c(levels_a, levels_b))
      if (!isTRUE(include_male_samples)) {
        # Unknown sex (NA) is kept: dropping male samples should not also drop
        # every sample whose sex was never recorded.
        meta_sub <- meta_sub |>
          dplyr::filter(is.na(sex) | !str_detect(sex, regex("^male$", ignore_case = TRUE)))
      }
      meta_sub <- meta_sub |>
        dplyr::mutate(comparison_group = factor(
          ifelse(as.character(.data[[group_column]]) %in% levels_a, group_a, group_b),
          levels = c(group_a, group_b)
        ))

      if (nrow(meta_sub) < 2 || any(table(meta_sub$comparison_group) == 0)) {
        warning(glue("Skipping {cmp_name} [{variant_label}]: one group has zero samples after filtering"))
        next
      }

      dds <- make_dds(rownames(meta_sub), meta_sub, ~comparison_group)

      if (isTRUE(collapse_replicates)) {
        collapse_group <- str_remove(meta_sub$replicates, "-[0-9]+$")
        dds <- collapseReplicates(dds, groupby = collapse_group)
        # collapseReplicates() keeps the first technical replicate's colData
        # (e.g. sample_id "Line 80-1") even though the column is now pooled
        # across all replicates for that group, so re-point sample_id at the
        # actual collapsed group name ("Line 80") to avoid a misleading label.
        colData(dds)$sample_id <- colnames(dds)
      }

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

      n_sig <- sum(!is.na(res_df$padj) & res_df$padj < 0.1)
      base_name <- safe_filename(glue("deg_{group_a}_v_{group_b}_{variant_label}_{n_sig}"))
      results_path <- file.path(combo_dir, glue("{base_name}.csv"))
      dds_path <- file.path(combo_dir, glue("{base_name}_dds.rds"))
      norm_counts_path <- file.path(combo_dir, glue("{base_name}_normalized_counts.csv"))
      pca_path <- file.path(combo_dir, glue("{base_name}_pca.pdf"))

      write.csv(res_df, file = results_path, row.names = FALSE)
      saveRDS(dds, file = dds_path)
      write.csv(counts(dds, normalized = TRUE), norm_counts_path)

      pca_plot <- tryCatch(
        {
          vsd <- vst(dds, blind = FALSE)
          p <- plotPCA(vsd, intgroup = "comparison_group", returnData = TRUE)
          p$sample_label <- as.character(colData(vsd)$sample_id)
          percentVar <- round(100 * attr(p, "percentVar"))
          ggplot(p, aes(PC1, PC2, color = comparison_group)) +
            geom_point(size = 3) +
            geom_text_repel(aes(label = sample_label), size = 3, max.overlaps = 20) +
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
        n_samples = ncol(dds),
        n_group_a = sum(as.character(colData(dds)$comparison_group) == group_a),
        n_group_b = sum(as.character(colData(dds)$comparison_group) == group_b),
        n_sig_padj_0_1 = n_sig,
        results_csv = results_path,
        dds_rds = dds_path
      )
    }
  }
}

manifest <- if (length(manifest_rows) > 0) dplyr::bind_rows(manifest_rows) else tibble::tibble()
readr::write_csv(manifest, file.path(out_dir, "deseq2_comparisons_manifest.csv"))

## ---- Legacy single dds.rds / deseq2_results.csv for downstream rules -----------------------------------------------------------
## Points at the first comparison, run without replicate collapsing and with all
## samples, since rule deseq2 declares one fixed dds.rds / deseq2_results.csv.

primary_cmp <- comparisons[[1]]$name
primary_row <- manifest |>
  dplyr::filter(comparison == primary_cmp, !collapse_replicates, include_male_samples)
if (nrow(primary_row) > 0) {
  file.copy(primary_row$dds_rds[1], file.path(out_dir, "dds.rds"), overwrite = TRUE)
  file.copy(primary_row$results_csv[1], file.path(out_dir, "deseq2_results.csv"), overwrite = TRUE)
} else {
  warning(glue("Primary combination not found in manifest for {primary_cmp}; dds.rds/deseq2_results.csv not written"))
}
