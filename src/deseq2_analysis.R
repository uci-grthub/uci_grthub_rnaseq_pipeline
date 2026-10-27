#!/usr/bin/env Rscript
# DESeq2 analysis for the MaccF HTR2A knockdown study.
# Usage: Rscript deseq2_analysis.R counts.txt metadata.csv output_dir comparisons.yaml [ensdb_sqlite] [gtf]
#
# Design (src/de_comparisons.yaml, metadata/notes.md): 4 animals x 8 regions,
# one library each, sequenced on two lanes that HISAT2 aligns together, so the
# count columns are libraries (metadata `library`). Each comparison fits
# `~ region + <factor>` on a region set -- region is the blocking factor -- and
# contrasts each level of <factor> against its reference.
#
# The animal is the unit of inference (1 control, 2 COG301, 1 CPG201), so next
# to the model this writes descriptive per-region, per-animal log2FC for the
# target and reference genes, and cell-composition marker scores to flag
# dissection differences that could mimic a knockdown.
#
# Gene symbols come from the EnsDb when it exists, otherwise from the gene
# records of the GTF the counts were made against.

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(glue)
  library(yaml)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
arg_or <- function(i, default) if (length(args) >= i && nzchar(args[i])) args[i] else default

counts_file <- arg_or(1, "output/feature_count/macaque_samples_counts.txt")
meta_file <- arg_or(2, "metadata/metadata.csv")
out_dir <- arg_or(3, "output/deseq2/macaque")
comparisons_config_path <- arg_or(4, "src/de_comparisons.yaml")
ensdb_path <- arg_or(5, "")
gtf_path <- arg_or(6, "")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Shared with count_matrix.R so every output labels genes the same way
local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  source(file.path(if (length(f) == 1) dirname(normalizePath(f)) else "src", "gene_annotation.R"))
})

cfg <- yaml::read_yaml(comparisons_config_path)
stopifnot(
  "comparisons config must define region_sets" = !is.null(cfg$region_sets),
  "comparisons config must define comparisons" = !is.null(cfg$comparisons)
)

## ---- Gene annotation --------------------------------------------------------

read_gtf_genes <- function(path) {
  # Gene records only; a full rtracklayer import of a 1.4M-line GTF is slow
  lines <- readr::read_lines(path) |> purrr::keep(\(x) grepl("\tgene\t", x, fixed = TRUE))
  attr_of <- function(key) str_match(lines, glue('{key} "([^"]+)"'))[, 2]
  tibble(gene_id = attr_of("gene_id"), gene_name = attr_of("gene_name"), gene_biotype = attr_of("gene_biotype"))
}

gene_ann <- if (nzchar(ensdb_path) && file.exists(ensdb_path)) {
  read_gene_annotation(ensdb_path)
} else if (nzchar(gtf_path) && file.exists(gtf_path)) {
  message(glue("EnsDb not available; reading gene names from {gtf_path}"))
  read_gtf_genes(gtf_path)
} else {
  stop("Need an EnsDb or GTF to resolve gene symbols (arguments 5 and 6)")
}

symbol_to_id <- function(symbols) {
  ids <- gene_ann$gene_id[match(symbols, gene_ann$gene_name)]
  missing <- symbols[is.na(ids)]
  if (length(missing) > 0) warning(glue("Symbols not in annotation: {paste(missing, collapse = ', ')}"))
  set_names(ids, symbols)[!is.na(ids)]
}

## ---- Counts and metadata ----------------------------------------------------

fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
annotation_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
count_matrix <- as.matrix(fc[, setdiff(colnames(fc), annotation_cols), drop = FALSE])
rownames(count_matrix) <- fc$Geneid
# Column names are BAM paths: output/hisat2_alignment/<library>_align_sorted_markdup.bam.
# Lanes were aligned together, so each column is already one library.
colnames(count_matrix) <- basename(colnames(count_matrix)) |> str_remove("_align_sorted_markdup\\.bam$")

meta <- read.csv(meta_file, stringsAsFactors = FALSE)
required <- c("library", "condition", "nhp", "region")
missing_cols <- setdiff(required, colnames(meta))
if (length(missing_cols) > 0) stop(glue("Metadata {meta_file} lacks: {paste(missing_cols, collapse = ', ')}"))

# One row per library; lane-specific columns (sample, lane, sample_id) are dropped
meta <- meta |>
  dplyr::select(-dplyr::any_of(c("sample", "lane", "sample_id"))) |>
  dplyr::distinct(library, .keep_all = TRUE)

unmatched <- setdiff(colnames(count_matrix), meta$library)
if (length(unmatched) > 0) stop(glue("No metadata row for: {paste(unmatched, collapse = ', ')}"))
meta <- meta |>
  dplyr::filter(library %in% colnames(count_matrix)) |>
  dplyr::mutate(sample_id = library, region = factor(region, levels = cfg$region_sets$all)) |>
  tibble::column_to_rownames("library")
count_matrix <- count_matrix[, rownames(meta)]

dds_full <- DESeqDataSetFromMatrix(count_matrix, colData = meta, design = ~1)

## ---- QC PCA across all libraries --------------------------------------------

dds_full <- estimateSizeFactors(dds_full)
vsd_full <- vst(dds_full, blind = TRUE)
pca_plot <- function(vsd, v, title) {
  p <- plotPCA(vsd, intgroup = v, returnData = TRUE, ntop = 1000)
  p$sample_label <- colData(vsd)$sample_id
  pv <- round(100 * attr(p, "percentVar"))
  ggplot(p, aes(PC1, PC2, color = .data[[v]])) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = sample_label), size = 2.5, max.overlaps = 30) +
    labs(x = glue("PC1: {pv[1]}% variance"), y = glue("PC2: {pv[2]}% variance"), title = glue("{title} - {v}")) +
    theme_bw()
}
pdf(file.path(out_dir, "pca_all_libraries.pdf"), width = 9, height = 7)
walk(c("region", "condition", "nhp"), \(v) print(pca_plot(vsd_full, v, "PCA - all libraries")))
dev.off()

## ---- Descriptive tables: target, reference and marker genes ------------------

norm_counts <- counts(dds_full, normalized = TRUE)
lib_meta <- as.data.frame(colData(dds_full)) |>
  tibble::rownames_to_column("library_id") |>
  dplyr::select(library_id, condition, nhp, region)

target_ids <- symbol_to_id(unlist(cfg$target_genes))
reference_ids <- symbol_to_id(unlist(cfg$reference_genes))
gene_long <- function(ids) {
  norm_counts[ids, , drop = FALSE] |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    dplyr::mutate(gene_symbol = names(ids)) |>
    pivot_longer(-c(gene, gene_symbol), names_to = "library_id", values_to = "norm_count") |>
    dplyr::left_join(lib_meta, by = "library_id")
}

reference_level <- function(factor_name) {
  cmp <- purrr::detect(cfg$comparisons, \(x) x$factor == factor_name)
  if (is.null(cmp)) stop(glue("No comparison uses factor {factor_name}"))
  cmp$reference
}
control_group <- reference_level("condition")

# Per region, each treated animal against the control animal. Pseudocount 0.5
# keeps a zero count from producing an infinite ratio.
per_region <- gene_long(c(target_ids, reference_ids)) |>
  dplyr::group_by(gene, gene_symbol, region) |>
  dplyr::mutate(control_norm_count = norm_count[condition == control_group][1]) |>
  dplyr::ungroup() |>
  dplyr::filter(condition != control_group) |>
  dplyr::mutate(log2fc_vs_control = log2((norm_count + 0.5) / (control_norm_count + 0.5)))

reference_mean <- per_region |>
  dplyr::filter(gene_symbol %in% names(reference_ids)) |>
  dplyr::group_by(nhp, region) |>
  dplyr::summarise(reference_mean_log2fc = mean(log2fc_vs_control), .groups = "drop")
per_region <- per_region |>
  dplyr::left_join(reference_mean, by = c("nhp", "region")) |>
  dplyr::mutate(
    region_tier = map_chr(as.character(region), \(r) {
      tier <- names(cfg$region_sets)[map_lgl(cfg$region_sets, \(s) r %in% s)]
      setdiff(tier, "all")[1]
    }),
    log2fc_relative_to_reference = log2fc_vs_control - reference_mean_log2fc
  )
readr::write_csv(per_region, file.path(out_dir, "target_per_region_log2fc.csv"))

# Per animal over the primary regions: how many go down, and the mean log2FC
per_animal_summary <- per_region |>
  dplyr::filter(gene_symbol %in% names(target_ids), region %in% cfg$region_sets$primary) |>
  dplyr::group_by(gene, gene_symbol, condition, nhp) |>
  dplyr::summarise(
    n_regions = dplyr::n(),
    n_decreased = sum(log2fc_vs_control < 0),
    mean_log2fc = mean(log2fc_vs_control),
    n_decreased_relative = sum(log2fc_relative_to_reference < 0),
    mean_log2fc_relative = mean(log2fc_relative_to_reference),
    .groups = "drop"
  )
readr::write_csv(per_animal_summary, file.path(out_dir, "target_per_animal_summary.csv"))

# Baseline check: target expression in every animal and region
pdf(file.path(out_dir, "target_expression_by_region.pdf"), width = 10, height = 5)
gene_long(target_ids) |>
  dplyr::group_by(gene_symbol) |>
  dplyr::group_walk(\(d, key) print(
    ggplot(d, aes(region, norm_count, color = condition, shape = nhp)) +
      geom_point(size = 3, position = position_dodge(width = 0.5)) +
      scale_y_log10() +
      labs(y = "DESeq2-normalised count (log scale)", title = glue("{key$gene_symbol} by region")) +
      theme_bw()
  ))
dev.off()

# Cell-composition scores: mean log2(normalised count + 1) of each marker set
composition <- imap_dfr(cfg$composition_markers, \(symbols, set) {
  ids <- symbol_to_id(unlist(symbols))
  tibble(
    library_id = colnames(norm_counts),
    marker_set = set,
    n_genes = length(ids),
    score = colMeans(log2(norm_counts[ids, , drop = FALSE] + 1))
  )
}) |>
  dplyr::left_join(lib_meta, by = "library_id") |>
  dplyr::group_by(marker_set, region) |>
  dplyr::mutate(
    delta_vs_control = score - score[condition == control_group][1],
    flagged = condition != control_group & abs(delta_vs_control) > cfg$composition_flag_log2
  ) |>
  dplyr::ungroup()
readr::write_csv(composition, file.path(out_dir, "composition_scores.csv"))
if (any(composition$flagged)) {
  flagged <- composition |> dplyr::filter(flagged) |> glue::glue_data("{library_id}:{marker_set}")
  warning(glue("Composition differs from control by > {cfg$composition_flag_log2} log2: {paste(flagged, collapse = ', ')}"))
}

## ---- Model fits -------------------------------------------------------------

safe_filename <- function(x) tolower(str_replace_all(x, "[^A-Za-z0-9_-]", "_"))

manifest_rows <- list()
target_rows <- list()

for (cmp in cfg$comparisons) {
  regions <- cfg$region_sets[[cmp$regions]]
  if (is.null(regions)) stop(glue("Comparison {cmp$name}: unknown region set {cmp$regions}"))
  message(glue("Fitting {cmp$name}: ~ region + {cmp$factor} on {paste(regions, collapse = ', ')}"))

  dds <- dds_full[, colData(dds_full)$region %in% regions]
  colData(dds)$region <- droplevels(colData(dds)$region)
  levels_present <- unique(as.character(colData(dds)[[cmp$factor]]))
  if (!cmp$reference %in% levels_present) stop(glue("Comparison {cmp$name}: reference {cmp$reference} has no samples"))
  colData(dds)[[cmp$factor]] <- relevel(factor(colData(dds)[[cmp$factor]]), ref = cmp$reference)
  design(dds) <- as.formula(glue("~ region + {cmp$factor}"))

  # Drop genes with almost no reads; DESeq2's independent filtering does the rest
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  dds <- DESeq(dds, quiet = TRUE)

  cmp_dir <- file.path(out_dir, cmp$name)
  dir.create(cmp_dir, showWarnings = FALSE, recursive = TRUE)
  dds_path <- file.path(cmp_dir, "dds.rds")
  saveRDS(dds, dds_path)

  for (level in setdiff(levels(colData(dds)[[cmp$factor]]), cmp$reference)) {
    res <- results(dds, contrast = c(cmp$factor, level, cmp$reference))
    res_df <- as.data.frame(res) |>
      tibble::rownames_to_column("gene") |>
      dplyr::mutate(
        ci_low = log2FoldChange - qnorm(0.975) * lfcSE,
        ci_high = log2FoldChange + qnorm(0.975) * lfcSE
      ) |>
      annotate_gene_symbols(gene_ann) |>
      dplyr::arrange(padj)

    contrast_label <- glue("{level}_vs_{cmp$reference}")
    results_path <- file.path(cmp_dir, glue("{safe_filename(contrast_label)}.csv"))
    readr::write_csv(res_df, results_path)

    # The target is one pre-specified gene, so its raw p-value is the test;
    # padj is the genome-wide adjustment and is kept only for reference.
    target_rows[[length(target_rows) + 1]] <- res_df |>
      dplyr::filter(gene %in% c(target_ids, reference_ids)) |>
      dplyr::mutate(
        comparison = cmp$name, contrast = contrast_label, regions = paste(regions, collapse = ";"),
        role = ifelse(gene %in% target_ids, "target", "reference")
      )

    n_sig <- sum(!is.na(res_df$padj) & res_df$padj < 0.1)
    manifest_rows[[length(manifest_rows) + 1]] <- tibble(
      comparison = cmp$name,
      contrast = contrast_label,
      design = glue("~ region + {cmp$factor}"),
      regions = paste(regions, collapse = ";"),
      n_libraries = ncol(dds),
      n_group = sum(colData(dds)[[cmp$factor]] == level),
      n_reference = sum(colData(dds)[[cmp$factor]] == cmp$reference),
      n_genes_tested = sum(!is.na(res_df$padj)),
      n_sig_padj_0_1 = n_sig,
      results_csv = results_path,
      dds_rds = dds_path
    )
  }

  pdf(file.path(cmp_dir, "pca.pdf"), width = 9, height = 7)
  vsd <- vst(dds, blind = FALSE)
  walk(c("region", cmp$factor), \(v) print(pca_plot(vsd, v, glue("PCA - {cmp$name}"))))
  dev.off()
}

manifest <- dplyr::bind_rows(manifest_rows)
readr::write_csv(manifest, file.path(out_dir, "deseq2_comparisons_manifest.csv"))

target_results <- dplyr::bind_rows(target_rows) |>
  dplyr::select(comparison, contrast, regions, role, gene, gene_symbol, baseMean,
                log2FoldChange, lfcSE, ci_low, ci_high, stat, pvalue, padj)
readr::write_csv(target_results, file.path(out_dir, "target_gene_results.csv"))

## ---- Fixed-name outputs for downstream rules --------------------------------
## The first comparison's first contrast (COG301 vs Control on the primary regions).

primary <- manifest[1, ]
invisible(file.copy(primary$dds_rds, file.path(out_dir, "dds.rds"), overwrite = TRUE))
invisible(file.copy(primary$results_csv, file.path(out_dir, "deseq2_results.csv"), overwrite = TRUE))
message(glue("Primary contrast: {primary$comparison} / {primary$contrast}"))
