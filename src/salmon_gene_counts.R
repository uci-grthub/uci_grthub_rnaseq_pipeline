#!/usr/bin/env Rscript
# Build a gene-level count matrix from Salmon quantification via tximport.
#
# This is the Salmon-based counterpart to count_matrix.R, which builds the same
# deliverables from featureCounts/HISAT2. The two are independent branches and
# their numbers are not expected to agree: Salmon distributes multi-mapping
# reads across transcripts by EM, while featureCounts discards them.
#
# Counts use countsFromAbundance="lengthScaledTPM". Raw tximport counts are not
# valid input to a plain count matrix -- they carry no correction for
# differential isoform usage changing a gene's average transcript length across
# samples. lengthScaledTPM bakes that correction into the counts themselves, so
# the matrix can be handed to DESeq2 without the separate length offset that
# DESeqDataSetFromTximport would otherwise supply.
#
# Usage: Rscript salmon_gene_counts.R <salmon_dir> <gtf> <metadata_csv> <out_dir>

suppressPackageStartupMessages({
  library(tximport)
  library(tidyverse)
  library(glue)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(script_path), "tx2gene.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript salmon_gene_counts.R <salmon_dir> <gtf> <metadata_csv> <out_dir>")
}
salmon_dir <- args[1]
gtf_file <- args[2]
metadata_path <- args[3]
out_dir <- args[4]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

quant_files <- find_quant_files(salmon_dir)
tx2gene <- build_tx2gene(gtf_file, quant_files[1])

message("Importing gene-level counts with tximport (lengthScaledTPM)...")
txi <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE,
  countsFromAbundance = "lengthScaledTPM"
)

## ---- Sample labelling ------------------------------------------------------
# Mirrors count_matrix.R so the two branches produce comparable column names.

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"sample" %in% colnames(metadata)) {
  stop(glue("Metadata {metadata_path} has no 'sample' column"))
}

fastq_names <- colnames(txi$counts)

label_for <- function(fastq_name) {
  row <- metadata[metadata$sample == fastq_name, , drop = FALSE]
  if (nrow(row) == 1 && "sample_id" %in% colnames(row) && nzchar(as.character(row$sample_id[1]))) {
    as.character(row$sample_id[1])
  } else {
    fastq_name
  }
}
sample_labels <- vapply(fastq_names, label_for, character(1), USE.NAMES = FALSE)

if (anyDuplicated(sample_labels)) {
  warning("Duplicate sample labels after metadata mapping; falling back to FASTQ names")
  sample_labels <- fastq_names
}

# lengthScaledTPM values are non-integer expectations; DESeq2 and the CPM below
# both want whole counts.
counts <- round(txi$counts)
mode(counts) <- "integer"
colnames(counts) <- sample_labels

gene_lengths <- txi$length
colnames(gene_lengths) <- sample_labels

# Order columns to match metadata order so downstream group blocks stay contiguous
metadata_order <- match(metadata$sample, fastq_names)
metadata_order <- metadata_order[!is.na(metadata_order)]
if (length(metadata_order) == ncol(counts)) {
  counts <- counts[, metadata_order, drop = FALSE]
  gene_lengths <- gene_lengths[, metadata_order, drop = FALSE]
}

## ---- Deliverables ----------------------------------------------------------

counts_df <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
counts_csv <- file.path(out_dir, "gene_counts.csv")
write.csv(counts_df, counts_csv, row.names = FALSE)
message(glue("Salmon gene-level counts written to {counts_csv} ({nrow(counts)} genes x {ncol(counts)} samples)"))

saveRDS(counts, file.path(out_dir, "gene_counts.rds"))
saveRDS(txi, file.path(out_dir, "txi_gene_lengthscaledtpm.rds"))

# Unlike featureCounts, Salmon has no single gene length: the effective length
# is per sample. Report the mean alongside the spread so a gene whose isoform
# usage shifts between samples is visible.
annotation_df <- data.frame(
  gene_id = rownames(gene_lengths),
  mean_effective_length = round(rowMeans(gene_lengths), 2),
  min_effective_length = round(apply(gene_lengths, 1, min), 2),
  max_effective_length = round(apply(gene_lengths, 1, max), 2),
  check.names = FALSE
)
write.csv(annotation_df, file.path(out_dir, "gene_annotation.csv"), row.names = FALSE)

# CPM for quick inspection; not a substitute for the DESeq2/TPM normalisations
lib_sizes <- colSums(counts)
cpm <- sweep(counts, 2, lib_sizes, "/") * 1e6
cpm_df <- data.frame(gene_id = rownames(cpm), round(cpm, 4), check.names = FALSE)
write.csv(cpm_df, file.path(out_dir, "gene_counts_cpm.csv"), row.names = FALSE)

## ---- Per-sample metrics ----------------------------------------------------
# Salmon's mapping rate is the analogue of the featureCounts assignment rate.

metrics <- data.frame(
  sample = colnames(counts),
  library_size = as.numeric(lib_sizes),
  genes_detected = as.integer(colSums(counts > 0)),
  check.names = FALSE
)

read_salmon_meta <- function(quant_file) {
  meta_path <- file.path(dirname(quant_file), "aux_info", "meta_info.json")
  if (!file.exists(meta_path)) {
    return(data.frame(num_processed = NA_real_, num_mapped = NA_real_, percent_mapped = NA_real_))
  }
  info <- jsonlite::fromJSON(meta_path)
  data.frame(
    num_processed = as.numeric(info$num_processed),
    num_mapped = as.numeric(info$num_mapped),
    percent_mapped = round(as.numeric(info$percent_mapped), 2)
  )
}

if (requireNamespace("jsonlite", quietly = TRUE)) {
  salmon_meta <- do.call(rbind, lapply(quant_files, read_salmon_meta))
  salmon_meta$sample <- sample_labels
  metrics <- merge(metrics, salmon_meta, by = "sample", all.x = TRUE, sort = FALSE)
} else {
  warning("jsonlite not available; skipping Salmon mapping-rate metrics")
}

write.csv(metrics, file.path(out_dir, "count_matrix_metrics.csv"), row.names = FALSE)

message("Salmon count matrix export complete.")
