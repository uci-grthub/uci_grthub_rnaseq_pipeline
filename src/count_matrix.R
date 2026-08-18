#!/usr/bin/env Rscript
# Build a clean gene-level count matrix from the featureCounts output.
#
# featureCounts writes a `#` command line, six annotation columns, and column
# names that are full BAM paths. None of that is usable as a deliverable, so
# this splits it into a plain counts matrix, a gene-length table, and a CPM
# matrix, relabelling columns with the metadata sample_id where one exists.
#
# Usage: Rscript count_matrix.R <counts_txt> <metadata_csv> <out_dir>

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript count_matrix.R <counts_txt> <metadata_csv> <out_dir>")
}
counts_path <- args[1]
metadata_path <- args[2]
out_dir <- args[3]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# comment.char="#" drops the featureCounts command line that precedes the header
fc <- read.delim(counts_path, comment.char = "#", check.names = FALSE)

annotation_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
missing_cols <- setdiff(annotation_cols, colnames(fc))
if (length(missing_cols) > 0) {
  stop(glue("featureCounts table is missing columns: {paste(missing_cols, collapse=', ')}"))
}

count_cols <- setdiff(colnames(fc), annotation_cols)
if (length(count_cols) == 0) stop(glue("No count columns found in {counts_path}"))

# Column names are BAM paths: output/hisat2_alignment/<sample>_align_sorted_markdup.bam
fastq_names <- basename(count_cols) |>
  sub("_align_sorted_markdup\\.bam$", "", x = _) |>
  sub("\\.bam$", "", x = _)

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"sample" %in% colnames(metadata)) {
  stop(glue("Metadata {metadata_path} has no 'sample' column"))
}

# Prefer the short experimental label (GAC1) over the sequencing-core filename
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

counts <- as.matrix(fc[, count_cols, drop = FALSE])
mode(counts) <- "integer"
rownames(counts) <- fc$Geneid
colnames(counts) <- sample_labels

# Order columns to match metadata order so downstream group blocks stay contiguous
metadata_order <- match(metadata$sample, fastq_names)
metadata_order <- metadata_order[!is.na(metadata_order)]
if (length(metadata_order) == ncol(counts)) {
  counts <- counts[, metadata_order, drop = FALSE]
}

counts_df <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
counts_csv <- file.path(out_dir, "gene_counts.csv")
write.csv(counts_df, counts_csv, row.names = FALSE)
message(glue("Raw gene-level counts written to {counts_csv} ({nrow(counts)} genes x {ncol(counts)} samples)"))

lengths_df <- data.frame(
  gene_id = fc$Geneid,
  chr = fc$Chr,
  start = fc$Start,
  end = fc$End,
  strand = fc$Strand,
  length = fc$Length,
  check.names = FALSE
)
write.csv(lengths_df, file.path(out_dir, "gene_annotation.csv"), row.names = FALSE)

# CPM for quick inspection; not a substitute for the DESeq2/TPM normalisations
lib_sizes <- colSums(counts)
cpm <- sweep(counts, 2, lib_sizes, "/") * 1e6
cpm_df <- data.frame(gene_id = rownames(cpm), round(cpm, 4), check.names = FALSE)
write.csv(cpm_df, file.path(out_dir, "gene_counts_cpm.csv"), row.names = FALSE)

saveRDS(counts, file.path(out_dir, "gene_counts.rds"))

# Per-sample assignment metrics, taken from the featureCounts .summary sidecar
summary_path <- paste0(counts_path, ".summary")
metrics <- data.frame(
  sample = colnames(counts),
  library_size = as.numeric(lib_sizes),
  genes_detected = as.integer(colSums(counts > 0)),
  check.names = FALSE
)
if (file.exists(summary_path)) {
  fc_summary <- read.delim(summary_path, check.names = FALSE)
  summary_labels <- basename(colnames(fc_summary)[-1]) |>
    sub("_align_sorted_markdup\\.bam$", "", x = _)
  assigned <- as.numeric(fc_summary[fc_summary$Status == "Assigned", -1])
  total <- colSums(fc_summary[, -1, drop = FALSE])
  assign_df <- data.frame(
    fastq_name = summary_labels,
    assigned_reads = assigned,
    total_reads_counted = as.numeric(total),
    percent_assigned = round(100 * assigned / as.numeric(total), 2),
    check.names = FALSE
  )
  assign_df$sample <- vapply(assign_df$fastq_name, label_for, character(1), USE.NAMES = FALSE)
  metrics <- merge(metrics, assign_df, by = "sample", all.x = TRUE, sort = FALSE)
}
write.csv(metrics, file.path(out_dir, "count_matrix_metrics.csv"), row.names = FALSE)

message("Count matrix export complete.")
