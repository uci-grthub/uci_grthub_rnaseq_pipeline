#!/usr/bin/env Rscript
# Build a gene-level count matrix from Salmon quantification via tximeta.
#
# This is the provenance-aware alternative to salmon_gene_counts.R. tximeta
# reads the index checksum Salmon recorded in aux_info/meta_info.json and
# matches it against its own table of known transcriptomes, so the annotation
# is identified from the quantification itself rather than from a GTF path
# passed in on the command line. The mouse index here matches GENCODE vM24
# (GRCm38) exactly, so no makeLinkedTxome() step is needed.
#
# What this buys over the tximport route:
#   - rowRanges and genome build travel with the counts
#   - no hand-rolled GTF parsing, so no tx2gene and no version-stripping
#   - the saved SummarizedExperiment can go straight into DESeqDataSet(),
#     which applies the average-transcript-length offsets natively
#
# Usage: Rscript salmon_gene_counts_tximeta.R <salmon_dir> <metadata_csv> <out_dir> [linked_txome_json]

suppressPackageStartupMessages({
  library(tximeta)
  library(SummarizedExperiment)
  library(tidyverse)
  library(glue)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(script_path), "tx2gene.R")) # for find_quant_files()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript salmon_gene_counts_tximeta.R <salmon_dir> <metadata_csv> <out_dir> [linked_txome_json]")
}
salmon_dir <- args[1]
metadata_path <- args[2]
out_dir <- args[3]
linked_txome_json <- if (length(args) >= 4 && nzchar(args[4])) args[4] else NULL

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# An explicit linkedTxome pins the annotation to a local GTF, for the case
# where the index is not in tximeta's table or the machine has no network.
if (!is.null(linked_txome_json) && file.exists(linked_txome_json)) {
  message(glue("Loading linkedTxome from {linked_txome_json}"))
  loadLinkedTxome(linked_txome_json)
}

quant_files <- find_quant_files(salmon_dir)

## ---- Sample labelling ------------------------------------------------------
# Mirrors salmon_gene_counts.R / count_matrix.R so all branches agree on names.

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"sample" %in% colnames(metadata)) {
  stop(glue("Metadata {metadata_path} has no 'sample' column"))
}

fastq_names <- names(quant_files)

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

# Order to match metadata so downstream group blocks stay contiguous
metadata_order <- match(metadata$sample, fastq_names)
metadata_order <- metadata_order[!is.na(metadata_order)]
if (length(metadata_order) == length(quant_files)) {
  quant_files <- quant_files[metadata_order]
  fastq_names <- fastq_names[metadata_order]
  sample_labels <- sample_labels[metadata_order]
}

coldata <- data.frame(
  files = as.character(quant_files),
  names = sample_labels,
  fastq_name = fastq_names,
  stringsAsFactors = FALSE
)

## ---- Import ----------------------------------------------------------------

message("Importing with tximeta (transcript level)...")
se <- tximeta(coldata)

meta <- metadata(se)
txome <- meta$txomeInfo
if (!is.null(txome)) {
  message(glue(
    "Detected transcriptome: {txome$source} {txome$organism} ",
    "release {txome$release} genome {txome$genome}"
  ))
} else {
  warning("tximeta could not identify the transcriptome; rowRanges will be absent")
}

message("Summarising to gene level (lengthScaledTPM)...")
gse <- summarizeToGene(se, countsFromAbundance = "lengthScaledTPM")

# The SummarizedExperiment is the point of this script: it carries rowRanges
# and the genome build, and DESeqDataSet() consumes it directly.
saveRDS(gse, file.path(out_dir, "gse_tximeta.rds"))
saveRDS(se, file.path(out_dir, "se_tximeta_transcript.rds"))

## ---- Deliverables ----------------------------------------------------------
# Same shapes as salmon_gene_counts.R so the branches stay comparable.

counts <- round(assay(gse, "counts"))
mode(counts) <- "integer"

counts_df <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
counts_csv <- file.path(out_dir, "gene_counts.csv")
write.csv(counts_df, counts_csv, row.names = FALSE)
message(glue("tximeta gene-level counts written to {counts_csv} ({nrow(counts)} genes x {ncol(counts)} samples)"))

saveRDS(counts, file.path(out_dir, "gene_counts.rds"))

gene_lengths <- assay(gse, "length")
annotation_df <- data.frame(
  gene_id = rownames(gene_lengths),
  mean_effective_length = round(rowMeans(gene_lengths), 2),
  min_effective_length = round(apply(gene_lengths, 1, min), 2),
  max_effective_length = round(apply(gene_lengths, 1, max), 2),
  check.names = FALSE
)

# Unlike the tximport route, tximeta knows where each gene lives; carry that
# through so the annotation table is usable without a second lookup.
rr <- rowRanges(gse)
if (!is.null(rr) && length(rr) == nrow(gene_lengths)) {
  annotation_df$seqnames <- as.character(GenomicRanges::seqnames(rr))
  annotation_df$start <- GenomicRanges::start(rr)
  annotation_df$end <- GenomicRanges::end(rr)
  annotation_df$strand <- as.character(GenomicRanges::strand(rr))
}
write.csv(annotation_df, file.path(out_dir, "gene_annotation.csv"), row.names = FALSE)

lib_sizes <- colSums(counts)
cpm <- sweep(counts, 2, lib_sizes, "/") * 1e6
cpm_df <- data.frame(gene_id = rownames(cpm), round(cpm, 4), check.names = FALSE)
write.csv(cpm_df, file.path(out_dir, "gene_counts_cpm.csv"), row.names = FALSE)

## ---- Per-sample metrics ----------------------------------------------------

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

salmon_meta <- do.call(rbind, lapply(coldata$files, read_salmon_meta))
salmon_meta$sample <- sample_labels
metrics <- merge(metrics, salmon_meta, by = "sample", all.x = TRUE, sort = FALSE)
write.csv(metrics, file.path(out_dir, "count_matrix_metrics.csv"), row.names = FALSE)

# Record what tximeta resolved, so a run's provenance is inspectable without
# opening the RDS.
if (!is.null(txome)) {
  prov <- data.frame(
    source = txome$source,
    organism = txome$organism,
    release = as.character(txome$release),
    genome = txome$genome,
    index_sha256 = paste(txome$sha256, collapse = ","),
    tximeta_version = as.character(packageVersion("tximeta")),
    stringsAsFactors = FALSE
  )
  write.csv(prov, file.path(out_dir, "txome_provenance.csv"), row.names = FALSE)
}

message("tximeta count matrix export complete.")
