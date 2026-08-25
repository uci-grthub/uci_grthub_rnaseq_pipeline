#!/usr/bin/env Rscript
# Calculate TPM using tximport from Salmon quantification
# Usage: Rscript tximport_tpm.R output_dir gtf_file

suppressPackageStartupMessages({
  library(tximport)
  library(tidyverse)
  library(glue)
})

# Resolve the helper relative to this script so it works from any working dir
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(script_path), "tx2gene.R"))

args <- commandArgs(trailingOnly=TRUE)
default_out <- "output/salmon"
default_gtf <- "reference/genome.gtf"

# Use provided output dir if given, otherwise fall back to default
out_dir <- default_out
if(length(args) >= 1 && nzchar(args[1])) out_dir <- args[1]

# GTF file for transcript-to-gene mapping
gtf_file <- default_gtf
if(length(args) >= 2 && nzchar(args[2])) gtf_file <- args[2]

# Optional explicit output directory for the TPM files (defaults to
# dirname(out_dir)/tpm, e.g. when splitting runs per species)
tpm_dir_override <- NULL
if(length(args) >= 3 && nzchar(args[3])) tpm_dir_override <- args[3]

quant_files <- find_quant_files(out_dir)

message("Sample files:")
print(quant_files)

tx2gene <- build_tx2gene(gtf_file, quant_files[1])

# Import quantifications using tximport
message("Importing quantifications with tximport...")
txi <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE,
  countsFromAbundance = "no"  # Don't convert TPM to counts
)

# Second pass at transcript resolution. txOut=TRUE returns before
# summarizeToGene, so tx2gene/ignoreTxVersion are not applied and rownames stay
# the raw quant.sf Name field.
message("Importing quantifications at transcript level (txOut = TRUE)...")
txi_tx <- tximport(
  quant_files,
  type = "salmon",
  txOut = TRUE,
  countsFromAbundance = "no"
)

# Extract TPM (abundance) matrix
tpm_matrix <- txi$abundance

# Convert to data frame with gene IDs
tpm_df <- as.data.frame(tpm_matrix)
tpm_df$gene <- rownames(tpm_df)
tpm_df <- tpm_df[, c("gene", setdiff(colnames(tpm_df), "gene"))]

# Create output directory if it doesn't exist
out_tpm_dir <- if(!is.null(tpm_dir_override)) tpm_dir_override else file.path(dirname(out_dir), "tpm")
dir.create(out_tpm_dir, showWarnings = FALSE, recursive = TRUE)

# Export TPM to CSV
tpm_csv_path <- file.path(out_tpm_dir, "tpm_salmon.csv")
write.csv(tpm_df, file = tpm_csv_path, row.names = FALSE)
message(glue("TPM matrix exported to: {tpm_csv_path}"))

# Export TPM to RDS
tpm_rds_path <- file.path(out_tpm_dir, "tpm_salmon.rds")
saveRDS(tpm_matrix, file = tpm_rds_path)
message(glue("TPM matrix (RDS) saved to: {tpm_rds_path}"))

# Export full tximport object
txi_rds_path <- file.path(out_tpm_dir, "txi_salmon.rds")
saveRDS(txi, file = txi_rds_path)
message(glue("Full tximport object saved to: {txi_rds_path}"))

# --- Transcript-level matrices ---------------------------------------------
# GENCODE FASTA headers are pipe-delimited (ENSMUST00000193812.1|ENSMUSG...|...),
# so keep only the leading transcript ID. A plain index without pipes is
# unaffected.
tx_ids_full <- rownames(txi_tx$counts)
tx_ids <- sub("\\|.*", "", tx_ids_full)
rownames(txi_tx$abundance) <- tx_ids
rownames(txi_tx$counts) <- tx_ids
rownames(txi_tx$length) <- tx_ids

# Map each transcript back to its gene using the same version-stripped keys
# tximport uses for gene summarisation.
tx_gene_ids <- tx2gene$GENEID[match(sub("\\..*", "", tx_ids), tx2gene$TXNAME)]
# The dummy tx2gene fallback keys on the raw quant.sf Name, so the stripped IDs
# above miss. Fall back to the transcript ID rather than emitting an NA column.
tx_gene_ids[is.na(tx_gene_ids)] <- tx_ids[is.na(tx_gene_ids)]

build_tx_df <- function(mat){
  df <- as.data.frame(mat)
  df <- cbind(
    transcript = tx_ids,
    gene = tx_gene_ids,
    df,
    stringsAsFactors = FALSE
  )
  rownames(df) <- NULL
  df
}

tx_counts_csv_path <- file.path(out_tpm_dir, "transcript_counts.csv")
write.csv(build_tx_df(txi_tx$counts), file = tx_counts_csv_path, row.names = FALSE)
message(glue("Transcript count matrix exported to: {tx_counts_csv_path}"))

tx_tpm_csv_path <- file.path(out_tpm_dir, "transcript_tpm.csv")
write.csv(build_tx_df(txi_tx$abundance), file = tx_tpm_csv_path, row.names = FALSE)
message(glue("Transcript TPM matrix exported to: {tx_tpm_csv_path}"))

txi_tx_rds_path <- file.path(out_tpm_dir, "txi_transcript_salmon.rds")
saveRDS(txi_tx, file = txi_tx_rds_path)
message(glue("Transcript-level tximport object saved to: {txi_tx_rds_path}"))

message("TPM calculation complete!")
