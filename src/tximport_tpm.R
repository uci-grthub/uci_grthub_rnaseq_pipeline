#!/usr/bin/env Rscript
# Calculate TPM using tximport from Salmon quantification
# Usage: Rscript tximport_tpm.R output_dir gtf_file

suppressPackageStartupMessages({
  library(tximport)
  library(tidyverse)
  library(glue)
})

args <- commandArgs(trailingOnly=TRUE)
default_out <- "output/salmon"
default_gtf <- "reference/genome.gtf"

# Use provided output dir if given, otherwise fall back to default
out_dir <- default_out
if(length(args) >= 1 && nzchar(args[1])) out_dir <- args[1]

# GTF file for transcript-to-gene mapping
gtf_file <- default_gtf
if(length(args) >= 2 && nzchar(args[2])) gtf_file <- args[2]

# Find all quant.sf files
quant_files <- list.files(
  path = out_dir,
  pattern = "quant.sf",
  recursive = TRUE,
  full.names = TRUE
)

if(length(quant_files) == 0){
  stop(glue("No quant.sf files found in {out_dir}"))
}

message(glue("Found {length(quant_files)} Salmon quantification files"))

# Extract sample names from file paths
# Assuming structure: output_dir/SAMPLE_salmon_quant/SAMPLE_quant.sf
sample_names <- sapply(quant_files, function(path){
  parts <- strsplit(path, "/")[[1]]
  # Get the directory name like "SAMPLE_salmon_quant"
  dir_name <- parts[length(parts) - 1]
  # Extract sample name
  sample_name <- gsub("_salmon_quant", "", dir_name)
  sample_name
}, USE.NAMES = FALSE)

names(quant_files) <- sample_names

message("Sample files:")
print(quant_files)

# Read transcript-to-gene mapping from GTF
# Simple approach: extract gene_id from GTF file
tx2gene <- NULL
if(file.exists(gtf_file)){
  message(glue("Reading transcript-to-gene mapping from {gtf_file}"))
  
  tx2gene <- tryCatch({
    gtf <- read.delim(gtf_file, header=FALSE, comment.char="#")
    
    # Extract attributes from GTF (column 9)
    attributes <- gtf$V9
    
    # Parse transcript_id and gene_id
    tx_ids <- stringr::str_extract(attributes, 'transcript_id "([^"]+)"')
    tx_ids <- stringr::str_remove(tx_ids, 'transcript_id "')
    tx_ids <- stringr::str_remove(tx_ids, '"')
    
    gene_ids <- stringr::str_extract(attributes, 'gene_id "([^"]+)"')
    gene_ids <- stringr::str_remove(gene_ids, 'gene_id "')
    gene_ids <- stringr::str_remove(gene_ids, '"')
    
    # Remove duplicates and create mapping
    mapping <- data.frame(
      TXNAME = tx_ids,
      GENEID = gene_ids,
      stringsAsFactors = FALSE
    )
    mapping <- mapping[!duplicated(mapping$TXNAME), ]
    mapping
  }, error = function(e){
    warning(glue("Failed to parse GTF: {e$message}. Using transcript IDs as gene IDs."))
    NULL
  })
} else {
  warning(glue("GTF file not found at {gtf_file}"))
}

# If tx2gene is NULL, create a dummy mapping using transcript IDs as gene IDs
if(is.null(tx2gene)){
  message("Creating transcript-to-gene mapping from quantification files")
  
  # Read one quant.sf file to get transcript names
  sample_quant <- read.delim(quant_files[1])
  tx2gene <- data.frame(
    TXNAME = sample_quant$Name,
    GENEID = sample_quant$Name,  # Use transcript ID as gene ID
    stringsAsFactors = FALSE
  )
}

message(glue("Transcript-to-gene mapping has {nrow(tx2gene)} entries"))

# Import quantifications using tximport
message("Importing quantifications with tximport...")
txi <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE,
  countsFromAbundance = "no"  # Don't convert TPM to counts
)

# Extract TPM (abundance) matrix
tpm_matrix <- txi$abundance

# Convert to data frame with gene IDs
tpm_df <- as.data.frame(tpm_matrix)
tpm_df$gene <- rownames(tpm_df)
tpm_df <- tpm_df[, c("gene", setdiff(colnames(tpm_df), "gene"))]

# Create output directory if it doesn't exist
out_tpm_dir <- file.path(dirname(out_dir), "tpm")
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

message("TPM calculation complete!")
