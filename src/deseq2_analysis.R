#!/usr/bin/env Rscript 
# DESeq2 analysis script
# Usage: Rscript deseq2_analysis.R counts.txt metadata.csv output_dir

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(glue)
})

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 3) stop("Usage: Rscript deseq2_analysis.R counts.txt metadata.csv output_dir condition group_a group_b")

counts_file <- args[1]
meta_file <- args[2]
out_dir <- args[3]
condition <- args[4]  # Condition column in metadata
group_a <- args[5]  # First group for comparison
group_b <- args[6]  # Second group for comparison

dir.create(out_dir, showWarnings=FALSE)

counts <- read.table(counts_file, header=TRUE, row.names=1)

# Extract the count matrix (assuming first 5 columns are metadata)
count_matrix <- counts[,6:ncol(counts)]

# format column name to match sample name
format_sample_id <- function(colname){
    colname |>
    str_remove("results.hisat2_alignment\\.") %>%
  # Remove the suffix
  str_remove("_align_sorted.bam") %>%
  # Replace all periods with hyphens
  str_replace_all("\\.", "-")
}

colnames(count_matrix) <- sapply(colnames(count_matrix), format_sample_id)

meta <- read.csv(meta_file, row.names=1)

print(rownames(meta))
print(colnames(count_matrix))
# Ensure sample names match
count_matrix <- count_matrix[, rownames(meta)]

dds <- DESeqDataSetFromMatrix(countData=count_matrix, 
  colData=meta, 
  design=~condition)

dds <- DESeq(dds)

res <- results(dds)

dds <- embedContrastResults(res, dds, name = glue("{condition}: {group_a} vs {group_b}")
contrastResults(dds)

saveRDS(dds, file=file.path(out_dir, "dds.rds"))

write.csv(as.data.frame(res), file=file.path(out_dir, "deseq2_results.csv"))
