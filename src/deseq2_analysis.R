#!/usr/bin/env Rscript 
# DESeq2 analysis script
# Usage: Rscript deseq2_analysis.R counts.txt metadata.csv output_dir

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(glue)
  library(AnnotationDbi)
  library(yaml)
  library(iSEEde)
  library(purrr)
})

args <- commandArgs(trailingOnly=TRUE)
default_counts <- "output/feature_count/all_samples_counts.txt"
## Defaults (allow running without providing all args)
default_meta <- "metadata/metadata.csv"
default_out <- "output/deseq2"

# Use provided counts file if given, otherwise fall back to default
counts_file <- default_counts
if(length(args) >= 1 && nzchar(args[1])) counts_file <- args[1]

# meta and out dir: use provided values if present, otherwise use defaults
meta_file <- default_meta
if(length(args) >= 2 && nzchar(args[2])) meta_file <- args[2]

out_dir <- default_out
if(length(args) >= 3 && nzchar(args[3])) out_dir <- args[3]
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

meta <- read.csv(meta_file)  |> 
janitor::clean_names() |>
mutate(index_pair = glue("{i7_barcode}-{i5_barcode_nova_seq_v1_5}"))  |> 
dplyr::right_join(
  tibble(index_pair = str_extract(colnames(count_matrix), "[ACGT]+-[ACGT]+")),
  by = c("index_pair")
) |>
dplyr::mutate(condition = group_name)  |> 
tibble::column_to_rownames("sample") |>
  identity()

colnames(count_matrix) <- rownames(meta)

print(rownames(meta))
print(colnames(count_matrix))
# Ensure sample names match
count_matrix <- count_matrix[, rownames(meta)]

dds <- DESeqDataSetFromMatrix(countData=count_matrix, 
  colData=meta, 
  design=~condition)

dds <- DESeq(dds)


# 6. Quality Control and Visualization
# ----------------------------------------------------------- #
# Optional: Transform data for visualization
vsd <- vst(dds, blind = FALSE) # Variance stabilizing transformation
rld <- rlog(dds, blind = FALSE) # Regularized log transformation


plot_var <- c("condition")  |> 
set_names()

pca_plots <- map(plot_var, ~{plotPCA(vsd, intgroup = .x) + labs(title = glue("PCA - {.x}"))})

pdf("results/pca_plots.pdf")
print(pca_plots)
dev.off()

comparisons_yaml <- file.path(dirname(meta_file), "comparisons.yaml")

# If a comparisons YAML exists, read it; otherwise create a default and write it.
if(file.exists(comparisons_yaml)){
  comparisons <- tryCatch({
    yaml::read_yaml(comparisons_yaml)
  }, error = function(e){
    warning(glue::glue("Failed to read {comparisons_yaml}: {e$message}. Falling back to defaults."))
    NULL
  })
  # validate structure
  if(is.null(comparisons) || !is.list(comparisons) || length(comparisons) == 0){
    warning(glue::glue("Invalid comparisons in {comparisons_yaml}; using defaults."))
    comparisons <- list(c("CHOW", "VV"), c("CHOW", "HFD"), c("VV", "HFD+VV"))
    try({ yaml::write_yaml(comparisons, comparisons_yaml) }, silent = TRUE)
  }
} else {
  comparisons <- list(c("CHOW", "VV"), c("CHOW", "HFD"), c("VV", "HFD+VV"))
  dir.create(dirname(meta_file), showWarnings = FALSE, recursive = TRUE)
  tryCatch({
    yaml::write_yaml(comparisons, comparisons_yaml)
    message(glue::glue("Wrote comparisons to {comparisons_yaml}"))
  }, error = function(e){
    warning(glue::glue("Failed to write comparisons YAML: {e$message}"))
  })
}

# Iterate over defined comparisons, extract results for each, and save RDS + CSV
safe_filename <- function(x){
  x <- stringr::str_replace_all(x, "\\+", "plus")
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9_-]", "_")
  tolower(x)
}

all_res_list <- list()
for(cmp in comparisons){
  group_a <- cmp[1]
  group_b <- cmp[2]

  message(glue::glue("Extracting results: {group_a} vs {group_b}"))

  # results() uses contrast = c("factorName","level1","level2")
  res_i <- tryCatch(
    results(dds, contrast = c("condition", group_a, group_b)),
    error = function(e){
      warning(glue::glue("Failed to get results for {group_a} vs {group_b}: {e$message}"))
      NULL
    }
  )

  if(is.null(res_i)) next
  # filter out rows where pvalue is NA
  res_df0 <- as.data.frame(res_i)
  if("pvalue" %in% colnames(res_df0)){
    res_df <- res_df0[!is.na(res_df0$pvalue), , drop = FALSE]
  } else {
    res_df <- res_df0
  }

  if(nrow(res_df) == 0){
    message(glue::glue("No rows with non-NA pvalue for {group_a} vs {group_b}; skipping."))
    next
  }

  nm <- glue::glue("{group_a}_vs_{group_b}") %>% safe_filename()
  rds_path <- file.path(out_dir, glue::glue("results_{nm}.rds"))
  csv_path <- file.path(out_dir, glue::glue("results_{nm}.csv"))

  # save filtered results: save RDS as DESeqResults subset and CSV as data.frame
  keep_idx <- rownames(res_df)
  res_i_filt <- tryCatch(res_i[keep_idx, ], error = function(e) res_i)
  saveRDS(res_i_filt, file = rds_path)

  # store for combined table (add contrast column)
  df_i <- res_df
  df_i$gene <- rownames(df_i)
  df_i$contrast <- glue::glue("{group_a}_vs_{group_b}")

  # map Ensembl IDs (trim version suffixes) to gene symbols using Bioconductor OrgDb
  gene_ids <- df_i$gene
  gene_ids_nover <- sub("\\.\\d+$", "", gene_ids)
  species <- NULL
  if(any(grepl('^ENSMUSG', gene_ids_nover))){
    species <- 'mouse'
  } else if(any(grepl('^ENSG', gene_ids_nover))){
    species <- 'human'
  }

  df_i$gene_symbol <- NA_character_
  if(!is.null(species)){
    if(!requireNamespace('AnnotationDbi', quietly=TRUE)){
      warning('AnnotationDbi not available; cannot map gene IDs to symbols')
    } else {
      OrgDb <- NULL
      if(species == 'mouse'){
        if(requireNamespace('org.Mm.eg.db', quietly=TRUE)){
          OrgDb <- get('org.Mm.eg.db', envir = asNamespace('org.Mm.eg.db'))
        } else {
          warning("Bioconductor package 'org.Mm.eg.db' not installed; install it to map mouse Ensembl IDs to symbols")
        }
      } else if(species == 'human'){
        if(requireNamespace('org.Hs.eg.db', quietly=TRUE)){
          OrgDb <- get('org.Hs.eg.db', envir = asNamespace('org.Hs.eg.db'))
        } else {
          warning("Bioconductor package 'org.Hs.eg.db' not installed; install it to map human Ensembl IDs to symbols")
        }
      }

      if(!is.null(OrgDb)){
        map_df <- tryCatch({
          AnnotationDbi::select(OrgDb,
                                keys = unique(gene_ids_nover),
                                keytype = 'ENSEMBL',
                                columns = c('SYMBOL'))
        }, error = function(e){
          warning(paste('AnnotationDbi::select failed:', e$message))
          NULL
        })

         if(!is.null(map_df) && nrow(map_df) > 0){
          map_df <- map_df[!duplicated(map_df$ENSEMBL), , drop = FALSE]
          map_df <- tibble::as_tibble(map_df)
          map_df <- dplyr::rename(map_df, gene_nover = ENSEMBL, gene_symbol = SYMBOL)
          df_i$gene_nover <- gene_ids_nover
          df_i <- dplyr::left_join(df_i, map_df, by = c('gene_nover'))

          # coalesce possible gene_symbol columns produced by joins (gene_symbol.x / gene_symbol.y)
          gs_cols <- c("gene_symbol.x", "gene_symbol.y", "gene_symbol")
          present_cols <- gs_cols[gs_cols %in% colnames(df_i)]
          if(length(present_cols) > 0){
            tmp <- NULL
            for(col in present_cols){
              vec <- df_i[[col]]
              if(is.null(tmp)){
                tmp <- vec
              } else {
                tmp <- dplyr::coalesce(tmp, vec)
              }
            }
            df_i$gene_symbol <- tmp
            df_i <- dplyr::select(df_i, -dplyr::any_of(c("gene_symbol.x", "gene_symbol.y")))
          }
        }
      }
    }
  }

  # remove temporary gene_nover column if present
  if("gene_nover" %in% colnames(df_i)) df_i$gene_nover <- NULL

  # ensure gene_symbol column exists
  if(!"gene_symbol" %in% colnames(df_i)) df_i$gene_symbol <- NA_character_

  # reorder columns: gene, gene_symbol, then everything else
  other_cols <- setdiff(colnames(df_i), c("gene", "gene_symbol"))
  df_i <- df_i[, c("gene", "gene_symbol", other_cols), drop = FALSE]

  # sort by increasing padj (NA last) if padj column exists
  if("padj" %in% colnames(df_i)){
    df_i <- df_i[order(df_i$padj, na.last = TRUE), , drop = FALSE]
  }

  write.csv(df_i, file = csv_path, row.names = FALSE)
  all_res_list[[nm]] <- df_i
}

# Save combined results if any
if(length(all_res_list) > 0){
  combined <- dplyr::bind_rows(all_res_list)
  if("padj" %in% colnames(combined)){
    combined <- combined[order(combined$padj, na.last = TRUE), , drop = FALSE]
  }
  write.csv(combined, file = file.path(out_dir, "results_all_contrasts.csv"), row.names = FALSE)
}

# Save the full DESeq dataset for downstream use
saveRDS(dds, file=file.path(out_dir, "dds.rds"))
