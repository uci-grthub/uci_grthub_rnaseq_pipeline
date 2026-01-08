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
  library(ggrepel)
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
mutate(index_pair = glue("{i7_barcode_sequence}-{i5_barcode_sequence}"))  |> 
dplyr::right_join(
  tibble(index_pair = str_extract(colnames(count_matrix), "[ACGT]+-[ACGT]+")),
  by = c("index_pair")
) |>
tidyr::separate(sample_name, into = c("experiment", "age", "sex", "condition", "replicate"), sep = "_", remove = FALSE) |>
tibble::column_to_rownames("sample_name") |>
  identity()

# Ensure factors and set reference levels for interpretable contrasts
meta$condition <- factor(meta$condition)
if("C" %in% levels(meta$condition)){
  meta$condition <- relevel(meta$condition, ref = "C")
}

meta$age <- factor(meta$age)
if("10" %in% levels(meta$age)){
  meta$age <- relevel(meta$age, ref = "10")
}

meta$sex <- factor(meta$sex)

colnames(count_matrix) <- rownames(meta)

count_matrix  |> 
tibble::rownames_to_column("ensgene") |>
  readr::write_csv(file.path(out_dir, "all_sample_counts.csv"))

print(rownames(meta))
print(colnames(count_matrix))
# Ensure sample names match
count_matrix <- count_matrix[, rownames(meta)]

dds <- DESeqDataSetFromMatrix(countData=count_matrix, 
  colData=meta, 
  design=~condition+age+condition:age)
  
## Run combined DESeq2 once for cross-sex QC PCA (contrasts are sex-specific below)
dds_all <- DESeq(dds)
vsd_all <- vst(dds_all, blind = FALSE)
plot_var_all <- c("condition", "age", "sex") |> set_names()
pca_plots_all <- map(plot_var_all, ~{
  p <- plotPCA(vsd_all, intgroup = .x, returnData = TRUE)
  percentVar <- round(100 * attr(p, "percentVar"))
  ggplot(p, aes(PC1, PC2, color = !!sym(.x))) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = name), size = 3, max.overlaps = 20) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    labs(title = glue("PCA - all - {.x}")) +
    theme_bw()
})
dir.create("results", showWarnings = FALSE, recursive = TRUE)
pdf(file.path("results", "pca_plots_all.pdf"), width = 12, height = 6)
print(pca_plots_all)
dev.off()

## Running DESeq2 for downstream contrasts will be done separately per sex below


## PCA plots will also be generated within each sex-specific analysis below

comparisons_yaml <- file.path(dirname(meta_file), "comparisons.yaml")

# If a comparisons YAML exists, read it; otherwise create a default and write it.
if(file.exists(comparisons_yaml)){
  comparisons_config <- tryCatch({
    yaml::read_yaml(comparisons_yaml)
  }, error = function(e){
    warning(glue::glue("Failed to read {comparisons_yaml}: {e$message}. Falling back to defaults."))
    NULL
  })
  # validate structure
  if(is.null(comparisons_config) || !is.list(comparisons_config) || length(comparisons_config) == 0){
    warning(glue::glue("Invalid comparisons in {comparisons_yaml}; using defaults."))
    comparisons_config <- list(condition = list(c("R", "C")))
    try({ yaml::write_yaml(comparisons_config, comparisons_yaml) }, silent = TRUE)
  }
} else {
  comparisons_config <- list(condition = list(c("R", "C")))
  dir.create(dirname(meta_file), showWarnings = FALSE, recursive = TRUE)
  tryCatch({
    yaml::write_yaml(comparisons_config, comparisons_yaml)
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

match_level <- function(x, levs){
  idx <- match(tolower(x), tolower(levs))
  if(is.na(idx)) x else levs[idx]
}

find_interaction <- function(cond_level, age_level, rnames){
  cand <- rnames[grepl("condition", rnames, ignore.case = TRUE) & grepl("\\.age", rnames, ignore.case = TRUE)]
  pick <- cand[grepl(cond_level, cand, ignore.case = TRUE) & grepl(age_level, cand, ignore.case = TRUE)]
  if(length(pick) > 0) pick[1] else NULL
}

## Run analysis separately for each sex
sex_levels <- sort(unique(as.character(meta$sex)))
for(sx in sex_levels){
  message(glue::glue("Running DESeq2 for sex = {sx}"))
  out_dir_sex <- file.path(out_dir, glue::glue("sex_{sx}"))
  dir.create(out_dir_sex, showWarnings = FALSE, recursive = TRUE)

  # subset DESeqDataSet and run DESeq2
  dds_sex <- dds[, colData(dds)$sex == sx]
  dds_sex <- DESeq(dds_sex)

  # PCA per sex (condition and age vary within sex)
  vsd_sex <- vst(dds_sex, blind = FALSE)
  plot_var <- c("condition", "age") |> set_names()
  pca_plots <- map(plot_var, ~{
    p <- plotPCA(vsd_sex, intgroup = .x, returnData = TRUE)
    percentVar <- round(100 * attr(p, "percentVar"))
    ggplot(p, aes(PC1, PC2, color = !!sym(.x))) +
      geom_point(size = 3) +
      geom_text_repel(aes(label = name), size = 3, max.overlaps = 20) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      labs(title = glue("PCA - {sx} - {.x}")) +
      theme_bw()
  })
  pdf(file.path("results", glue::glue("pca_plots_{sx}.pdf")))
  print(pca_plots)
  dev.off()

  all_res_list <- list()

  # Process condition comparisons stratified by age (compare R vs C within each age)
  if(!is.null(comparisons_config$condition)){
    age_levels <- levels(colData(dds_sex)$age)
    for(age_level in age_levels){
      for(cmp in comparisons_config$condition){
        group_a <- cmp[1]
        group_b <- cmp[2]

        message(glue::glue("Condition comparison: {group_a} vs {group_b} at age={age_level} (sex={sx})"))

        res_i <- tryCatch({
          # Subset to specific age, refit, and compare condition within that age
          dds_age <- dds_sex[, colData(dds_sex)$age == age_level]
          if(ncol(dds_age) >= 2 && nlevels(droplevels(colData(dds_age)$condition)) >= 2){
            colData(dds_age)$condition <- droplevels(colData(dds_age)$condition)
            if("C" %in% levels(colData(dds_age)$condition)){
              colData(dds_age)$condition <- relevel(colData(dds_age)$condition, ref = "C")
            }
            # Use simplified design since age is constant within this subset
            design(dds_age) <- ~condition
            dds_age <- DESeq(dds_age)
            results(dds_age, contrast = c("condition", group_a, group_b))
          } else {
            warning(glue::glue("Skipping {group_a} vs {group_b} at age={age_level} (sex={sx}): insufficient samples or condition levels."))
            NULL
          }
        }, error = function(e){
          warning(glue::glue("Failed to get results for {group_a} vs {group_b} at age={age_level} (sex={sx}): {e$message}"))
          NULL
        })

        if(!is.null(res_i)){
          nm <- glue::glue("sex_{sx}_age_{age_level}_condition_{group_a}_vs_{group_b}") %>% safe_filename()
          all_res_list[[nm]] <- list(results = res_i, name = nm,
                                       contrast_type = "condition_within_age",
                                       age = age_level, group_a = group_a, group_b = group_b)
        }
      }
    }
  }

  # Process interaction effect: condition:age (tests if condition effect varies across ages)
  if(!is.null(comparisons_config$condition) && !is.null(comparisons_config$age)){
    message(glue::glue("Interaction effect - condition:age (sex={sx})"))
    
    # Get available interaction term names from the model
    available_names <- resultsNames(dds_sex)
    interaction_names <- available_names[grepl("condition.*age", available_names, ignore.case = TRUE)]
    
    if(length(interaction_names) > 0){
      for(int_name in interaction_names){
        res_i <- tryCatch(
          results(dds_sex, name = int_name),
          error = function(e){
            warning(glue::glue("Failed to get interaction results for {int_name} (sex={sx}): {e$message}"))
            NULL
          }
        )
        
        if(!is.null(res_i)){
          nm <- glue::glue("sex_{sx}_interaction_{int_name}") %>% safe_filename()
          all_res_list[[nm]] <- list(results = res_i, name = nm,
                                       contrast_type = "interaction",
                                       description = int_name)
        }
      }
    } else {
      message(glue::glue("No interaction terms found in resultsNames for sex={sx}"))
    }
  }

  # Process all collected results for this sex
  final_res_list <- list()
  for(nm in names(all_res_list)){
    res_item <- all_res_list[[nm]]
    res_i <- res_item$results

    if(is.null(res_i)) next

    # filter out rows where pvalue is NA
    res_df0 <- as.data.frame(res_i)
    if("pvalue" %in% colnames(res_df0)){
      res_df <- res_df0[!is.na(res_df0$pvalue), , drop = FALSE]
    } else {
      res_df <- res_df0
    }

    if(nrow(res_df) == 0){
      message(glue::glue("No rows with non-NA pvalue for {nm}; skipping."))
      next
    }

    rds_path <- file.path(out_dir_sex, glue::glue("results_{nm}.rds"))
    csv_path <- file.path(out_dir_sex, glue::glue("results_{nm}.csv"))

    # save filtered results: save RDS as DESeqResults subset and CSV as data.frame
    keep_idx <- rownames(res_df)
    res_i_filt <- tryCatch(res_i[keep_idx, ], error = function(e) res_i)
    saveRDS(res_i_filt, file = rds_path)

    # store for combined table (add contrast column)
    df_i <- res_df
    df_i$gene <- rownames(df_i)
    df_i$contrast <- nm
    df_i$sex <- sx

  # map gene IDs to gene symbols using Bioconductor OrgDb
  gene_ids <- df_i$gene
  gene_ids_nover <- sub("\\.\\d+$", "", gene_ids)
  species <- NULL
  keytype <- 'ENSEMBL'
  
  if(any(grepl('^ENSMUSG', gene_ids_nover))){
    species <- 'mouse'
    keytype <- 'ENSEMBL'
  } else if(any(grepl('^ENSG', gene_ids_nover))){
    species <- 'human'
    keytype <- 'ENSEMBL'
  } else if(any(grepl('^FBgn', gene_ids_nover))){
    species <- 'drosophila'
    keytype <- 'FLYBASE'
  } else if(any(grepl('^ENSMUSDM', gene_ids_nover))){
    # Alternative pattern for Drosophila Ensembl IDs
    species <- 'drosophila'
    keytype <- 'ENSEMBL'
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
      } else if(species == 'drosophila'){
        if(requireNamespace('org.Dm.eg.db', quietly=TRUE)){
          OrgDb <- get('org.Dm.eg.db', envir = asNamespace('org.Dm.eg.db'))
        } else {
          warning("Bioconductor package 'org.Dm.eg.db' not installed; install it to map Drosophila Ensembl IDs to symbols")
        }
      }

      if(!is.null(OrgDb)){
        map_df <- tryCatch({
          AnnotationDbi::select(OrgDb,
                                keys = unique(gene_ids_nover),
                                keytype = keytype,
                                columns = c('SYMBOL'))
        }, error = function(e){
          warning(paste('AnnotationDbi::select failed:', e$message))
          NULL
        })

         if(!is.null(map_df) && nrow(map_df) > 0){
          map_df <- map_df[!duplicated(map_df[[keytype]]), , drop = FALSE]
          map_df <- tibble::as_tibble(map_df)
          
          # Rename the key column to gene_nover for consistent joining
          key_col_name <- keytype
          if(keytype == 'FLYBASE'){
            map_df <- dplyr::rename(map_df, gene_nover = FLYBASE, gene_symbol = SYMBOL)
          } else if(keytype == 'ENSEMBL'){
            map_df <- dplyr::rename(map_df, gene_nover = ENSEMBL, gene_symbol = SYMBOL)
          }
          
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
    final_res_list[[nm]] <- df_i
  }

  # Save combined results for this sex
  if(length(final_res_list) > 0){
    combined <- dplyr::bind_rows(final_res_list)
    if("padj" %in% colnames(combined)){
      combined <- combined[order(combined$padj, na.last = TRUE), , drop = FALSE]
    }
    write.csv(combined, file = file.path(out_dir_sex, "results_all_contrasts.csv"), row.names = FALSE)
  }

  # Save the DESeq dataset for this sex
  saveRDS(dds_sex, file=file.path(out_dir_sex, "dds.rds"))
}
