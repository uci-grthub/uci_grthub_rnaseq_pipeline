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
  library(tximport)
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

# Extract gene lengths (typically in column 5 from featureCounts)
gene_lengths <- counts[,5]

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
dplyr::mutate(condition = condition)  |> 
tibble::column_to_rownames("sample_id") |>
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

# Ensure defaults for condition and age contrasts if not present in YAML
cond_levels_all <- levels(meta$condition)
age_levels_all <- levels(meta$age)

if(is.null(comparisons_config$condition) || length(comparisons_config$condition) == 0){
  if(length(cond_levels_all) >= 2){
    comparisons_config$condition <- list(c(cond_levels_all[2], cond_levels_all[1]))
  } else {
    comparisons_config$condition <- list()
  }
}

if(is.null(comparisons_config$age) || length(comparisons_config$age) == 0){
  if(length(age_levels_all) >= 2){
    comparisons_config$age <- list(c(age_levels_all[1], age_levels_all[2]))
  } else {
    comparisons_config$age <- list()
  }
}

safe_filename <- function(x){
  x <- stringr::str_replace_all(x, "\\+", "plus")
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9_-]", "_")
  tolower(x)
}

write_normalized_counts <- function(dds_obj, out_dir_path, label){
  counts_normalized <- counts(dds_obj, normalized = TRUE)
  write.csv(
    counts_normalized,
    file.path(out_dir_path, glue::glue("{label}_normalized_counts_wo_transformation.csv"))
  )
}

match_level <- function(x, levs){
  idx <- match(tolower(x), tolower(levs))
  if(is.na(idx)) x else levs[idx]
}

comparison_manifest_path <- file.path(getwd(), "deseq2_comparisons.csv")
comparison_rows <- list()

append_comparison_row <- function(sex, contrast_name, contrast_type, comparison_text,
                                  subset_factor = NA_character_, subset_level = NA_character_,
                                  level_a = NA_character_, level_b = NA_character_,
                                  interaction_term = NA_character_, design_formula = NA_character_){
  comparison_rows[[length(comparison_rows) + 1]] <<- tibble::tibble(
    sex = as.character(sex),
    contrast_name = as.character(contrast_name),
    comparison_type = as.character(contrast_type),
    subset_factor = as.character(subset_factor),
    subset_level = as.character(subset_level),
    level_a = as.character(level_a),
    level_b = as.character(level_b),
    interaction_term = as.character(interaction_term),
    comparison_text = as.character(comparison_text),
    design_formula = as.character(design_formula),
    results_csv = file.path(out_dir, glue::glue("sex_{sex}"), glue::glue("results_{contrast_name}.csv"))
  )
}

full_design_formula <- deparse(design(dds))

build_comparison_manifest <- function(){
  sex_values <- sort(unique(as.character(meta$sex)))
  cond_levels <- levels(meta$condition)
  age_levels <- levels(meta$age)

  for(sx in sex_values){
    if(!is.null(comparisons_config$condition)){
      for(cmp in comparisons_config$condition){
        group_a <- match_level(cmp[1], cond_levels)
        group_b <- match_level(cmp[2], cond_levels)
        if(group_a %in% cond_levels && group_b %in% cond_levels){
          nm <- glue::glue("sex_{sx}_main_condition_{group_a}_vs_{group_b}") %>% safe_filename()
          append_comparison_row(
            sex = sx,
            contrast_name = nm,
            contrast_type = "main_condition",
            comparison_text = glue::glue("Condition main effect: {group_a} vs {group_b}"),
            level_a = group_a,
            level_b = group_b,
            design_formula = full_design_formula
          )
        }
      }
    }

    if(!is.null(comparisons_config$age)){
      for(cmp in comparisons_config$age){
        age_a <- match_level(cmp[1], age_levels)
        age_b <- match_level(cmp[2], age_levels)
        if(age_a %in% age_levels && age_b %in% age_levels){
          nm <- glue::glue("sex_{sx}_main_age_{age_a}_vs_{age_b}") %>% safe_filename()
          append_comparison_row(
            sex = sx,
            contrast_name = nm,
            contrast_type = "main_age",
            comparison_text = glue::glue("Age main effect: {age_a} vs {age_b}"),
            level_a = age_a,
            level_b = age_b,
            design_formula = full_design_formula
          )
        }
      }
    }

    if(!is.null(comparisons_config$condition)){
      for(age_level in age_levels){
        for(cmp in comparisons_config$condition){
          group_a <- match_level(cmp[1], cond_levels)
          group_b <- match_level(cmp[2], cond_levels)
          if(group_a %in% cond_levels && group_b %in% cond_levels){
            nm <- glue::glue("sex_{sx}_age_{age_level}_condition_{group_a}_vs_{group_b}") %>% safe_filename()
            append_comparison_row(
              sex = sx,
              contrast_name = nm,
              contrast_type = "condition_within_age",
              comparison_text = glue::glue("Condition within age {age_level}: {group_a} vs {group_b}"),
              subset_factor = "age",
              subset_level = age_level,
              level_a = group_a,
              level_b = group_b,
              design_formula = "~condition"
            )
          }
        }
      }
    }

    if(!is.null(comparisons_config$age)){
      for(cond_level in cond_levels){
        for(cmp in comparisons_config$age){
          age_a <- match_level(cmp[1], age_levels)
          age_b <- match_level(cmp[2], age_levels)
          if(age_a %in% age_levels && age_b %in% age_levels){
            nm <- glue::glue("sex_{sx}_condition_{cond_level}_age_{age_a}_vs_{age_b}") %>% safe_filename()
            append_comparison_row(
              sex = sx,
              contrast_name = nm,
              contrast_type = "age_within_condition",
              comparison_text = glue::glue("Age within condition {cond_level}: {age_a} vs {age_b}"),
              subset_factor = "condition",
              subset_level = cond_level,
              level_a = age_a,
              level_b = age_b,
              design_formula = "~age"
            )
          }
        }
      }
    }

    if(!is.null(comparisons_config$condition) && !is.null(comparisons_config$age) && length(cond_levels) >= 2 && length(age_levels) >= 2){
      for(cond_level in cond_levels[-1]){
        for(age_level in age_levels[-1]){
          interaction_name <- glue::glue("condition{cond_level}.age{age_level}")
          nm <- glue::glue("sex_{sx}_interaction_{interaction_name}") %>% safe_filename()
          append_comparison_row(
            sex = sx,
            contrast_name = nm,
            contrast_type = "interaction",
            comparison_text = glue::glue("Interaction term: {interaction_name}"),
            interaction_term = interaction_name,
            design_formula = full_design_formula
          )
        }
      }
    }
  }
}

build_comparison_manifest()
if(length(comparison_rows) > 0){
  readr::write_csv(dplyr::bind_rows(comparison_rows), comparison_manifest_path)
  message(glue::glue("Wrote comparison manifest to {comparison_manifest_path}"))
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
  
  # normalized counts without transformation (for QC and potential use in downstream analyses)
  write_normalized_counts(dds_sex, out_dir_sex, sx)
  # PCA per sex (condition and age vary within sex)
  vsd_sex <- vst(dds_sex, blind = FALSE)
  write.csv(assay(vsd_sex), file.path(out_dir_sex, glue::glue("{sx}_normalized_counts.csv")))
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

  # Main effect of condition (averaged across ages within the model baseline)
  if(!is.null(comparisons_config$condition)){
    for(cmp in comparisons_config$condition){
      cond_levels <- levels(colData(dds_sex)$condition)
      group_a <- match_level(cmp[1], cond_levels)
      group_b <- match_level(cmp[2], cond_levels)
      if(!(group_a %in% cond_levels && group_b %in% cond_levels)){
        warning(glue::glue("Skipping main condition contrast; levels not found for {paste(cmp, collapse=' vs ')} (sex={sx})"))
        next
      }
      message(glue::glue("Main effect (condition): {group_a} vs {group_b} (sex={sx})"))
      res_i <- tryCatch({
        results(dds_sex, contrast = c("condition", group_a, group_b))
      }, error = function(e){
        warning(glue::glue("Failed main effect (condition) for {group_a} vs {group_b} (sex={sx}): {e$message}"))
        NULL
      })
      if(!is.null(res_i)){
        nm <- glue::glue("sex_{sx}_main_condition_{group_a}_vs_{group_b}") %>% safe_filename()
        write_normalized_counts(dds_sex, out_dir_sex, nm)
        all_res_list[[nm]] <- list(
          results = res_i,
          name = nm,
          contrast_type = "main_condition",
          group_a = group_a,
          group_b = group_b,
          design_formula = deparse(design(dds_sex))
        )
      }
    }
  }

  # Main effect of age (averaged across conditions within the model baseline)
  if(!is.null(comparisons_config$age)){
    for(cmp in comparisons_config$age){
      age_levels <- levels(colData(dds_sex)$age)
      age_a <- match_level(cmp[1], age_levels)
      age_b <- match_level(cmp[2], age_levels)
      if(!(age_a %in% age_levels && age_b %in% age_levels)){
        warning(glue::glue("Skipping main age contrast; levels not found for {paste(cmp, collapse=' vs ')} (sex={sx})"))
        next
      }
      message(glue::glue("Main effect (age): {age_a} vs {age_b} (sex={sx})"))
      res_i <- tryCatch({
        results(dds_sex, contrast = c("age", age_a, age_b))
      }, error = function(e){
        warning(glue::glue("Failed main effect (age) for {age_a} vs {age_b} (sex={sx}): {e$message}"))
        NULL
      })
      if(!is.null(res_i)){
        nm <- glue::glue("sex_{sx}_main_age_{age_a}_vs_{age_b}") %>% safe_filename()
        write_normalized_counts(dds_sex, out_dir_sex, nm)
        all_res_list[[nm]] <- list(
          results = res_i,
          name = nm,
          contrast_type = "main_age",
          age_a = age_a,
          age_b = age_b,
          design_formula = deparse(design(dds_sex))
        )
      }
    }
  }

  # Process condition comparisons stratified by age (compare R vs C within each age)
  if(!is.null(comparisons_config$condition)){
    age_levels <- levels(colData(dds_sex)$age)
    cond_levels <- levels(colData(dds_sex)$condition)
    for(age_level in age_levels){
      for(cmp in comparisons_config$condition){
        group_a <- match_level(cmp[1], cond_levels)
        group_b <- match_level(cmp[2], cond_levels)
        if(!(group_a %in% cond_levels && group_b %in% cond_levels)){
          warning(glue::glue("Skipping condition contrast at age={age_level}; levels not found for {paste(cmp, collapse=' vs ')} (sex={sx})"))
          next
        }

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
          write_normalized_counts(dds_age, out_dir_sex, nm)
          all_res_list[[nm]] <- list(results = res_i, name = nm,
                                       contrast_type = "condition_within_age",
                                       age = age_level, group_a = group_a, group_b = group_b,
                                       design_formula = deparse(design(dds_age)))
        }
      }
    }
  }

  # Process age comparisons stratified by condition (compare ages within each condition)
  if(!is.null(comparisons_config$age)){
    age_levels <- levels(colData(dds_sex)$age)
    cond_levels <- levels(colData(dds_sex)$condition)
    for(cond_level in cond_levels){
      for(cmp in comparisons_config$age){
        age_a <- match_level(cmp[1], age_levels)
        age_b <- match_level(cmp[2], age_levels)
        if(!(age_a %in% age_levels && age_b %in% age_levels)){
          warning(glue::glue("Skipping age contrast at condition={cond_level}; levels not found for {paste(cmp, collapse=' vs ')} (sex={sx})"))
          next
        }

        message(glue::glue("Age comparison: {age_a} vs {age_b} at condition={cond_level} (sex={sx})"))

        res_i <- tryCatch({
          dds_cond <- dds_sex[, colData(dds_sex)$condition == cond_level]
          if(ncol(dds_cond) >= 2 && nlevels(droplevels(colData(dds_cond)$age)) >= 2){
            colData(dds_cond)$age <- droplevels(colData(dds_cond)$age)
            if(length(age_levels_all) >= 1 && age_levels_all[1] %in% levels(colData(dds_cond)$age)){
              colData(dds_cond)$age <- relevel(colData(dds_cond)$age, ref = age_levels_all[1])
            }
            design(dds_cond) <- ~age
            dds_cond <- DESeq(dds_cond)
            results(dds_cond, contrast = c("age", age_a, age_b))
          } else {
            warning(glue::glue("Skipping {age_a} vs {age_b} at condition={cond_level} (sex={sx}): insufficient samples or age levels."))
            NULL
          }
        }, error = function(e){
          warning(glue::glue("Failed to get results for {age_a} vs {age_b} at condition={cond_level} (sex={sx}): {e$message}"))
          NULL
        })

        if(!is.null(res_i)){
          nm <- glue::glue("sex_{sx}_condition_{cond_level}_age_{age_a}_vs_{age_b}") %>% safe_filename()
          write_normalized_counts(dds_cond, out_dir_sex, nm)
          all_res_list[[nm]] <- list(results = res_i, name = nm,
                                       contrast_type = "age_within_condition",
                                       condition = cond_level, age_a = age_a, age_b = age_b,
                                       design_formula = deparse(design(dds_cond)))
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
          write_normalized_counts(dds_sex, out_dir_sex, nm)
          all_res_list[[nm]] <- list(results = res_i, name = nm,
                                       contrast_type = "interaction",
                                       description = int_name,
                                       design_formula = deparse(design(dds_sex)))
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
    res_df <- as.data.frame(res_i)
    if(nrow(res_df) == 0) next
    rds_path <- file.path(out_dir_sex, glue::glue("results_{nm}.rds"))
    csv_path <- file.path(out_dir_sex, glue::glue("results_{nm}.csv"))
    keep_idx <- rownames(res_df)
    res_i_filt <- tryCatch(res_i[keep_idx, ], error = function(e) res_i)
    saveRDS(res_i_filt, file = rds_path)
    df_i <- res_df
    df_i$gene <- rownames(df_i)
    df_i$contrast <- nm
    df_i$sex <- sx

    # map gene IDs to gene symbols using Bioconductor OrgDb
    gene_ids <- df_i$gene
    gene_ids_nover <- sub("\\.\\d+$", "", gene_ids)
    species <- NULL
    keytype <- 'ENSEMBL'
    if(any(grepl('^ENSMUSG', gene_ids_nover))){ species <- 'mouse'; keytype <- 'ENSEMBL' }
    else if(any(grepl('^ENSG', gene_ids_nover))){ species <- 'human'; keytype <- 'ENSEMBL' }
    else if(any(grepl('^FBgn', gene_ids_nover))){ species <- 'drosophila'; keytype <- 'FLYBASE' }
    else if(any(grepl('^ENSMUSDM', gene_ids_nover))){ species <- 'drosophila'; keytype <- 'ENSEMBL' }
    df_i$gene_symbol <- NA_character_
    if(!is.null(species)){
      OrgDb <- NULL
      if(species == 'mouse' && requireNamespace('org.Mm.eg.db', quietly=TRUE)) OrgDb <- get('org.Mm.eg.db', envir = asNamespace('org.Mm.eg.db'))
      if(species == 'human' && requireNamespace('org.Hs.eg.db', quietly=TRUE)) OrgDb <- get('org.Hs.eg.db', envir = asNamespace('org.Hs.eg.db'))
      if(species == 'drosophila' && requireNamespace('org.Dm.eg.db', quietly=TRUE)) OrgDb <- get('org.Dm.eg.db', envir = asNamespace('org.Dm.eg.db'))
      if(!is.null(OrgDb)) {
        map_df <- AnnotationDbi::select(OrgDb, keys = unique(gene_ids_nover), keytype = keytype, columns = c('SYMBOL'))
        if(!is.null(map_df) && nrow(map_df) > 0) {
          if(keytype == 'FLYBASE') map_df <- dplyr::rename(map_df, gene_nover = FLYBASE, gene_symbol = SYMBOL)
          if(keytype == 'ENSEMBL') map_df <- dplyr::rename(map_df, gene_nover = ENSEMBL, gene_symbol = SYMBOL)
          df_i$gene_nover <- gene_ids_nover
          df_i <- dplyr::left_join(df_i, map_df, by = c('gene_nover'))
          gs_cols <- c("gene_symbol.x", "gene_symbol.y", "gene_symbol")
          present_cols <- gs_cols[gs_cols %in% colnames(df_i)]
          if(length(present_cols) > 0) {
            tmp <- NULL
            for(col in present_cols) tmp <- if(is.null(tmp)) df_i[[col]] else dplyr::coalesce(tmp, df_i[[col]])
            df_i$gene_symbol <- tmp
            df_i <- dplyr::select(df_i, -dplyr::any_of(c("gene_symbol.x", "gene_symbol.y")))
          }
        }
      }
    }
    if("gene_nover" %in% colnames(df_i)) df_i$gene_nover <- NULL
    if(!"gene_symbol" %in% colnames(df_i)) df_i$gene_symbol <- NA_character_
    other_cols <- setdiff(colnames(df_i), c("gene", "gene_symbol"))
    df_i <- df_i[, c("gene", "gene_symbol", other_cols), drop = FALSE]
    if("padj" %in% colnames(df_i)) df_i <- df_i[order(df_i$padj, na.last = TRUE), , drop = FALSE]
    write.csv(df_i, file = csv_path, row.names = FALSE)
    final_res_list[[nm]] <- df_i

    # --- Write filterThreshold report ---
    filter_threshold <- tryCatch({
      mt <- metadata(res_i)
      if(!is.null(mt$filterThreshold)) mt$filterThreshold else NA
    }, error = function(e) NA)
    design_formula <- if(!is.null(res_item$design_formula)) res_item$design_formula else NA
    report_path <- file.path(out_dir_sex, glue::glue("results_{nm}_filterThreshold.txt"))
    report_lines <- c(
      paste0("design: ", toString(design_formula)),
      paste0("filterThreshold: ", toString(names(filter_threshold)), " ", toString(filter_threshold))
    )
    writeLines(report_lines, con = report_path)
  }

  # Save the DESeq dataset for this sex
  saveRDS(dds_sex, file=file.path(out_dir_sex, "dds.rds"))
}
