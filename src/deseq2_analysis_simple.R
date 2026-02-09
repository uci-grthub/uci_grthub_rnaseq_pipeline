suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(glue)
  library(AnnotationDbi)
  library(yaml)
  library(iSEEde)
  library(ggrepel)
})

dir.create(out_dir, showWarnings=FALSE)
dplyr::right_join(
  tibble(index_pair = str_extract(colnames(count_matrix), "[ACGT]+-[ACGT]+")),
  by = c("index_pair")
) |>
tidyr::separate(sample_name, into = c("experiment", "age", "sex", "condition", "replicate"), sep = "_", remove = FALSE) |>
tibble::column_to_rownames("sample_name") |>
  identity()
tibble::rownames_to_column("ensgene") |>
  readr::write_csv(file.path(out_dir, "all_sample_counts.csv"))
dds <- DESeqDataSetFromMatrix(countData=count_matrix, 
  colData=meta, 
  design=~condition+age+condition:age)
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
age_levels_all <- levels(meta$age)

# --- Argument parsing and defaults ---
args <- commandArgs(trailingOnly=TRUE)
counts_file <- if(length(args) >= 1 && nzchar(args[1])) args[1] else "output/feature_count/all_samples_counts.txt"
meta_file   <- if(length(args) >= 2 && nzchar(args[2])) args[2] else "metadata/metadata.csv"
out_dir     <- if(length(args) >= 3 && nzchar(args[3])) args[3] else "output/deseq2"
dir.create(out_dir, showWarnings=FALSE)

# --- Read and format data ---
counts <- read.table(counts_file, header=TRUE, row.names=1)
count_matrix <- counts[,6:ncol(counts)]
format_sample_id <- function(colname) {
  colname %>%
    str_remove("results.hisat2_alignment\\.") %>%
    str_remove("_align_sorted.bam") %>%
    str_replace_all("\\.", "-")
}
colnames(count_matrix) <- sapply(colnames(count_matrix), format_sample_id)

meta <- read.csv(meta_file) %>%
  janitor::clean_names() %>%
  mutate(index_pair = glue("{i7_barcode_sequence}-{i5_barcode_sequence}")) %>%
  right_join(tibble(index_pair = str_extract(colnames(count_matrix), "[ACGT]+-[ACGT]+")), by = c("index_pair")) %>%
  tidyr::separate(sample_name, into = c("experiment", "age", "sex", "condition", "replicate"), sep = "_", remove = FALSE) %>%
  tibble::column_to_rownames("sample_name")

meta$condition <- factor(meta$condition)
if("C" %in% levels(meta$condition)) meta$condition <- relevel(meta$condition, ref = "C")
meta$age <- factor(meta$age)
if("10" %in% levels(meta$age)) meta$age <- relevel(meta$age, ref = "10")
meta$sex <- factor(meta$sex)
colnames(count_matrix) <- rownames(meta)
count_matrix <- count_matrix[, rownames(meta)]

# --- Save all sample counts ---
count_matrix %>%
  tibble::rownames_to_column("ensgene") %>%
  readr::write_csv(file.path(out_dir, "all_sample_counts.csv"))

# --- DESeq2 setup and PCA (all samples) ---
dds <- DESeqDataSetFromMatrix(countData=count_matrix, colData=meta, design=~condition+age+condition:age)
dds_all <- DESeq(dds)
vsd_all <- vst(dds_all, blind = FALSE)
plot_var_all <- c("condition", "age", "sex")
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

# --- Comparison config ---
comparisons_yaml <- file.path(dirname(meta_file), "comparisons.yaml")
if(file.exists(comparisons_yaml)) {
  comparisons_config <- yaml::read_yaml(comparisons_yaml)
} else {
  comparisons_config <- list(condition = list(c("R", "C")))
  dir.create(dirname(meta_file), showWarnings = FALSE, recursive = TRUE)
  yaml::write_yaml(comparisons_config, comparisons_yaml)
}

cond_levels_all <- levels(meta$condition)
age_levels_all <- levels(meta$age)
if(is.null(comparisons_config$condition) || length(comparisons_config$condition) == 0) {
  if(length(cond_levels_all) >= 2) comparisons_config$condition <- list(c(cond_levels_all[2], cond_levels_all[1]))
}
if(is.null(comparisons_config$age) || length(comparisons_config$age) == 0) {
  if(length(age_levels_all) >= 2) comparisons_config$age <- list(c(age_levels_all[1], age_levels_all[2]))
}

safe_filename <- function(x) {
  x %>% stringr::str_replace_all("\\+", "plus") %>% stringr::str_replace_all("[^A-Za-z0-9_-]", "_") %>% tolower()
}
match_level <- function(x, levs) {
  idx <- match(tolower(x), tolower(levs)); if(is.na(idx)) x else levs[idx]
}

# --- Run DESeq2 per sex ---
sex_levels <- sort(unique(as.character(meta$sex)))
for(sx in sex_levels) {
  out_dir_sex <- file.path(out_dir, glue("sex_{sx}"))
  dir.create(out_dir_sex, showWarnings = FALSE, recursive = TRUE)
  dds_sex <- dds[, colData(dds)$sex == sx]
  dds_sex <- DESeq(dds_sex)

  # PCA per sex
  vsd_sex <- vst(dds_sex, blind = FALSE)
  plot_var <- c("condition", "age")
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
  pdf(file.path("results", glue("pca_plots_{sx}.pdf")))
  print(pca_plots)
  dev.off()

  all_res_list <- list()

  # Main effect of condition
  if(!is.null(comparisons_config$condition)) {
    for(cmp in comparisons_config$condition) {
      cond_levels <- levels(colData(dds_sex)$condition)
      group_a <- match_level(cmp[1], cond_levels)
      group_b <- match_level(cmp[2], cond_levels)
      if(group_a %in% cond_levels && group_b %in% cond_levels) {
        res_i <- results(dds_sex, contrast = c("condition", group_a, group_b))
        nm <- glue("sex_{sx}_main_condition_{group_a}_vs_{group_b}") %>% safe_filename()
        all_res_list[[nm]] <- list(results = res_i, name = nm, contrast_type = "main_condition", group_a = group_a, group_b = group_b)
      }
    }
  }

  # Main effect of age
  if(!is.null(comparisons_config$age)) {
    for(cmp in comparisons_config$age) {
      age_levels <- levels(colData(dds_sex)$age)
      age_a <- match_level(cmp[1], age_levels)
      age_b <- match_level(cmp[2], age_levels)
      if(age_a %in% age_levels && age_b %in% age_levels) {
        res_i <- results(dds_sex, contrast = c("age", age_a, age_b))
        nm <- glue("sex_{sx}_main_age_{age_a}_vs_{age_b}") %>% safe_filename()
        all_res_list[[nm]] <- list(results = res_i, name = nm, contrast_type = "main_age", age_a = age_a, age_b = age_b)
      }
    }
  }

  # Condition comparisons within each age
  if(!is.null(comparisons_config$condition)) {
    age_levels <- levels(colData(dds_sex)$age)
    cond_levels <- levels(colData(dds_sex)$condition)
    for(age_level in age_levels) {
      for(cmp in comparisons_config$condition) {
        group_a <- match_level(cmp[1], cond_levels)
        group_b <- match_level(cmp[2], cond_levels)
        dds_age <- dds_sex[, colData(dds_sex)$age == age_level]
        if(ncol(dds_age) >= 2 && nlevels(droplevels(colData(dds_age)$condition)) >= 2) {
          colData(dds_age)$condition <- droplevels(colData(dds_age)$condition)
          if("C" %in% levels(colData(dds_age)$condition)) colData(dds_age)$condition <- relevel(colData(dds_age)$condition, ref = "C")
          design(dds_age) <- ~condition
          dds_age <- DESeq(dds_age)
          res_i <- results(dds_age, contrast = c("condition", group_a, group_b))
          nm <- glue("sex_{sx}_age_{age_level}_condition_{group_a}_vs_{group_b}") %>% safe_filename()
          all_res_list[[nm]] <- list(results = res_i, name = nm, contrast_type = "condition_within_age", age = age_level, group_a = group_a, group_b = group_b)
        }
      }
    }
  }

  # Age comparisons within each condition
  if(!is.null(comparisons_config$age)) {
    age_levels <- levels(colData(dds_sex)$age)
    cond_levels <- levels(colData(dds_sex)$condition)
    for(cond_level in cond_levels) {
      for(cmp in comparisons_config$age) {
        age_a <- match_level(cmp[1], age_levels)
        age_b <- match_level(cmp[2], age_levels)
        dds_cond <- dds_sex[, colData(dds_sex)$condition == cond_level]
        if(ncol(dds_cond) >= 2 && nlevels(droplevels(colData(dds_cond)$age)) >= 2) {
          colData(dds_cond)$age <- droplevels(colData(dds_cond)$age)
          if(length(age_levels_all) >= 1 && age_levels_all[1] %in% levels(colData(dds_cond)$age)) colData(dds_cond)$age <- relevel(colData(dds_cond)$age, ref = age_levels_all[1])
          design(dds_cond) <- ~age
          dds_cond <- DESeq(dds_cond)
          res_i <- results(dds_cond, contrast = c("age", age_a, age_b))
          nm <- glue("sex_{sx}_condition_{cond_level}_age_{age_a}_vs_{age_b}") %>% safe_filename()
          all_res_list[[nm]] <- list(results = res_i, name = nm, contrast_type = "age_within_condition", condition = cond_level, age_a = age_a, age_b = age_b)
        }
      }
    }
  }

  # Interaction effect: condition:age
  if(!is.null(comparisons_config$condition) && !is.null(comparisons_config$age)) {
    available_names <- resultsNames(dds_sex)
    interaction_names <- available_names[grepl("condition.*age", available_names, ignore.case = TRUE)]
    for(int_name in interaction_names) {
      res_i <- results(dds_sex, name = int_name)
      nm <- glue("sex_{sx}_interaction_{int_name}") %>% safe_filename()
      all_res_list[[nm]] <- list(results = res_i, name = nm, contrast_type = "interaction", description = int_name)
    }
  }

  # --- Save results ---
  for(nm in names(all_res_list)) {
    res_i <- all_res_list[[nm]]$results
    res_df <- as.data.frame(res_i)
    if("pvalue" %in% colnames(res_df)) res_df <- res_df[!is.na(res_df$pvalue), , drop = FALSE]
    if(nrow(res_df) == 0) next
    rds_path <- file.path(out_dir_sex, glue("results_{nm}.rds"))
    csv_path <- file.path(out_dir_sex, glue("results_{nm}.csv"))
    saveRDS(res_i, file = rds_path)
    df_i <- res_df
    df_i$gene <- rownames(df_i)
    df_i$contrast <- nm
    df_i$sex <- sx

    # --- Gene symbol mapping ---
    gene_ids <- df_i$gene
    gene_ids_nover <- sub("\\.\\d+$", "", gene_ids)
    species <- NULL
    keytype <- 'ENSEMBL'
    if(any(grepl('^ENSMUSG', gene_ids_nover))) { species <- 'mouse'; keytype <- 'ENSEMBL' }
    else if(any(grepl('^ENSG', gene_ids_nover))) { species <- 'human'; keytype <- 'ENSEMBL' }
    else if(any(grepl('^FBgn', gene_ids_nover))) { species <- 'drosophila'; keytype <- 'FLYBASE' }
    else if(any(grepl('^ENSMUSDM', gene_ids_nover))) { species <- 'drosophila'; keytype <- 'ENSEMBL' }
    df_i$gene_symbol <- NA_character_
    if(!is.null(species)) {
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
  }
  saveRDS(dds_sex, file=file.path(out_dir_sex, "dds.rds"))
}
