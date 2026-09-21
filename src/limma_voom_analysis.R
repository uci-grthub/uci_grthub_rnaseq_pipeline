#!/usr/bin/env Rscript
# limma-voom differential expression analysis
#
# Usage: Rscript limma_voom_analysis.R counts.txt metadata.csv output_dir \
#            comparisons.yaml [block_var] [quality_weights]
#
# Design. Every animal (`mouse_id`, e.g. "SOM1") contributed one sample from
# each of three brain regions, so region is a within-animal factor and the
# samples are not independent. The comparisons YAML assigns animals to named
# line_groups; each `comparisons` entry optionally subsets the metadata
# (`subset:`, e.g. to one region) and then contrasts two line_groups.
#
# Model. ONE global fit over all samples: filterByExpr -> TMM -> voom -> lmFit
# on a no-intercept design over the age_group x region cells, with the animal as
# a random intercept via duplicateCorrelation. Each comparison becomes a
# contrast of that single fit -- (mean of the cells its group_a side selects)
# minus (mean of the group_b cells) -- so a positive logFC still means "up in
# group_a".
#
# Why pooled rather than one fit per comparison. The earlier version subset to
# one region per comparison and fitted 3 vs 3 independently, which left ~4
# residual df and made blocking impossible (one sample per animal per subset).
# Pooling raises that to ~18 residual df. The blocking is what makes spending
# those df legitimate: once the three regions are pooled, 27 samples are still
# only 9 independent animals, so an unblocked pooled fit would be
# pseudo-replication. See the caveats above the fit for what pooling costs.

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(tidyverse)
  library(glue)
  library(AnnotationDbi)
  library(yaml)
  library(ggrepel)
})

# Sex inference, the sample-tracking checks and the counts/metadata loaders are
# shared with the standalone QC entry point, so they live in one file. Resolve
# it relative to this script rather than the working directory, which Snakemake
# does not guarantee.
local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  source(file.path(if (length(f) == 1) dirname(normalizePath(f)) else "src", "infer_sex.R"))
})

args <- commandArgs(trailingOnly = TRUE)
arg_or <- function(i, default) if (length(args) >= i && nzchar(args[i])) args[i] else default

counts_file <- arg_or(1, "output/feature_count/all_samples_counts.txt")
meta_file <- arg_or(2, "metadata/metadata.csv")
out_dir <- arg_or(3, "output/limma_voom")
comparisons_config_path <- arg_or(4, "src/de_comparisons.yaml")

# duplicateCorrelation blocking factor, or "none" to fit samples as independent.
block_var <- arg_or(5, "mouse_id")
use_dupcor <- !identical(tolower(block_var), "none")

# Sample weighting. "per_sample" gives every sample its own arrayWeight,
# "per_region" pools the weight over samples sharing a brain region, "none"
# uses plain voom. This is a knob rather than a constant because the three modes
# give materially different gene counts on this data -- see the caveats above
# the fit.
weight_mode <- tolower(arg_or(6, "per_sample"))
if (!weight_mode %in% c("per_sample", "per_region", "none")) {
  stop(glue("quality_weights must be per_sample, per_region or none; got '{weight_mode}'"))
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

## ---- Counts and metadata ---------------------------------------------------

loaded <- load_sample_metadata(meta_file, load_counts(counts_file))
meta <- loaded$meta
count_matrix <- loaded$counts

## ---- Helpers ---------------------------------------------------------------

safe_filename <- function(x) tolower(str_replace_all(x, "[^A-Za-z0-9_-]", "_"))

pca_from_logcpm <- function(logcpm, coldata, colour_var, labels, title) {
  rv <- matrixStats::rowVars(as.matrix(logcpm))
  top <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
  pca <- prcomp(t(logcpm[top, , drop = FALSE]))
  percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2))

  tibble::tibble(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    group = as.character(coldata[[colour_var]]),
    sample_label = as.character(labels)
  ) |>
    ggplot(aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = sample_label), size = 3, max.overlaps = 20) +
    xlab(glue("PC1: {percent_var[1]}% variance")) +
    ylab(glue("PC2: {percent_var[2]}% variance")) +
    labs(title = title, color = colour_var) +
    theme_bw()
}

annotate_gene_symbols <- function(df) {
  gene_ids <- sub("\\.\\d+$", "", df$gene)
  df$gene_symbol <- NA_character_
  org_db <- c(ENSMUSG = "org.Mm.eg.db", ENSG = "org.Hs.eg.db", FBgn = "org.Dm.eg.db")
  prefix <- names(org_db)[vapply(names(org_db), \(p) any(startsWith(gene_ids, p)), logical(1))]

  if (length(prefix) == 1 && requireNamespace(org_db[[prefix]], quietly = TRUE)) {
    keytype <- if (prefix == "FBgn") "FLYBASE" else "ENSEMBL"
    mapped <- AnnotationDbi::select(
      get(org_db[[prefix]], envir = asNamespace(org_db[[prefix]])),
      keys = unique(gene_ids), keytype = keytype, columns = "SYMBOL"
    )
    if (nrow(mapped) > 0) {
      df$gene_symbol <- mapped$SYMBOL[match(gene_ids, mapped[[keytype]])]
    }
  }

  df[, c("gene", "gene_symbol", setdiff(colnames(df), c("gene", "gene_symbol"))), drop = FALSE]
}

# duplicateCorrelation is only worth fitting when the blocking factor is both
# replicated and separable from the design. Unlike the old per-comparison
# subsets -- where each animal contributed exactly one sample and blocking was
# always (correctly) dropped -- this runs once against the full sample set and
# the global cell-means design. Returns the block vector, or NULL having said
# why.
resolve_block <- function(meta_all, design) {
  if (!use_dupcor) {
    message("block_var = 'none': fitting all samples as independent")
    return(NULL)
  }
  if (!block_var %in% colnames(meta_all)) {
    stop(glue("Blocking column '{block_var}' not found in metadata. ",
              "Available: {paste(colnames(meta_all), collapse = ', ')}"))
  }

  block <- as.character(meta_all[[block_var]])

  if (anyNA(block) || any(!nzchar(block))) {
    warning(glue("Blocking column '{block_var}' has ",
                 "{sum(is.na(block) | !nzchar(block))} missing value(s); fitting unblocked"))
    return(NULL)
  }

  n_blocks <- length(unique(block))
  if (n_blocks < 2) {
    warning(glue("'{block_var}' gives {n_blocks} block(s); fitting unblocked"))
    return(NULL)
  }

  # Without at least one block holding more than one sample there is no
  # within-block information to estimate a correlation from.
  if (max(table(block)) < 2) {
    warning(glue("'{block_var}' gives one sample per block across {length(block)} ",
                 "samples; nothing to correlate, fitting unblocked"))
    return(NULL)
  }

  # Separability. If every block indicator already lies in the column space of
  # the design, the random intercept is confounded with the fixed effects and
  # would absorb the very contrasts being tested. Here each animal spans three
  # different cells, so it adds rank and the check passes.
  joint <- cbind(design, model.matrix(~ 0 + factor(block)))
  if (qr(joint)$rank <= qr(design)$rank) {
    warning(glue("'{block_var}' is confounded with the design (adds no rank beyond ",
                 "the {qr(design)$rank} design columns); fitting unblocked"))
    return(NULL)
  }

  message(glue("Blocking on '{block_var}': {n_blocks} blocks, ",
               "{paste(range(table(block)), collapse = '-')} samples each"))
  block
}

# voom, with the weighting mode the caller asked for. `...` carries block= and
# correlation= through to voom() in every mode: limma's voomWithQualityWeights
# declares neither argument itself but forwards its dots into both of its
# internal voom() calls.
run_voom <- function(dge, design, plot = FALSE, ...) {
  switch(weight_mode,
    none = voom(dge, design, plot = plot, ...),
    per_region = voomWithQualityWeights(dge, design, var.group = meta$condition,
                                        plot = plot, ...),
    per_sample = voomWithQualityWeights(dge, design, plot = plot, ...)
  )
}

## ---- Inferred sex ----------------------------------------------------------
## Runs on the unfiltered matrix, before filterByExpr below: in a single-sex
## experiment the Y genes would not survive the filter, and their absence is
## the signal.

meta <- annotate_and_write_sex(meta, count_matrix, out_dir)
invisible(check_sex_coherence(meta, block_var = block_var))

## ---- Global design, filter and normalisation -------------------------------

meta$cell <- factor(paste(meta$age_group, meta$condition, sep = "_"))
design <- model.matrix(~ 0 + meta$cell)
colnames(design) <- levels(meta$cell)
rownames(design) <- rownames(meta)

if (qr(design)$rank < ncol(design)) {
  stop(glue("Global design is rank-deficient ({qr(design)$rank} < {ncol(design)} cells); ",
            "some age_group x condition cell has no samples"))
}

dge <- DGEList(counts = count_matrix, samples = meta)
# Filtering on the cell factor keeps genes expressed in only one region, which
# a filter on age group alone would drop. Note the consequence: all contrasts
# now share one gene set and one BH universe, where previously each comparison
# filtered separately and their adjusted p-values were not on a common footing.
dge <- dge[filterByExpr(dge, group = meta$cell), , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
message(glue("Kept {nrow(dge)} genes across {ncol(dge)} samples ",
             "({ncol(design)} age_group x condition cells)"))

## ---- Ordination across all samples -----------------------------------------
## On the same TMM log-CPM the contrasts are fitted on.

logcpm_all <- cpm(dge, log = TRUE, prior.count = 3)

# Points are labelled "<age group>-<condition>", so whichever variable a panel
# is coloured by, the label still carries the other.
pca_labels <- paste(meta$age_group, meta$condition, sep = "-")
pca_plots <- c("condition", "age_group", "line", "inferred_sex") |>
  set_names() |>
  map(\(v) pca_from_logcpm(logcpm_all, meta, v, pca_labels, glue("PCA (log-CPM) - all samples - {v}")))

pdf(file.path("results", "pca_plots_all_limma_voom.pdf"), width = 12, height = 6)
print(pca_plots)
dev.off()

# PNGs for the project report; the PDF above keeps every panel in vector form.
for (v in c("condition", "age_group", "inferred_sex")) {
  ggsave(file.path(out_dir, glue("pca_all_samples_{v}.png")), pca_plots[[v]],
         width = 7, height = 5, dpi = 200)
}

## ---- Global blocked fit ----------------------------------------------------
##
## Caveats worth knowing before reading the results:
##
##  * The contrasts of interest are BETWEEN animals (age group varies between
##    animals only), so a positive intra-animal correlation INFLATES their
##    variance rather than shrinking it. Blocking here is a correctness fix
##    against pseudo-replication introduced by pooling, not a power gain; the
##    power comes from the pooled residual df, and blocking is what makes
##    spending those df legitimate. This is the easiest thing to get backwards.
##  * eBayes now shrinks every gene toward a single prior fitted across three
##    biologically distinct regions. A gene with genuinely region-dependent
##    dispersion is served worse than by the old per-region fits.
##  * duplicateCorrelation estimates ONE scalar correlation for the whole
##    matrix, not a per-gene value.
##  * Sample weighting is a real modelling choice, not a technical detail. With
##    only 3 samples per cell, a per-sample weight cannot distinguish a
##    technically noisy sample from a biological outlier, so down-weighting one
##    removes within-group variance that may be real. weight_mode = per_region
##    is the conservative alternative and "none" reproduces plain voom.

block <- resolve_block(meta, design)
voom_path <- file.path(out_dir, "voom_mean_variance.pdf")
# voomWithQualityWeights draws two panels (trend plus a sample-weight barplot);
# plain voom draws one.
voom_dims <- if (weight_mode == "none") c(7, 6) else c(11, 5)
message(glue("Sample weighting: {weight_mode}"))

if (is.null(block)) {
  pdf(voom_path, width = voom_dims[1], height = voom_dims[2])
  v <- run_voom(dge, design, plot = TRUE)
  dev.off()
  fit <- lmFit(v, design)
  consensus_cor <- NA_real_
} else {
  # The precision weights and the consensus correlation each depend on the
  # other, so the limma User's Guide iterates the pair once: unblocked voom ->
  # correlation -> blocked voom -> correlation, and fits on the second estimate.
  v <- run_voom(dge, design)
  dc <- duplicateCorrelation(v, design, block = block)
  pdf(voom_path, width = voom_dims[1], height = voom_dims[2])
  v <- run_voom(dge, design, plot = TRUE, block = block,
                correlation = dc$consensus.correlation)
  dev.off()
  consensus_cor <- duplicateCorrelation(v, design, block = block)$consensus.correlation
  message(glue("duplicateCorrelation consensus ({block_var}): {round(consensus_cor, 4)}"))
  fit <- lmFit(v, design, block = block, correlation = consensus_cor)
}

if (!is.null(v$targets$sample.weights)) {
  message(glue("Sample weights range {round(min(v$targets$sample.weights), 3)} - ",
               "{round(max(v$targets$sample.weights), 3)}"))
}

## ---- Contrasts -------------------------------------------------------------

config <- yaml::read_yaml(comparisons_config_path)
stopifnot(
  "comparisons config must define line_groups" = !is.null(config$line_groups),
  "comparisons config must define comparisons" = !is.null(config$comparisons)
)

# A comparison names two line_groups and optionally subsets the metadata
# (`subset:`, e.g. to one brain region). The samples a side selects fall into
# one or more cells of the global design, and the side's estimate is the mean of
# those cell means -- weighting cells equally rather than samples, so an
# unbalanced subset cannot let the larger cell dominate.
cells_for_side <- function(lines, subset_spec, cmp_name, side_label) {
  keep <- rep(TRUE, nrow(meta))
  for (subset_col in names(subset_spec)) {
    if (!subset_col %in% colnames(meta)) {
      stop(glue("{cmp_name} subsets on '{subset_col}', which is not a column in {meta_file}"))
    }
    wanted <- as.character(unlist(subset_spec[[subset_col]]))
    present <- unique(as.character(meta[[subset_col]]))
    unknown <- setdiff(wanted, present)
    if (length(unknown) > 0) {
      warning(glue("{cmp_name}: subset {subset_col} = {paste(unknown, collapse = ', ')} ",
                   "matches no sample (present: {paste(present, collapse = ', ')})"))
    }
    keep <- keep & as.character(meta[[subset_col]]) %in% wanted
  }
  keep <- keep & meta$animal %in% lines

  if (!any(keep)) {
    warning(glue("{cmp_name}: {side_label} selects no samples; skipping comparison"))
    return(NULL)
  }
  list(cells = unique(as.character(meta$cell[keep])), samples = rownames(meta)[keep])
}

build_contrast <- function(cmp) {
  lines_a <- config$line_groups[[cmp$group_a]]
  lines_b <- config$line_groups[[cmp$group_b]]
  if (is.null(lines_a) || is.null(lines_b)) {
    warning(glue("Skipping {cmp$name}: unknown line group(s) {cmp$group_a}/{cmp$group_b}"))
    return(NULL)
  }

  a <- cells_for_side(lines_a, cmp$subset, cmp$name, glue("group_a ({cmp$group_a})"))
  b <- cells_for_side(lines_b, cmp$subset, cmp$name, glue("group_b ({cmp$group_b})"))
  if (is.null(a) || is.null(b)) return(NULL)

  # A cell on both sides would silently cancel out of the contrast, which is
  # never what a comparison means -- it signals overlapping line_groups.
  shared <- intersect(a$cells, b$cells)
  if (length(shared) > 0) {
    warning(glue("Skipping {cmp$name}: cell(s) {paste(shared, collapse = ', ')} appear on ",
                 "both sides (are {cmp$group_a} and {cmp$group_b} disjoint?)"))
    return(NULL)
  }

  cvec <- setNames(numeric(ncol(design)), colnames(design))
  cvec[a$cells] <- 1 / length(a$cells)
  cvec[b$cells] <- -1 / length(b$cells)

  list(name = cmp$name, group_a = cmp$group_a, group_b = cmp$group_b, contrast = cvec,
       samples_a = a$samples, samples_b = b$samples)
}

contrast_specs <- Filter(Negate(is.null), lapply(config$comparisons, build_contrast))
if (length(contrast_specs) == 0) {
  stop(glue("No usable comparisons in {comparisons_config_path}"))
}

contrast_matrix <- do.call(cbind, lapply(contrast_specs, `[[`, "contrast"))
colnames(contrast_matrix) <- vapply(contrast_specs, `[[`, character(1), "name")
rownames(contrast_matrix) <- colnames(design)

efit <- eBayes(contrasts.fit(fit, contrast_matrix))
residual_df <- unique(fit$df.residual)[1]
message(glue("Fitted {ncol(contrast_matrix)} contrast(s); residual df {residual_df}, ",
             "df.total {round(unique(efit$df.total)[1], 1)}"))

## ---- Per-comparison outputs ------------------------------------------------

manifest_rows <- list()

for (spec in contrast_specs) {
  cmp_name <- spec$name
  message(glue("Writing {cmp_name}: {spec$group_a} vs {spec$group_b}"))

  res_df <- topTable(efit, coef = cmp_name, number = Inf, sort.by = "P")
  res_df$gene <- rownames(res_df)
  res_df <- annotate_gene_symbols(res_df)

  cmp_dir <- file.path(out_dir, cmp_name)
  dir.create(cmp_dir, showWarnings = FALSE, recursive = TRUE)

  # 0.1 is an exploratory cut, used here only to stamp the file names; the
  # project report recounts significance at its own threshold.
  n_sig <- sum(res_df$adj.P.Val < 0.1, na.rm = TRUE)
  base_name <- safe_filename(glue("deg_{spec$group_a}_v_{spec$group_b}_{n_sig}"))
  results_path <- file.path(cmp_dir, glue("{base_name}.csv"))
  fit_path <- file.path(cmp_dir, glue("{base_name}_efit.rds"))

  write.csv(res_df, results_path, row.names = FALSE)
  # The single-contrast slice of the global fit, so a per-comparison RDS still
  # means what it used to while carrying the global df.prior/s2.prior.
  saveRDS(efit[, cmp_name], fit_path)

  # Log-CPM for this contrast's samples, taken from the global voom rather than
  # renormalised, so every comparison is on one common scale.
  cmp_samples <- c(spec$samples_a, spec$samples_b)
  write.csv(v$E[, cmp_samples, drop = FALSE],
            file.path(cmp_dir, glue("{base_name}_logcpm.csv")))

  # Within a comparison both age group and condition are constant per side, so
  # points are labelled by animal instead.
  side <- data.frame(
    line_group = factor(ifelse(cmp_samples %in% spec$samples_a, spec$group_a, spec$group_b),
                        levels = c(spec$group_a, spec$group_b)),
    row.names = cmp_samples
  )
  pca_plot <- tryCatch(
    pca_from_logcpm(v$E[, cmp_samples, drop = FALSE], side, "line_group",
                    meta[cmp_samples, "animal"], glue("PCA (log-CPM) - {cmp_name}")),
    error = function(e) {
      warning(glue("PCA failed for {cmp_name}: {e$message}"))
      NULL
    }
  )
  if (!is.null(pca_plot)) {
    pdf(file.path(cmp_dir, glue("{base_name}_pca.pdf")), width = 8, height = 6)
    print(pca_plot)
    dev.off()
  }

  manifest_rows[[length(manifest_rows) + 1]] <- tibble::tibble(
    comparison = cmp_name,
    group_a = spec$group_a,
    group_b = spec$group_b,
    # n_samples is this contrast's samples; n_samples_in_fit is the pooled model
    # every contrast is drawn from.
    n_samples = length(cmp_samples),
    n_samples_in_fit = ncol(dge),
    n_group_a = length(spec$samples_a),
    n_group_b = length(spec$samples_b),
    n_genes_after_filter = nrow(dge),
    block_var = if (is.null(block)) NA_character_ else block_var,
    consensus_correlation = consensus_cor,
    quality_weights = weight_mode,
    residual_df = residual_df,
    df_total = unique(efit$df.total)[1],
    n_sig_adj_p_0_1 = n_sig,
    results_csv = results_path,
    efit_rds = fit_path
  )
}

manifest <- if (length(manifest_rows) > 0) dplyr::bind_rows(manifest_rows) else tibble::tibble()
readr::write_csv(manifest, file.path(out_dir, "limma_voom_comparisons_manifest.csv"))

## ---- Fixed-path outputs for downstream rules -------------------------------

# The whole multi-contrast fit, which is what any interactive follow-up wants.
saveRDS(efit, file.path(out_dir, "efit.rds"))

primary_cmp <- config$comparisons[[1]]$name
primary_row <- dplyr::filter(manifest, comparison == primary_cmp)
if (nrow(primary_row) > 0) {
  invisible(file.copy(primary_row$results_csv[1],
                      file.path(out_dir, "limma_voom_results.csv"), overwrite = TRUE))
} else {
  warning(glue("{primary_cmp} not in manifest; limma_voom_results.csv not written"))
}
