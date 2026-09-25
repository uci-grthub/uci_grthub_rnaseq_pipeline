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
  d <- if (length(f) == 1) dirname(normalizePath(f)) else "src"
  source(file.path(d, "infer_sex.R"))
  # Gene symbols and biotypes, from the EnsDb built off the counting GTF.
  source(file.path(d, "gene_annotation.R"))
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

# Gene-annotation EnsDb, built once per reference by
# src/build_gene_annotation.sh and shared across projects.
ensdb_path <- arg_or(7, "")
if (!nzchar(ensdb_path) || !file.exists(ensdb_path)) {
  stop(glue("Gene annotation database not found: '{ensdb_path}'. ",
            "Build it with: bash src/build_gene_annotation.sh <species>"))
}
gene_ann <- read_gene_annotation(ensdb_path)
message(glue("Gene annotation: {nrow(gene_ann)} genes from {ensdb_path}"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

## ---- Counts and metadata ---------------------------------------------------

loaded <- load_sample_metadata(meta_file, load_counts(counts_file))
meta <- loaded$meta
count_matrix <- loaded$counts

## ---- Helpers ---------------------------------------------------------------

safe_filename <- function(x) tolower(str_replace_all(x, "[^A-Za-z0-9_-]", "_"))

# Multidimensional scaling as limma defines it, not PCA. limma::plotMDS uses
# the leading log-fold-change distance: for EVERY PAIR of samples it takes that
# pair's own top 500 most different genes and computes the root-mean-square of
# their log-fold changes, then scales the resulting distance matrix. PCA instead
# picks one global set of top-variable genes and projects every sample onto the
# same axes.
#
# The difference matters when samples differ from each other for different
# reasons -- a pair separated by a handful of genes no other pair cares about
# still shows up as distant here, where PCA can miss it because those genes
# never made the global top-variable list. This is the standard ordination in
# the limma-voom workflow, which is why it replaces the PCA panels.
#
# Coordinates come from plot = FALSE so the figure is drawn with ggplot and
# keeps the repel labelling the rest of this script uses.
mds_from_logcpm <- function(logcpm, coldata, colour_var, labels, title,
                            legend_title = colour_var) {
  mds <- limma::plotMDS(as.matrix(logcpm), top = 500,
                        gene.selection = "pairwise", plot = FALSE)

  # var.explained is of the underlying eigen-decomposition. limma reports it for
  # both gene.selection modes but it is only a true variance fraction for
  # "common"; under "pairwise" the axes come from a non-Euclidean distance
  # matrix, so it is a rough guide to relative axis importance rather than a
  # percentage of variance. Labelled "scaling" for that reason.
  pct <- round(100 * mds$var.explained[mds$dim.plot])
  axis_label <- function(i) {
    if (length(pct) >= i && is.finite(pct[i])) {
      glue("{mds$axislabel} {mds$dim.plot[i]} ({pct[i]}% scaling)")
    } else {
      glue("{mds$axislabel} {mds$dim.plot[i]}")
    }
  }

  tibble::tibble(
    dim1 = mds$x,
    dim2 = mds$y,
    group = as.character(coldata[[colour_var]]),
    sample_label = as.character(labels)
  ) |>
    ggplot(aes(dim1, dim2, color = group)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = sample_label), size = 3, max.overlaps = 20) +
    xlab(axis_label(1)) +
    ylab(axis_label(2)) +
    labs(title = title, color = legend_title) +
    theme_bw()
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
sex_findings <- check_sex_coherence(meta, block_var = block_var)
# The marker scatter alongside the DE outputs. The sex_qc rule already wrote one
# from the same function, but a reader of this directory should not have to go
# looking for it in another to see which libraries the fit was built on.
invisible(plot_sex_inference(meta, out_dir, sex_findings))

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

## ---- Ordination across all samples -----------------------------------------
## On the voom log-CPM the model is actually fitted on -- v$E, not a separate
## cpm(log = TRUE) call. voom computes log2((count + 0.5) / (lib.size + 1) * 1e6)
## against the TMM norm factors, which is a different offset from the
## prior.count = 3 used before, so the ordination now sits on exactly the values
## the contrasts are estimated from rather than on a parallel transformation of
## the same counts. This is why the block runs after the fit: v does not exist
## until voom has been called.
##
## The voom precision weights are NOT used here. plotMDS has no weighted mode,
## so the ordination is unweighted even under quality_weights = per_sample; a
## down-weighted sample still gets equal say in the distances.

voom_logcpm <- v$E

# Points are labelled "<age group>-<condition>", so whichever variable a panel
# is coloured by, the label still carries the other. The two sex panels are the
# exception: they are read to find WHICH library sits on the wrong side, and
# neither age group nor region identifies an animal, so those use the sample
# name (animal + region) instead.
mds_labels <- paste(meta$age_group, meta$condition, sep = "-")
sex_labels <- as.character(meta$sample_name)

mds_titles <- c(condition = "brain region", age_group = "age group",
                line = "animal", inferred_sex = "sex inferred from expression",
                stated_sex = "sex as provided by the client")
mds_vars <- list(condition = mds_labels, age_group = mds_labels, line = mds_labels,
                 inferred_sex = sex_labels)

# stated_sex only earns a panel when the metadata actually states a sex. It is
# all NA when the column is blank, which would plot 27 points in one "NA" group
# and imply the labels were checked when they were never supplied.
if ("stated_sex" %in% colnames(meta) && any(!is.na(meta$stated_sex))) {
  mds_vars$stated_sex <- sex_labels
} else {
  message("No sex stated in the metadata; skipping the client-sex MDS panel")
}

# colour_var, not v: `v` is the voom EList and these panels now run after the
# fit, so reusing it as a loop variable would overwrite the object the contrast
# loop below reads v$E from.
mds_plots <- imap(mds_vars, \(lab, colour_var) mds_from_logcpm(
  voom_logcpm, meta, colour_var, lab,
  glue("MDS (voom log2-CPM) - all samples - {mds_titles[[colour_var]]}"),
  legend_title = mds_titles[[colour_var]]
))

pdf(file.path("results", "mds_plots_all_limma_voom.pdf"), width = 12, height = 6)
print(mds_plots)
dev.off()

# PNGs for the project report; the PDF above keeps every panel in vector form.
for (colour_var in intersect(c("condition", "age_group", "inferred_sex", "stated_sex"),
                             names(mds_plots))) {
  ggsave(file.path(out_dir, glue("mds_all_samples_{colour_var}.png")),
         mds_plots[[colour_var]], width = 7, height = 5, dpi = 200)
}

## ---- Ordination within one brain region ------------------------------------
## The all-samples panels above are dominated by region (condition), which can
## hide structure within a region. Here each region gets its own MDS over its
## own samples (9 per region: 3 age groups x 3 animals), coloured by age group
## and labelled by animal, so age-related structure is visible without the
## between-region variance swamping it. The pairwise gene selection and the
## axes are recomputed within the subset; the
## values themselves are the same v$E the global fit uses, so the panel stays on
## the same scale as everything else in this directory.

region_mds_plots <- list()
for (region in unique(as.character(meta$condition))) {
  idx <- which(meta$condition == region)
  region_meta <- droplevels(meta[idx, , drop = FALSE])
  region_voom_logcpm <- voom_logcpm[, idx, drop = FALSE]

  # Both sex panels are here for completeness, not because MDS resolves sex:
  # the markers are five genes out of ~19k, so they never dominate a leading
  # logFC dimension. These answer "is sex a major axis of variation here?"
  # (it is not, in any region) --
  # the per-library sex calls come from the marker scatter in inferred_sex.png.
  for (colour_var in intersect(names(mds_titles), names(mds_vars)) |>
       intersect(c("age_group", "inferred_sex", "stated_sex"))) {
    p <- mds_from_logcpm(
      region_voom_logcpm, region_meta, colour_var,
      as.character(region_meta$animal),
      glue("MDS (voom log2-CPM) - {region} only - {mds_titles[[colour_var]]}"),
      legend_title = mds_titles[[colour_var]]
    )
    region_mds_plots[[glue("{region}_{colour_var}")]] <- p
    ggsave(file.path(out_dir, glue("mds_region_{safe_filename(region)}_{colour_var}.png")),
           p, width = 7, height = 5, dpi = 200)
  }
}

pdf(file.path("results", "mds_plots_by_region_limma_voom.pdf"), width = 12, height = 6)
print(region_mds_plots)
dev.off()

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
  res_df <- annotate_gene_symbols(res_df, gene_ann)

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
  mds_plot <- tryCatch(
    mds_from_logcpm(v$E[, cmp_samples, drop = FALSE], side, "line_group",
                    meta[cmp_samples, "animal"], glue("MDS (voom log2-CPM) - {cmp_name}")),
    error = function(e) {
      warning(glue("MDS failed for {cmp_name}: {e$message}"))
      NULL
    }
  )
  if (!is.null(mds_plot)) {
    pdf(file.path(cmp_dir, glue("{base_name}_mds.pdf")), width = 8, height = 6)
    print(mds_plot)
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
