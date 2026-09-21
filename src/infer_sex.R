#!/usr/bin/env Rscript
# Sex inference from expression, and the sample-tracking checks built on it.
#
# Dual purpose. Sourced, it provides the marker definitions, infer_sex() and
# check_sex_coherence() to limma_voom_analysis.R. Run directly, it is a
# standalone QC gate:
#
#   Rscript src/infer_sex.R counts.txt metadata.csv [out_dir] [--strict]
#
# with --strict exiting non-zero when anything is flagged, so it can gate a
# pipeline before any model is fitted.
#
# Why infer sex at all. The `sex` column in metadata is routinely empty, and
# even when filled it records what was intended rather than what was
# sequenced -- so it cannot catch a sample swap. Expression can: Y-linked genes
# and Xist are near-binary and mutually exclusive, and an animal whose brain
# regions disagree about which of them is expressed has a labelling problem that
# no amount of modelling will fix.

suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
  library(glue)
  library(readr)
  library(stringr)
  library(tibble)
})

## ---- Markers and cutoffs ---------------------------------------------------

# Y-linked genes and Xist are near-binary. In this project clean females carry
# 0.0-0.3 cpm summed over the Y genes against 81-132 cpm in males, while Xist
# runs 229-655 cpm in females against 0.0-3.9 in males. A library showing BOTH
# marker sets is not a valid genotype and is called "mixed" rather than forced
# into one or the other.
SEX_MARKERS <- list(
  ENSMUSG = list(y = c("ENSMUSG00000069045",   # Ddx3y
                       "ENSMUSG00000068457",   # Uty
                       "ENSMUSG00000069049",   # Eif2s3y
                       "ENSMUSG00000056673"),  # Kdm5d
                 x = "ENSMUSG00000086503"),    # Xist
  ENSG = list(y = c("ENSG00000067048",         # DDX3Y
                    "ENSG00000183878",         # UTY
                    "ENSG00000198692",         # EIF1AY
                    "ENSG00000012817"),        # KDM5D
              x = "ENSG00000229807")           # XIST
)

# Cutoffs sit in the empty space between the two populations, which spans two
# orders of magnitude in both directions, so their exact value is not critical.
Y_CPM_CUTOFF <- 5
X_CPM_CUTOFF <- 20

## ---- Inference -------------------------------------------------------------

# Takes the RAW, UNFILTERED count matrix. In a single-sex experiment the Y genes
# would not survive filterByExpr, and their absence is precisely the signal, so
# this must run before any expression filter.
infer_sex <- function(counts) {
  ids <- sub("\\.\\d+$", "", rownames(counts))
  blank <- tibble::tibble(sample = colnames(counts), y_cpm = NA_real_,
                          xist_cpm = NA_real_, inferred_sex = "unknown")
  prefix <- names(SEX_MARKERS)[
    vapply(names(SEX_MARKERS), \(p) any(startsWith(ids, p)), logical(1))]
  if (length(prefix) != 1) {
    warning("No sex markers defined for these gene IDs; inferred_sex = unknown")
    return(blank)
  }
  mk <- SEX_MARKERS[[prefix]]

  cpm_all <- cpm(counts)
  y_rows <- which(ids %in% mk$y)
  x_rows <- which(ids %in% mk$x)
  if (length(y_rows) == 0) {
    warning("No Y-linked marker genes in the count matrix; inferred_sex = unknown")
    return(blank)
  }
  y_cpm <- colSums(cpm_all[y_rows, , drop = FALSE])
  xist_cpm <- if (length(x_rows) > 0) colSums(cpm_all[x_rows, , drop = FALSE]) else NA_real_

  tibble::tibble(
    sample = colnames(counts),
    y_cpm = round(y_cpm, 2),
    xist_cpm = round(xist_cpm, 2),
    inferred_sex = dplyr::case_when(
      y_cpm > Y_CPM_CUTOFF & !is.na(xist_cpm) & xist_cpm > X_CPM_CUTOFF ~ "mixed",
      y_cpm > Y_CPM_CUTOFF ~ "male",
      !is.na(xist_cpm) & xist_cpm > X_CPM_CUTOFF ~ "female",
      TRUE ~ "unknown"
    )
  )
}

## ---- Coherence checks ------------------------------------------------------

# `meta` must carry `inferred_sex`, `y_cpm`, `animal` and `sample_name`, with
# rownames matching the count-matrix columns. Emits a warning per finding and
# returns them, so a caller can decide whether to stop.
check_sex_coherence <- function(meta, block_var = "animal") {
  findings <- list()
  # `%in%` rather than `==` throughout: an inferred_sex of NA (a count column
  # that never matched a metadata row) would otherwise index rownames() with NA
  # and put literal NAs in the finding.
  names_for <- function(ids) meta$sample_name[match(ids, rownames(meta))]

  # A library carrying both marker sets is a mixture, not a genotype.
  mixed <- rownames(meta)[meta$inferred_sex %in% "mixed"]
  if (length(mixed) > 0) {
    msg <- glue("{length(mixed)} librar(ies) express both Y-linked genes and Xist, ",
                "indicating contaminated or pooled material: ",
                "{paste(names_for(mixed), collapse = ', ')}")
    warning(msg)
    findings$mixed <- mixed
  }

  # "unknown" means the inference itself produced nothing for that library --
  # no markers in the matrix, an ambiguous ID scheme, or a library clearing
  # neither cutoff. It has to be a finding in its own right: when the marker
  # lookup fails wholesale EVERY call is "unknown", which makes all the checks
  # below pass vacuously, and a --strict gate would then exit 0 while reporting
  # that nothing is wrong. Silence from a check that never ran is not a pass.
  unknown <- rownames(meta)[meta$inferred_sex %in% "unknown" | is.na(meta$inferred_sex)]
  if (length(unknown) > 0) {
    msg <- glue("{length(unknown)} of {nrow(meta)} librar(ies) could not be assigned a ",
                "sex; the inference produced no call for them, so the checks below say ",
                "nothing about them: {paste(names_for(unknown), collapse = ', ')}")
    warning(msg)
    findings$unknown <- unknown
  }

  # Sex must be constant within an animal. If it is not, the animal labels are
  # wrong, and blocking on them -- which assumes one block is one animal -- is
  # modelling something that does not exist.
  #
  # Test on Y-bearing material rather than on the categorical call. A "mixed"
  # library still has a dominant genotype, so collapsing it to a single label
  # (or discarding it) hides real contradictions: a contaminated library whose Y
  # genes sit at full male level, in an animal whose other regions carry zero Y,
  # is a swap regardless of how much Xist rode along with it.
  per_animal <- tapply(meta$y_cpm > Y_CPM_CUTOFF, meta$animal,
                       \(has_y) length(unique(has_y[!is.na(has_y)])))
  incoherent <- names(which(per_animal > 1))
  if (length(incoherent) > 0) {
    msg <- glue("Inferred sex is not constant within animal(s): ",
                "{paste(incoherent, collapse = ', ')}. Blocking on '{block_var}' assumes each ",
                "block is one animal, so check these for sample swaps before trusting the fit.")
    warning(msg)
    findings$incoherent_animals <- incoherent
  }

  # If the metadata states a sex, the inference should agree with it. Compare
  # the resolved words, not first letters: substr(inferred_sex, 1, 1) maps
  # "mixed" onto "m", which would report every contaminated library as a
  # male-vs-stated disagreement and hide the contamination behind the wrong
  # finding. Only a confident "male"/"female" call can disagree with metadata;
  # "mixed" and "unknown" are reported by their own checks above.
  if ("sex" %in% names(meta) && any(nzchar(as.character(meta$sex)) & !is.na(meta$sex))) {
    stated_raw <- tolower(str_trim(as.character(meta$sex)))
    stated <- dplyr::case_when(
      startsWith(stated_raw, "m") ~ "male",
      startsWith(stated_raw, "f") ~ "female",
      TRUE ~ NA_character_
    )
    bad <- which(!is.na(stated) & meta$inferred_sex %in% c("male", "female") &
                   stated != meta$inferred_sex)
    if (length(bad) > 0) {
      msg <- glue("Metadata sex disagrees with inferred sex for: ",
                  "{paste(meta$sample_name[bad], collapse = ', ')}")
      warning(msg)
      findings$metadata_disagreement <- meta$sample_name[bad]
    }
  }

  findings
}

## ---- Shared loaders --------------------------------------------------------
## Here rather than duplicated in each caller: the barcode join in particular is
## fiddly enough that two copies would drift.

load_counts <- function(counts_file) {
  counts <- read.table(counts_file, header = TRUE, row.names = 1)
  count_matrix <- counts[, 6:ncol(counts)] # cols 1-5 are featureCounts annotation

  # read.table() mangles the BAM paths featureCounts writes as column headers.
  colnames(count_matrix) <- colnames(count_matrix) |>
    str_remove("^output\\.hisat2_alignment\\.") |>
    str_remove("_align_sorted_markdup\\.bam$") |>
    str_replace_all("\\.", "-")
  count_matrix
}

# Returns list(meta, counts) with the count columns reordered to match the
# metadata rows, so callers can rely on the two lining up.
load_sample_metadata <- function(meta_file, count_matrix) {
  meta_raw <- read.csv(meta_file) |> janitor::clean_names()

  # clean_names() renders "i5barcode_NovaSeqV1.5" differently across janitor
  # versions (2.2.1 splits the camel case), so match the column by prefix.
  i5_col <- grep("^i5barcode", names(meta_raw), value = TRUE)
  if (length(i5_col) != 1) {
    stop(glue("Expected one i5barcode column in {meta_file}, ",
              "found: {paste(i5_col, collapse = ', ')}"))
  }

  missing_cols <- setdiff(c("sample_name", "age_group", "replicate", "condition"),
                          names(meta_raw))
  if (length(missing_cols) > 0) {
    stop(glue("{meta_file} is missing column(s): {paste(missing_cols, collapse = ', ')}"))
  }

  # Samples are matched to count-matrix columns on the i7-i5 barcode pair, which
  # both the metadata and the FASTQ-derived column names carry.
  sample_map <- tibble::tibble(
    sample_col = colnames(count_matrix),
    index_pair = str_extract(colnames(count_matrix), "[ACGT]+-[ACGT]+")
  )

  # The animal. `mouse_id` states it explicitly; metadata files predating that
  # column encode the same thing as age group + replicate, so fall back to it
  # rather than failing. `line` is kept as an alias so the PCA panels and any
  # config naming it keep working.
  has_mouse_id <- "mouse_id" %in% names(meta_raw) && !all(is.na(meta_raw$mouse_id))
  animal_source <- if (has_mouse_id) "mouse_id" else "derived (age_group + replicate)"

  meta <- meta_raw |>
    dplyr::mutate(
      index_pair = paste0(i7barcode, "-", .data[[i5_col]]),
      animal = if (has_mouse_id) as.character(mouse_id) else paste0(age_group, replicate),
      line = animal,
      condition = factor(condition)
    ) |>
    dplyr::inner_join(sample_map, by = "index_pair") |>
    tibble::column_to_rownames("sample_col")

  if (nrow(meta) == 0) {
    stop(glue("No samples in {meta_file} matched count matrix columns"))
  }

  message(glue("Animal identity from {animal_source}: ",
               "{length(unique(meta$animal))} animals over {nrow(meta)} samples"))

  # When both are available they must agree; a disagreement means the metadata
  # is internally inconsistent and every blocked fit would be wrong.
  if (has_mouse_id) {
    derived <- paste0(meta$age_group, meta$replicate)
    if (!identical(as.character(meta$animal), derived)) {
      warning(glue("mouse_id disagrees with age_group + replicate for: ",
                   "{paste(unique(meta$animal[meta$animal != derived]), collapse = ', ')}"))
    }
  }

  list(meta = meta, counts = count_matrix[, rownames(meta), drop = FALSE])
}

# Attaches the calls to `meta` and writes the table. Returns the annotated meta.
annotate_and_write_sex <- function(meta, count_matrix, out_dir) {
  sex_calls <- infer_sex(count_matrix)
  meta$inferred_sex <- sex_calls$inferred_sex[match(rownames(meta), sex_calls$sample)]
  meta$y_cpm <- sex_calls$y_cpm[match(rownames(meta), sex_calls$sample)]

  readr::write_csv(
    dplyr::mutate(sex_calls,
                  animal = meta$animal[match(sample, rownames(meta))],
                  condition = as.character(meta$condition[match(sample, rownames(meta))]),
                  .after = sample),
    file.path(out_dir, "inferred_sex.csv")
  )
  message(glue("Inferred sex: {paste(names(table(meta$inferred_sex)), ",
               "table(meta$inferred_sex), sep = '=', collapse = ', ')}"))
  meta
}

## ---- CLI -------------------------------------------------------------------

# True only when this file is the script Rscript was pointed at, so sourcing it
# from limma_voom_analysis.R does not trigger the CLI.
invoked_directly <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  length(f) == 1 && identical(basename(f), "infer_sex.R")
}

if (invoked_directly()) {
  suppressPackageStartupMessages(library(janitor))
  args <- commandArgs(trailingOnly = TRUE)
  strict <- "--strict" %in% args
  args <- args[args != "--strict"]

  if (length(args) < 2) {
    stop("Usage: Rscript src/infer_sex.R counts.txt metadata.csv [out_dir] [--strict]")
  }
  counts_file <- args[1]
  meta_file <- args[2]
  out_dir <- if (length(args) >= 3) args[3] else "."
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  loaded <- load_sample_metadata(meta_file, load_counts(counts_file))
  meta <- annotate_and_write_sex(loaded$meta, loaded$counts, out_dir)
  findings <- check_sex_coherence(meta, block_var = "animal")

  message(glue("Wrote {file.path(out_dir, 'inferred_sex.csv')}"))
  if (length(findings) > 0) {
    message(glue("QC findings: {paste(names(findings), collapse = ', ')}"))
    # Echo the whole table into the log. Under --strict this run exits non-zero,
    # and a workflow runner deletes a failed job's output files, so the CSV that
    # would otherwise carry the evidence is gone by the time anyone reads this.
    message("\nPer-library calls:")
    message(paste(capture.output(print(
      data.frame(sample_name = meta$sample_name, animal = meta$animal,
                 condition = as.character(meta$condition), y_cpm = meta$y_cpm,
                 inferred_sex = meta$inferred_sex),
      row.names = FALSE
    )), collapse = "\n"))
    if (strict) {
      quit(status = 1)
    }
  } else {
    message("No sex-tracking problems found.")
  }
}
