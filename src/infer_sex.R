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
  library(ggplot2)
  library(ggrepel)
  library(scales)
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

# The metadata `sex` column as a resolved word, or NA where it is blank or
# unrecognised. Defined once and used by both the coherence check and the plot,
# so the two can never disagree about what the client actually stated.
resolve_stated_sex <- function(sex) {
  s <- tolower(str_trim(as.character(sex)))
  dplyr::case_when(
    startsWith(s, "m") ~ "male",
    startsWith(s, "f") ~ "female",
    TRUE ~ NA_character_
  )
}

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
  if ("sex" %in% names(meta) && any(!is.na(resolve_stated_sex(meta$sex)))) {
    stated <- resolve_stated_sex(meta$sex)
    bad <- which(!is.na(stated) & meta$inferred_sex %in% c("male", "female") &
                   stated != meta$inferred_sex)
    if (length(bad) > 0) {
      # Roll the per-library hits up to the animal. A stated sex is a property
      # of the animal, not of one dissection, so "2 of 3 SOM1 libraries contradict
      # the label" is the actionable form -- it separates one odd library from a
      # whole animal being labelled wrong, and those call for different fixes.
      per_animal <- vapply(split(seq_len(nrow(meta)), meta$animal), function(ix) {
        hits <- intersect(ix, bad)
        if (length(hits) == 0) return(NA_character_)
        glue("{meta$animal[ix][1]} (stated {stated[ix][1]}): ",
             "{length(hits)}/{length(ix)} librar(ies) inferred ",
             "{paste(sort(unique(meta$inferred_sex[hits])), collapse = '/')}")
      }, character(1))
      msg <- glue("Metadata sex disagrees with inferred sex in ",
                  "{length(bad)} librar(ies) across {sum(!is.na(per_animal))} animal(s):\n  ",
                  "{paste(stats::na.omit(per_animal), collapse = '\n  ')}")
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
  i <- match(rownames(meta), sex_calls$sample)
  meta$inferred_sex <- sex_calls$inferred_sex[i]
  meta$y_cpm <- sex_calls$y_cpm[i]
  meta$xist_cpm <- sex_calls$xist_cpm[i]
  meta$stated_sex <- resolve_stated_sex(meta$sex)

  readr::write_csv(
    dplyr::mutate(sex_calls,
                  animal = meta$animal[match(sample, rownames(meta))],
                  condition = as.character(meta$condition[match(sample, rownames(meta))]),
                  stated_sex = meta$stated_sex[match(sample, rownames(meta))],
                  .after = sample) |>
      dplyr::mutate(agrees_with_metadata = dplyr::case_when(
        is.na(stated_sex) ~ NA,
        inferred_sex %in% c("male", "female") ~ inferred_sex == stated_sex,
        TRUE ~ NA
      )),
    file.path(out_dir, "inferred_sex.csv")
  )
  message(glue("Inferred sex: {paste(names(table(meta$inferred_sex)), ",
               "table(meta$inferred_sex), sep = '=', collapse = ', ')}"))
  meta
}

## ---- Findings table --------------------------------------------------------

# Machine-readable findings, always written -- header-only when clean. This is
# what lets the gate be a separate step from the report: infer_sex.R can always
# succeed and keep its CSV, plot and log, and a downstream rule decides whether
# a non-empty findings table should stop the pipeline. Deciding both here would
# mean the artifacts explaining the failure are deleted along with it.
write_sex_qc_findings <- function(findings, meta, out_dir) {
  label_for <- function(ids) {
    hit <- match(ids, rownames(meta))
    ifelse(is.na(hit), ids, as.character(meta$sample_name)[hit])
  }
  rows <- if (length(findings) == 0) {
    tibble::tibble(finding = character(), subject = character())
  } else {
    dplyr::bind_rows(lapply(names(findings), function(nm) {
      ids <- as.character(findings[[nm]])
      tibble::tibble(
        finding = nm,
        # incoherent_animals holds animal IDs; the rest hold sample identifiers.
        subject = if (nm == "incoherent_animals") ids else label_for(ids)
      )
    }))
  }
  path <- file.path(out_dir, "sex_qc_findings.csv")
  readr::write_csv(rows, path)
  message(glue("Wrote {path} ({nrow(rows)} finding row(s))"))
  invisible(rows)
}

## ---- Plot ------------------------------------------------------------------

# The two markers against each other, which is the only view that separates the
# three outcomes at a glance: clean libraries sit hard against one axis, while a
# contaminated one leaves the axis entirely and lands in the interior. A bar
# chart of the categorical call cannot show that, because the call is exactly
# what is in doubt.
#
# Both axes are pseudo-log: the marker signal spans zero to several hundred cpm
# and the informative structure is at the low end, where a linear axis would
# stack every clean library on the origin. pseudo_log is linear near zero, so
# the true zeros (y_cpm is exactly 0.00 for four libraries here) still plot.
plot_sex_inference <- function(meta, out_dir, findings = list()) {
  incoherent <- findings$incoherent_animals %||% character(0)

  df <- tibble::tibble(
    sample_name = as.character(meta$sample_name),
    animal = as.character(meta$animal),
    y_cpm = meta$y_cpm,
    xist_cpm = meta$xist_cpm,
    inferred_sex = factor(meta$inferred_sex,
                          levels = c("female", "male", "mixed", "unknown")),
    # Colour is what expression says, shape is what the client's sheet says, so
    # a disagreement is visible as a mismatch between the two channels without
    # having to cross-reference anything. The ring then calls it out explicitly.
    stated_sex = factor(dplyr::coalesce(resolve_stated_sex(meta$sex), "not stated"),
                        levels = c("female", "male", "not stated"))
  ) |>
    dplyr::mutate(
      disagrees = !is.na(stated_sex) & stated_sex != "not stated" &
        inferred_sex %in% c("male", "female") &
        as.character(stated_sex) != as.character(inferred_sex)
    )

  if (all(is.na(df$xist_cpm))) {
    warning("No Xist marker in the count matrix; skipping the sex-inference plot")
    return(invisible(NULL))
  }

  # Built here rather than inline in the subtitle: a glue() nested inside
  # another glue()'s braces re-parses the inner braces and dropped the paste(),
  # which silently printed only the first animal.
  incoherent_note <- if (length(incoherent) > 0) {
    paste0("; animal(s) ", paste(incoherent, collapse = ", "),
           " are internally inconsistent")
  } else ""

  ptrans <- scales::pseudo_log_trans(sigma = 0.1, base = 10)
  brk <- c(0, 0.1, 1, 10, 100, 1000)

  p <- ggplot2::ggplot(df, ggplot2::aes(y_cpm, xist_cpm)) +
    ggplot2::annotate("rect", xmin = Y_CPM_CUTOFF, xmax = Inf,
                      ymin = X_CPM_CUTOFF, ymax = Inf,
                      fill = "#b2182b", alpha = 0.06) +
    ggplot2::geom_vline(xintercept = Y_CPM_CUTOFF, linetype = "dashed",
                        colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = X_CPM_CUTOFF, linetype = "dashed",
                        colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = inferred_sex, shape = stated_sex),
                        size = 2.6, stroke = 1.1) +
    # Ringed separately rather than as another shape or colour: a disagreement
    # is a relation between the two channels already in use, not a third value
    # of either, and drawing it as one would imply it is an alternative call.
    ggplot2::geom_point(data = ~ dplyr::filter(.x, disagrees),
                        shape = 21, size = 5.5, stroke = 0.7,
                        colour = "grey15", fill = NA) +
    # seed fixes the repel layout: without it every rerun nudges the labels
    # somewhere new and the figure differs byte-for-byte on identical input,
    # which makes a pipeline artifact impossible to diff.
    ggrepel::geom_text_repel(ggplot2::aes(label = sample_name, colour = inferred_sex),
                             size = 2.5, max.overlaps = Inf, min.segment.length = 0,
                             box.padding = 0.5, point.padding = 0.35, force = 3,
                             seed = 1, segment.colour = "grey70",
                             segment.size = 0.25, show.legend = FALSE) +
    ggplot2::scale_x_continuous(transform = ptrans, breaks = brk,
                                labels = scales::label_number(drop0trailing = TRUE)) +
    ggplot2::scale_y_continuous(transform = ptrans, breaks = brk,
                                labels = scales::label_number(drop0trailing = TRUE)) +
    ggplot2::scale_colour_manual(
      values = c(female = "#2166ac", male = "#1a9850",
                 mixed = "#b2182b", unknown = "grey50")
    ) +
    ggplot2::scale_shape_manual(
      values = c(female = 16, male = 15, `not stated` = 4), drop = FALSE) +
    ggplot2::labs(
      title = "Sex inference per library",
      subtitle = glue(
        "Colour is inferred from expression, shape is the sex stated in the ",
        "metadata; a ring marks a library where the two disagree.\nDashed lines ",
        "are the calling cutoffs (Y {Y_CPM_CUTOFF} cpm, Xist {X_CPM_CUTOFF} cpm); ",
        "the shaded corner expresses BOTH marker sets, which is not a genotype. ",
        "\n{sum(df$inferred_sex == 'mixed', na.rm = TRUE)} of {nrow(df)} librar(ies) ",
        "are mixed, {sum(df$disagrees)} contradict the stated sex{incoherent_note}."
      ),
      x = "Y-linked genes (summed cpm)",
      y = "Xist (cpm)", colour = "Inferred sex", shape = "Stated sex"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 7.5, colour = "grey30"),
                   legend.position = "right")

  png_path <- file.path(out_dir, "inferred_sex.png")
  ggplot2::ggsave(png_path, p, width = 9, height = 6, dpi = 200)
  ggplot2::ggsave(file.path(out_dir, "inferred_sex.pdf"), p, width = 9, height = 6)
  message(glue("Wrote {png_path}"))
  invisible(p)
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
  plot_sex_inference(meta, out_dir, findings)
  write_sex_qc_findings(findings, meta, out_dir)

  message(glue("Wrote {file.path(out_dir, 'inferred_sex.csv')}"))
  if (length(findings) > 0) {
    message(glue("QC findings: {paste(names(findings), collapse = ', ')}"))
    # Echo the whole table into the log as well as the CSV: --strict exits
    # non-zero, and a workflow runner deletes a failed job's output files.
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
