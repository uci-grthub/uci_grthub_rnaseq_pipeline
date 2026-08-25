#!/usr/bin/env Rscript
# Shared transcript-to-gene mapping, sourced by the tximport entry points.
#
# Kept in one place because the GTF parse below has two non-obvious details
# (the quote="" argument and the version stripping) that must stay identical
# across every tximport call, or the gene- and transcript-level matrices stop
# agreeing with each other.

suppressPackageStartupMessages({
  library(glue)
})

# Returns a data.frame(TXNAME, GENEID) keyed on version-stripped transcript IDs.
# Falls back to a transcript-as-gene mapping read from `fallback_quant_file`
# when the GTF is missing or unparseable.
build_tx2gene <- function(gtf_file, fallback_quant_file) {
  tx2gene <- NULL

  if (file.exists(gtf_file)) {
    message(glue("Reading transcript-to-gene mapping from {gtf_file}"))

    tx2gene <- tryCatch({
      # quote="" is required: GTF column 9 contains double-quotes (e.g. gene_id "ENSG...").
      # With the default quote='"', read.delim treats those as field quoting, swallows
      # tab/newline boundaries across the file, and collapses the mapping to a single NA row.
      gtf <- read.delim(gtf_file, header = FALSE, comment.char = "#", quote = "")

      attributes <- gtf$V9

      tx_ids <- stringr::str_extract(attributes, 'transcript_id "([^"]+)"')
      tx_ids <- stringr::str_remove(tx_ids, 'transcript_id "')
      tx_ids <- stringr::str_remove(tx_ids, '"')

      gene_ids <- stringr::str_extract(attributes, 'gene_id "([^"]+)"')
      gene_ids <- stringr::str_remove(gene_ids, 'gene_id "')
      gene_ids <- stringr::str_remove(gene_ids, '"')

      mapping <- data.frame(
        TXNAME = tx_ids,
        GENEID = gene_ids,
        stringsAsFactors = FALSE
      )
      # Salmon index here is built without transcript versions, but the GENCODE GTF
      # carries them (ENSMUST00000193812.1). Strip versions on this side: older
      # tximport applies ignoreTxVersion only to the quantification IDs, so the
      # mismatch survives otherwise.
      mapping$TXNAME <- sub("\\..*", "", mapping$TXNAME)
      mapping <- mapping[!is.na(mapping$TXNAME) & !is.na(mapping$GENEID), ]
      mapping <- mapping[!duplicated(mapping$TXNAME), ]
      mapping
    }, error = function(e) {
      warning(glue("Failed to parse GTF: {e$message}. Using transcript IDs as gene IDs."))
      NULL
    })
  } else {
    warning(glue("GTF file not found at {gtf_file}"))
  }

  if (is.null(tx2gene)) {
    message("Creating transcript-to-gene mapping from quantification files")
    sample_quant <- read.delim(fallback_quant_file)
    tx2gene <- data.frame(
      TXNAME = sample_quant$Name,
      GENEID = sample_quant$Name, # Use transcript ID as gene ID
      stringsAsFactors = FALSE
    )
  }

  message(glue("Transcript-to-gene mapping has {nrow(tx2gene)} entries"))
  tx2gene
}

# Locate the per-sample quant.sf files under a staged salmon directory and name
# them by sample, i.e. <dir>/<SAMPLE>_salmon_quant/<SAMPLE>_quant.sf.
find_quant_files <- function(salmon_dir) {
  quant_files <- list.files(
    path = salmon_dir,
    pattern = "quant.sf",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(quant_files) == 0) {
    stop(glue("No quant.sf files found in {salmon_dir}"))
  }

  sample_names <- sapply(quant_files, function(path) {
    parts <- strsplit(path, "/")[[1]]
    dir_name <- parts[length(parts) - 1]
    gsub("_salmon_quant", "", dir_name)
  }, USE.NAMES = FALSE)

  names(quant_files) <- sample_names
  message(glue("Found {length(quant_files)} Salmon quantification files"))
  quant_files
}
