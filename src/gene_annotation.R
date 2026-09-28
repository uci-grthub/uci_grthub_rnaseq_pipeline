#!/usr/bin/env Rscript
# Gene symbols and biotypes, taken from the same GTF the counts were made
# against, via an ensembldb EnsDb.
#
# Why not org.*.eg.db. It is keyed on Entrez records, so it has no symbol for
# lncRNA, TEC or pseudogene entries -- and it also misses real protein-coding
# genes (ENSMUSG00000021745 = Ptprg is one of 99 it fails on in GENCODE vM24).
# On this project it left 2467 of 19254 genes (12.8%) unlabelled, and its
# one-Ensembl-to-many-Entrez mappings produced a "'select()' returned 1:many
# mapping" warning on every contrast.
#
# Why not a TxDb. makeTxDbFromGFF() keeps transcript structure and DISCARDS the
# gene_name attribute: a TxDb built from this GTF exposes GENEID, TXNAME, TXTYPE
# and exon/CDS coordinates, and no SYMBOL or GENENAME column at all. It cannot
# answer this question. An EnsDb can -- its schema carries GENENAME and
# GENEBIOTYPE.
#
# Why the GTF. It is the file featureCounts counted against, so the mapping is
# complete (100% of the previously-missing genes have a gene_name), matched to
# the same annotation release, and strictly one-to-one.

suppressPackageStartupMessages({
  library(glue)
})

## ---- Build -----------------------------------------------------------------

# Builds the EnsDb SQLite from a GTF. `organism` is the ensembldb spelling with
# an underscore ("Mus_musculus"), `genome_version` a bare build name ("GRCm38"),
# `version` the annotation release as a number.
build_gene_ensdb <- function(gtf_file, sqlite_out, organism, genome_version, version) {
  suppressPackageStartupMessages({
    library(ensembldb)
    library(rtracklayer)
  })

  if (!file.exists(gtf_file)) {
    stop(glue("GTF not found: {gtf_file}"))
  }
  message(glue("Reading {gtf_file}"))
  gr <- rtracklayer::import(gtf_file, format = "gtf")

  # ensembldb expects Ensembl's attribute names. GENCODE writes gene_type and
  # transcript_type instead, and ensembldb does not translate: handed a GENCODE
  # GTF unrenamed it warns "I'm missing column(s): 'entrezid','gene_biotype'"
  # and stores every biotype as NA. Ensembl-sourced GTFs (the human reference
  # here) already use the expected names, so this is a no-op for them.
  mc <- names(S4Vectors::mcols(gr))
  if (!"gene_biotype" %in% mc && "gene_type" %in% mc) {
    S4Vectors::mcols(gr)$gene_biotype <- S4Vectors::mcols(gr)$gene_type
    message("Renamed GENCODE 'gene_type' to 'gene_biotype' for ensembldb")
  }
  if (!"tx_biotype" %in% mc && "transcript_type" %in% mc) {
    S4Vectors::mcols(gr)$tx_biotype <- S4Vectors::mcols(gr)$transcript_type
    message("Renamed GENCODE 'transcript_type' to 'tx_biotype' for ensembldb")
  }

  if (file.exists(sqlite_out)) file.remove(sqlite_out)

  # ensembldb tries to fetch sequence lengths from Ensembl over the network.
  # Compute nodes here have no outbound access, and it degrades gracefully, so
  # muffle only those two warnings and let anything else through.
  known_noise <- "Could not determine length for all seqnames|Unable to retrieve sequence lengths"
  withCallingHandlers(
    ensembldb::ensDbFromGRanges(
      gr, outfile = sqlite_out, organism = organism,
      genomeVersion = genome_version, version = version
    ),
    warning = function(w) {
      if (grepl(known_noise, conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )

  edb <- ensembldb::EnsDb(sqlite_out)
  ann <- ensembldb::genes(edb, return.type = "data.frame")
  n_biotype_na <- sum(is.na(ann$gene_biotype))
  message(glue("Wrote {sqlite_out}: {nrow(ann)} genes, ",
               "{sum(is.na(ann$gene_name) | !nzchar(ann$gene_name))} without a name, ",
               "{n_biotype_na} without a biotype"))
  if (n_biotype_na == nrow(ann)) {
    warning("Every gene_biotype is NA; the GTF's biotype attribute was not recognised")
  }
  invisible(sqlite_out)
}

## ---- Annotate ---------------------------------------------------------------

# Reads the gene table once. Kept separate so a caller annotating several frames
# pays the SQLite read only once.
read_gene_annotation <- function(ensdb_path) {
  suppressPackageStartupMessages(library(ensembldb))
  if (!file.exists(ensdb_path)) {
    stop(glue("Gene annotation database not found: {ensdb_path}"))
  }
  ann <- ensembldb::genes(ensembldb::EnsDb(ensdb_path), return.type = "data.frame")
  ann[, intersect(c("gene_id", "gene_name", "gene_biotype"), names(ann)), drop = FALSE]
}

# Attaches gene_symbol and gene_biotype to a results frame keyed on `gene`, and
# moves them to the front. `ann` may be a path or the frame read_gene_annotation
# returned.
#
# Matching is exact first, then version-stripped on both sides. Both are needed:
# the mouse GENCODE GTF carries versioned gene IDs (ENSMUSG00000102693.1) while
# the human Ensembl GTF does not (ENSG00000223972), and the count matrix
# inherits whichever form its own GTF used.
annotate_gene_symbols <- function(df, ann) {
  if (is.character(ann)) ann <- read_gene_annotation(ann)

  idx <- match(df$gene, ann$gene_id)
  unmatched <- is.na(idx)
  if (any(unmatched)) {
    strip <- function(x) sub("\\.\\d+$", "", x)
    idx[unmatched] <- match(strip(df$gene[unmatched]), strip(ann$gene_id))
  }

  df$gene_symbol <- ann$gene_name[idx]
  df$gene_symbol[!is.na(df$gene_symbol) & !nzchar(df$gene_symbol)] <- NA_character_
  df$gene_biotype <- if ("gene_biotype" %in% names(ann)) ann$gene_biotype[idx] else NA_character_

  n_missing <- sum(is.na(df$gene_symbol))
  if (n_missing > 0) {
    warning(glue("{n_missing} of {nrow(df)} genes have no symbol in the annotation database; ",
                 "check that it was built from the same GTF as the counts"))
  }

  lead <- c("gene", "gene_symbol", "gene_biotype")
  df[, c(lead, setdiff(colnames(df), lead)), drop = FALSE]
}

## ---- CLI --------------------------------------------------------------------

# True only when this file is the script Rscript was pointed at, so sourcing it
# from a DE script does not trigger the build.
if (identical(basename(sub("^--file=", "",
                           commandArgs(trailingOnly = FALSE)[
                             grep("^--file=", commandArgs(trailingOnly = FALSE))][1])),
              "gene_annotation.R")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 5) {
    stop("Usage: Rscript src/gene_annotation.R <gtf> <out.sqlite> <organism> <genome> <version>")
  }
  dir.create(dirname(args[2]), showWarnings = FALSE, recursive = TRUE)
  build_gene_ensdb(args[1], args[2], args[3], args[4], as.numeric(args[5]))
}
