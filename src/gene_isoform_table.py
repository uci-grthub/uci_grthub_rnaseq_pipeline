#!/usr/bin/env python3
"""
Per-library isoform table for one gene from the tximport matrices.

Reads transcript_counts.csv, transcript_tpm.csv and tpm_salmon.csv written by
tximport_tpm.R, and writes one row per library: the metadata animal, condition
and region, Salmon estimated count and TPM for each transcript of the gene, the
gene-level TPM, and the featureCounts count of the same gene for comparison.

USAGE:
    python gene_isoform_table.py --gene HTR2A --tpm-dir output/tpm/macaque_refseq \
        --gene-counts output/counts/macaque/gene_counts.csv \
        --metadata metadata/metadata.csv --output HTR2A_isoforms.csv

The gene is matched on the `gene` column of the transcript tables (RefSeq GTFs
use the symbol as gene_id) or on `symbol` when present (Ensembl), and on
`gene_id` or `symbol` in gene_counts.csv.
"""

import argparse
import csv
import os
import sys

# Metadata columns carried into the table, renamed for the client
META_COLUMNS = (('nhp', 'animal'), ('condition', 'condition'), ('region', 'region'))


def read_rows(path):
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh)
        return reader.fieldnames, list(reader)


def gene_matches(row, gene):
    return row.get('gene') == gene or row.get('symbol') == gene


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--gene', required=True)
    parser.add_argument('--tpm-dir', required=True)
    parser.add_argument('--gene-counts', required=True)
    parser.add_argument('--metadata', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    header, counts = read_rows(os.path.join(args.tpm_dir, 'transcript_counts.csv'))
    _, tpms = read_rows(os.path.join(args.tpm_dir, 'transcript_tpm.csv'))
    _, gene_tpms = read_rows(os.path.join(args.tpm_dir, 'tpm_salmon.csv'))
    _, fc_rows = read_rows(args.gene_counts)

    libraries = [h for h in header if h not in ('transcript', 'gene', 'symbol')]
    tx_counts = {r['transcript']: r for r in counts if gene_matches(r, args.gene)}
    tx_tpms = {r['transcript']: r for r in tpms if gene_matches(r, args.gene)}
    transcripts = sorted(tx_counts)
    if not transcripts:
        sys.exit(f"No transcripts of {args.gene} in {args.tpm_dir}/transcript_counts.csv")

    gene_tpm = next((r for r in gene_tpms if gene_matches(r, args.gene)), None)
    fc = next((r for r in fc_rows if args.gene in (r.get('gene_id'), r.get('symbol'))), None)

    # Libraries in metadata order, which groups them by animal
    meta = {}
    with open(args.metadata, newline='') as fh:
        for row in csv.DictReader(fh):
            meta.setdefault(row.get('library') or row['sample'], row)
    order = {lib: i for i, lib in enumerate(meta)}
    libraries.sort(key=lambda lib: order.get(lib, len(order)))
    meta_cols = [(src, dst) for src, dst in META_COLUMNS
                 if any(src in row for row in meta.values())]

    out_header = ['library'] + [dst for _, dst in meta_cols]
    out_header += [f'{tx}_count' for tx in transcripts]
    out_header += [f'{tx}_TPM' for tx in transcripts]
    out_header += [f'{args.gene}_gene_TPM', f'{args.gene}_featureCounts']

    with open(args.output, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(out_header)
        for lib in libraries:
            row = meta.get(lib, {})
            writer.writerow(
                [lib]
                + [row.get(src, '') for src, _ in meta_cols]
                + [f"{float(tx_counts[tx][lib]):.1f}" for tx in transcripts]
                + [f"{float(tx_tpms[tx][lib]):.2f}" for tx in transcripts]
                + [f"{float(gene_tpm[lib]):.2f}" if gene_tpm else '']
                + [fc.get(lib, '') if fc else '']
            )
    print(f"Wrote {len(libraries)} libraries x {len(transcripts)} {args.gene} transcripts to {args.output}")


if __name__ == '__main__':
    main()
