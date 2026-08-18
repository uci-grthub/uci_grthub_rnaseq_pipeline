#!/usr/bin/env python3
"""Assemble the metadata and file manifests needed for an NCBI GEO/SRA submission.

GEO wants a sample sheet describing each library plus md5 checksums for every
raw and processed file it will host; SRA wants a parallel sheet describing the
sequencing run itself. This builds both from the project metadata, the FASTQ
directory, and the processed outputs, and reuses a pre-existing md5sums.txt in
the FASTQ directory rather than re-hashing tens of gigabytes.

Usage:
    prepare_ncbi_submission.py --metadata metadata/metadata.csv \
        --fastq-dir data/FASTQ --processed-dir output/counts/mouse \
        --tpm-dir output/tpm/mouse --species mouse --output-dir output/ncbi_submission/mouse
"""

import argparse
import csv
import hashlib
import os
import sys

# GEO/SRA controlled vocabulary for this library prep. Every sample in this
# project shares them, so they are constants rather than metadata columns.
LIBRARY_STRATEGY = "RNA-Seq"
LIBRARY_SOURCE = "TRANSCRIPTOMIC"
LIBRARY_SELECTION = "cDNA"
LIBRARY_LAYOUT = "paired"
PLATFORM = "ILLUMINA"
FILETYPE = "fastq"

ORGANISM_BY_SPECIES = {
    "mouse": "Mus musculus",
    "human": "Homo sapiens",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, help="Project metadata CSV")
    parser.add_argument("--fastq-dir", required=True, help="Directory holding raw FASTQ files")
    parser.add_argument("--processed-dir", default=None, help="Directory with the count matrix")
    parser.add_argument("--tpm-dir", default=None, help="Directory with the TPM matrix")
    parser.add_argument("--species", default="mouse", help="Species key used for this submission")
    parser.add_argument("--output-dir", required=True, help="Where to write the submission sheets")
    parser.add_argument("--instrument-model", default="Illumina NovaSeq 6000")
    parser.add_argument("--tissue", default="colon")
    parser.add_argument("--genome-build", default="", help="Genome build the processed files use")
    parser.add_argument("--annotation", default="", help="Annotation file the processed files use")
    return parser.parse_args()


def load_metadata(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"No rows found in metadata file {path}")
    return rows


def load_known_md5s(fastq_dir):
    """Read a checksum file written alongside the FASTQ files, if there is one."""
    md5_path = os.path.join(fastq_dir, "md5sums.txt")
    known = {}
    if not os.path.isfile(md5_path):
        return known
    with open(md5_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) != 2:
                continue
            checksum, name = parts
            known[os.path.basename(name)] = checksum
    return known


def md5_of(path, chunk_size=1 << 20):
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def find_read_files(fastq_dir, sample):
    """Return (R1, R2) basenames for a sample across the layouts the pipeline accepts."""
    candidates = [
        (f"{sample}-R1.fastq.gz", f"{sample}-R2.fastq.gz"),
        (f"{sample}_R1.fq.gz", f"{sample}_R2.fq.gz"),
        (f"{sample}_r1.fq.gz", f"{sample}_r2.fq.gz"),
        (f"{sample}-READ1-Sequences.txt.gz", f"{sample}-READ2-Sequences.txt.gz"),
    ]
    for r1, r2 in candidates:
        if os.path.isfile(os.path.join(fastq_dir, r1)):
            return r1, r2
        nested = os.path.join(fastq_dir, sample)
        if os.path.isfile(os.path.join(nested, r1)):
            return os.path.join(sample, r1), os.path.join(sample, r2)
    return "", ""


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    rows = load_metadata(args.metadata)
    species_rows = [r for r in rows if (r.get("species") or args.species).strip() == args.species]
    if not species_rows:
        species_rows = rows

    organism = ORGANISM_BY_SPECIES.get(args.species, args.species)
    known_md5s = load_known_md5s(args.fastq_dir)

    processed_files = []
    for directory in (args.processed_dir, args.tpm_dir):
        if not directory or not os.path.isdir(directory):
            continue
        for name in sorted(os.listdir(directory)):
            # RDS files are R-specific; GEO expects portable text/tabular formats
            if name.endswith(".csv") or name.endswith(".txt"):
                processed_files.append(os.path.join(directory, name))

    processed_names = [os.path.basename(p) for p in processed_files]
    shared_processed = ";".join(processed_names)

    geo_path = os.path.join(args.output_dir, "geo_samples.csv")
    sra_path = os.path.join(args.output_dir, "sra_metadata.csv")
    md5_path = os.path.join(args.output_dir, "md5sums.txt")

    missing_md5 = []
    with open(geo_path, "w", newline="") as geo_fh, \
            open(sra_path, "w", newline="") as sra_fh, \
            open(md5_path, "w") as md5_fh:
        geo_writer = csv.writer(geo_fh)
        geo_writer.writerow([
            "library name", "title", "organism", "tissue", "strain", "sex",
            "condition", "molecule", "single or paired-end", "instrument model",
            "description", "processed data file", "raw file 1", "raw file 2",
        ])
        sra_writer = csv.writer(sra_fh)
        sra_writer.writerow([
            "library_ID", "title", "library_strategy", "library_source",
            "library_selection", "library_layout", "platform", "instrument_model",
            "design_description", "filetype", "filename", "filename2",
        ])

        for row in species_rows:
            sample = (row.get("sample") or "").strip()
            if not sample:
                continue
            label = (row.get("sample_id") or "").strip() or sample
            condition = (row.get("condition") or "").strip() or "unspecified"
            strain = (row.get("strain") or "").strip()
            sex = (row.get("sex") or "").strip() or "not collected"
            r1, r2 = find_read_files(args.fastq_dir, sample)
            title = f"{organism} {args.tissue} RNA-seq, {condition}, {label}"

            geo_writer.writerow([
                label, title, organism, args.tissue, strain, sex, condition,
                "polyA RNA", "paired-end", args.instrument_model,
                f"Stranded paired-end mRNA-seq of {args.tissue} from {condition} animals",
                shared_processed, r1, r2,
            ])
            sra_writer.writerow([
                label, title, LIBRARY_STRATEGY, LIBRARY_SOURCE, LIBRARY_SELECTION,
                LIBRARY_LAYOUT, PLATFORM, args.instrument_model,
                "Stranded paired-end mRNA library, TruSeq-style dual-indexed",
                FILETYPE, r1, r2,
            ])

            for read_file in (r1, r2):
                if not read_file:
                    continue
                base = os.path.basename(read_file)
                checksum = known_md5s.get(base)
                if checksum is None:
                    full = os.path.join(args.fastq_dir, read_file)
                    if os.path.isfile(full):
                        checksum = md5_of(full)
                    else:
                        missing_md5.append(base)
                        continue
                md5_fh.write(f"{checksum}  {base}\n")

        for processed in processed_files:
            md5_fh.write(f"{md5_of(processed)}  {os.path.basename(processed)}\n")

    readme_path = os.path.join(args.output_dir, "SUBMISSION_README.txt")
    with open(readme_path, "w") as fh:
        fh.write("NCBI GEO / SRA submission package\n")
        fh.write("=" * 34 + "\n\n")
        fh.write(f"Species: {organism}\n")
        fh.write(f"Samples: {len(species_rows)}\n")
        if args.genome_build:
            fh.write(f"Genome build: {args.genome_build}\n")
        if args.annotation:
            fh.write(f"Annotation: {args.annotation}\n")
        fh.write("\nFiles in this directory\n")
        fh.write("  geo_samples.csv   - paste into the SAMPLES section of the GEO metadata workbook\n")
        fh.write("  sra_metadata.csv  - paste into the SRA metadata workbook\n")
        fh.write("  md5sums.txt       - checksums for every raw and processed file to upload\n")
        fh.write("\nRaw files to upload (from {}):\n".format(args.fastq_dir))
        fh.write("  one gzipped FASTQ per read per sample\n")
        fh.write("\nProcessed files to upload:\n")
        for processed in processed_files:
            fh.write(f"  {processed}\n")
        if missing_md5:
            fh.write("\nWARNING - no checksum available for:\n")
            for name in missing_md5:
                fh.write(f"  {name}\n")

    if missing_md5:
        print(f"Warning: {len(missing_md5)} raw files had no checksum and were not found on disk",
              file=sys.stderr)

    print(f"GEO sample sheet:  {geo_path}")
    print(f"SRA metadata:      {sra_path}")
    print(f"Checksums:         {md5_path}")
    print(f"Instructions:      {readme_path}")


if __name__ == "__main__":
    main()
