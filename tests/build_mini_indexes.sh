#!/bin/bash
# Build the HISAT2 and Salmon indexes for the test reference.
#   tests/build_mini_indexes.sh [OUT_DIR]   (default: tests/data/ref/index)
# Needs hisat2-build, salmon and gffread on PATH (the pixi `pipeline` env, or
# the HPC3 modules hisat2/2.2.1 + salmon/1.8.0 plus pixi `refbuild`).
# Indexes are not committed; CI caches them keyed on the FASTA/GTF hashes.
set -euo pipefail
REF=$(cd "$(dirname "$0")/data/ref" && pwd)
OUT=${1:-$REF/index}
THREADS=${THREADS:-2}
mkdir -p "$OUT/hisat2" "$OUT/salmon"

gunzip -c "$REF/genome.fa.gz" > "$OUT/genome.fa"
gunzip -c "$REF/annotation.gtf.gz" > "$OUT/annotation.gtf"

hisat2_extract_splice_sites.py "$OUT/annotation.gtf" > "$OUT/splicesites.txt"
hisat2_extract_exons.py "$OUT/annotation.gtf" > "$OUT/exons.txt"
hisat2-build -q -p "$THREADS" --seed 42 \
    --ss "$OUT/splicesites.txt" --exon "$OUT/exons.txt" \
    "$OUT/genome.fa" "$OUT/hisat2/genome_tran"

gffread -w "$OUT/transcripts.fa" -g "$OUT/genome.fa" "$OUT/annotation.gtf"
salmon index -t "$OUT/transcripts.fa" -i "$OUT/salmon/salmon_index" -p "$THREADS" -k 31 > "$OUT/salmon_index.log" 2>&1
echo "indexes in $OUT"
