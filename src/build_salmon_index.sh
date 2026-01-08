#!/bin/bash

##############################################################################
# Build a decoy-aware salmon index for Drosophila melanogaster
# 
# This script creates a salmon index using a transcriptome and genome
# as decoys.
#
# Usage: bash build_salmon_index.sh [OPTIONS]
#
# Options:
#   -t, --transcriptome PATH   Path to transcriptome FASTA file (gzipped or not)
#                              Default: ../data/dmel-all-transcript-r6.53.fasta.gz
#   -g, --genome PATH          Path to genome FASTA file (for decoys)
#                              Default: ../data/dmel-all-chromosome-r6.53.fasta
#   -o, --output DIR           Output directory for salmon index
#                              Default: ../salmon_index
#   -h, --help                 Show this help message
##############################################################################

set -euo pipefail

# Display usage
usage() {
    grep "^# " "$0" | grep -E "^\#\s+(Options|Usage)" -A 100 | sed 's/^# //' | sed 's/^##*//'
    exit 0
}

# Define default paths relative to script location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "${SCRIPT_DIR}")"

# Default values
TRANSCRIPTOME_GZ_DEFAULT="${PROJECT_ROOT}/data/dmel-all-transcript-r6.53.fasta.gz"
GENOME_DEFAULT="${PROJECT_ROOT}/data/dmel-all-chromosome-r6.53.fasta"
INDEX_DIR_DEFAULT="${PROJECT_ROOT}/salmon_index"

# Parse command-line arguments
TRANSCRIPTOME_GZ="${TRANSCRIPTOME_GZ_DEFAULT}"
GENOME="${GENOME_DEFAULT}"
INDEX_DIR="${INDEX_DIR_DEFAULT}"

while [[ $# -gt 0 ]]; do
    case $1 in
        -t|--transcriptome)
            TRANSCRIPTOME_GZ="$2"
            shift 2
            ;;
        -g|--genome)
            GENOME="$2"
            shift 2
            ;;
        -o|--output)
            INDEX_DIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            usage
            ;;
    esac
done

# Set up other paths
TEMP_DIR="${INDEX_DIR%/*}/salmon_temp"
TRANSCRIPTOME="${TEMP_DIR}/$(basename "${TRANSCRIPTOME_GZ}" .gz)"
DECOYS_FILE="${TEMP_DIR}/decoys.txt"
GENTROME="${TEMP_DIR}/gentrome.fasta"

echo "================================"
echo "Building Decoy-Aware Salmon Index"
echo "================================"
echo ""
echo "Configuration:"
echo "  Transcriptome: ${TRANSCRIPTOME_GZ}"
echo "  Genome:        ${GENOME}"
echo "  Index output:  ${INDEX_DIR}"
echo ""

# Check if required files exist
echo "[1/6] Checking for required files..."
if [[ ! -f "${TRANSCRIPTOME_GZ}" ]]; then
    echo "ERROR: Transcriptome file not found: ${TRANSCRIPTOME_GZ}"
    exit 1
fi

if [[ ! -f "${GENOME}" ]]; then
    echo "ERROR: Genome file not found: ${GENOME}"
    exit 1
fi
echo "✓ Required files found"
echo ""

# Create directories
echo "[2/6] Creating working directories..."
mkdir -p "${TEMP_DIR}" "${INDEX_DIR}"
echo "✓ Directories created: ${TEMP_DIR}, ${INDEX_DIR}"
echo ""

# Decompress transcriptome
echo "[3/6] Decompressing transcriptome..."
if [[ ! -f "${TRANSCRIPTOME}" ]]; then
    gunzip -c "${TRANSCRIPTOME_GZ}" > "${TRANSCRIPTOME}"
    echo "✓ Transcriptome decompressed to ${TRANSCRIPTOME}"
else
    echo "✓ Transcriptome already decompressed"
fi
echo ""

# Generate decoys file (names of all sequences in genome)
echo "[4/6] Generating decoys file..."
grep "^>" "${GENOME}" | cut -d " " -f 1 | sed 's/^>//' > "${DECOYS_FILE}"
NUM_DECOYS=$(wc -l < "${DECOYS_FILE}")
echo "✓ Decoys file created with ${NUM_DECOYS} entries"
echo ""

# Combine transcriptome and genome into gentrome
echo "[5/6] Creating gentrome (transcriptome + genome)..."
cat "${TRANSCRIPTOME}" "${GENOME}" > "${GENTROME}"
echo "✓ Gentrome created: ${GENTROME}"
echo ""

# Build salmon index
echo "[6/6] Building salmon index..."
echo "   Command: salmon index -t ${GENTROME} -d ${DECOYS_FILE} -i ${INDEX_DIR} -k 31"
salmon index \
    --transcripts "${GENTROME}" \
    --decoys "${DECOYS_FILE}" \
    --index "${INDEX_DIR}" \
    --threads 8

echo ""
echo "================================"
echo "✓ Salmon index build complete!"
echo "================================"
echo ""
echo "Index location: ${INDEX_DIR}"
echo "To use this index for quantification:"
echo "  salmon quant -i ${INDEX_DIR} -l A -1 reads_R1.fq -2 reads_R2.fq -o output_quant"
echo ""
