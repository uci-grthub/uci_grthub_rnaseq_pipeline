#!/bin/bash
#SBATCH --job-name=rsync_fastq
#SBATCH --account=sbsandme_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/rsync_fastq/%x_%j.log

# Copy *.fastq.gz from SRC_DIR into DEST_DIR, then verify them against SRC_DIR's
# md5sums.txt (if present), which is kept as DEST_DIR/<SRC_DIR name>.md5sums.txt so
# several source dirs can share one DEST_DIR. Subdirectories and dot-files are skipped.
#
# Usage (from the project root):
#   sbatch src/rsync_fastq.sh SRC_DIR [DEST_DIR]
# DEST_DIR defaults to data/ and must be group ucightf with setgid set.
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: sbatch src/rsync_fastq.sh SRC_DIR [DEST_DIR]" >&2
    exit 1
fi

SRC=$(realpath "$1")
DEST=$(realpath "${2:-data}")

if ! compgen -G "${SRC}/*.fastq.gz" > /dev/null; then
    echo "ERROR: no *.fastq.gz in ${SRC}" >&2
    exit 1
fi

# -rlt, not -a: never copy perms/owner/group. Source dirs can be group
# ucightf_lab_share (1-byte quota on dfs9), and -a with "SRC/" also rewrites DEST's
# own mode and group, dropping its setgid bit. DEST is setgid ucightf, so new files
# inherit ucightf.
if [[ "$(stat -c '%G' "${DEST}")" != ucightf || ! -g "${DEST}" ]]; then
    echo "ERROR: ${DEST} must be group ucightf with setgid set" >&2
    exit 1
fi

echo "Copying ${SRC} -> ${DEST}"
rsync -rltv --partial --include='*.fastq.gz' --exclude='*' "${SRC}/" "${DEST}/"

if [[ -f "${SRC}/md5sums.txt" ]]; then
    MD5="${DEST}/$(basename "${SRC}").md5sums.txt"
    cp "${SRC}/md5sums.txt" "${MD5}"
    cd "${DEST}"
    md5sum -c "${MD5}"
else
    echo "WARNING: no md5sums.txt in ${SRC}; checksums not verified" >&2
fi
