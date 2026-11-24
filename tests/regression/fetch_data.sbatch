#!/bin/bash
#SBATCH --job-name=rnaseq_regression_fetch
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=regression_fetch-%j.out
# Download the full PRJNA929683 FASTQs into REG/data/FASTQ and write REG/samples.csv.
#   sbatch tests/regression/fetch_data.sbatch REG
set -euo pipefail
REG=${1:?usage: fetch_data.sbatch REG_DIR}
mkdir -p "$REG/data/FASTQ"
declare -A RUNS=(
    [HepG2_ctrl_rep1]=SRR23280625 [HepG2_ctrl_rep2]=SRR23280624
    [HepG2_ctrl_rep3]=SRR23280623 [HepG2_ctrl_rep4]=SRR23280622
    [HepG2_AT7519_rep1]=SRR23280621 [HepG2_AT7519_rep2]=SRR23280620
    [HepG2_AT7519_rep3]=SRR23280619 [HepG2_AT7519_rep4]=SRR23280618
)
echo "sample,i7barcode,i5barcode_NovaSeqV1.5,sex,condition,species" > "$REG/samples.csv"
for sample in $(printf '%s\n' "${!RUNS[@]}" | sort); do
    run=${RUNS[$sample]}
    base="https://ftp.sra.ebi.ac.uk/vol1/fastq/${run:0:6}/0${run: -2}/$run"
    for mate in 1 2; do
        out="$REG/data/FASTQ/${sample}-R$mate.fastq.gz"
        [ -s "$out" ] || curl -fsSL --retry 5 -o "$out" "$base/${run}_$mate.fastq.gz"
    done
    condition=${sample#HepG2_}; condition=${condition%_rep*}
    echo "$sample,,,,$condition,human" >> "$REG/samples.csv"
done
(cd "$REG/data/FASTQ" && md5sum ./*.fastq.gz > md5sums.txt)
