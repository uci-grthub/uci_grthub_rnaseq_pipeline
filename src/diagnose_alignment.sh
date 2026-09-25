#!/bin/bash
#SBATCH --job-name=diag_align
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=1:00:00
#SBATCH --error=logs/diag_align-%j.err
#SBATCH --output=logs/diag_align-%j.out

# Why: the production run finished 109/109 but every one of the 27 samples came
# back with a 1.4-1.9% hisat2 overall alignment rate. The reads themselves look
# healthy (97.7% survive trimming, full-length 151bp, high Q), so this isolates
# whether the cause is the pipeline's hisat2 flags, the grch38_snp_tran index,
# or the reads not being human at all. Vary one thing at a time, then confirm
# against a completely independent human reference (Salmon GRCh38 v33).

set -euo pipefail

cd /dfs9/ucightf-lab/projects/SotoJ/260909_SotoJ_mRNAseq-limmavoom

SAMPLE=xR106-L5-G3-P041-GATGCGTC-GCAGCCTC
IDX=/dfs8/commondata/hisat2/references/grch38_snp_tran/genome_snp_tran
SALMON_IDX=/dfs9/ucightf-lab/projects/REFERENCES/human/Salmon/GRCh38_v33_salmon_index
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# 50k pairs is plenty to separate 1% from 90%. `head` closes the pipe early,
# so zcat takes SIGPIPE -- suspend pipefail here or set -e kills the script.
set +o pipefail
zcat "output/trimmed/${SAMPLE}_trimmed_1P.fq.gz" | head -200000 > "$WORK/sub_1P.fq"
zcat "output/trimmed/${SAMPLE}_trimmed_2P.fq.gz" | head -200000 > "$WORK/sub_2P.fq"
zcat "data/FASTQ/${SAMPLE}-R1.fastq.gz"          | head -200000 > "$WORK/raw_1.fq"
zcat "data/FASTQ/${SAMPLE}-R2.fastq.gz"          | head -200000 > "$WORK/raw_2.fq"
set -o pipefail

wc -l "$WORK"/*.fq

module load hisat2/2.2.1

echo "############ 1. trimmed, pipeline flags ############"
hisat2 -p 8 --qc-filter --rna-strandness RF --dta-cufflinks \
    -x "$IDX" -1 "$WORK/sub_1P.fq" -2 "$WORK/sub_2P.fq" -S /dev/null

echo "############ 2. trimmed, default flags ############"
hisat2 -p 8 -x "$IDX" -1 "$WORK/sub_1P.fq" -2 "$WORK/sub_2P.fq" -S /dev/null

echo "############ 3. raw untrimmed, default flags ############"
hisat2 -p 8 -x "$IDX" -1 "$WORK/raw_1.fq" -2 "$WORK/raw_2.fq" -S /dev/null

echo "############ 4. R1 only, unpaired, default flags ############"
hisat2 -p 8 -x "$IDX" -U "$WORK/sub_1P.fq" -S /dev/null

module unload hisat2/2.2.1

echo "############ 5. Salmon vs independent human reference (GRCh38 v33) ############"
module load salmon/1.8.0
salmon quant -i "$SALMON_IDX" -l A -p 8 --validateMappings \
    -1 "$WORK/sub_1P.fq" -2 "$WORK/sub_2P.fq" -o "$WORK/salmon_out"
echo "--- salmon mapping rate ---"
grep -i "percent_mapped\|num_mapped\|num_processed" "$WORK/salmon_out/aux_info/meta_info.json"
module unload salmon/1.8.0
