#!/bin/bash
#SBATCH --job-name=diag_species
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=1:00:00
#SBATCH --error=logs/diag_species-%j.err
#SBATCH --output=logs/diag_species-%j.out

# Why: reads map 1.5% to the human genome (GRCh38) but 17.6% to the human
# transcriptome via Salmon. That gap is what a RELATED species looks like --
# conserved coding sequence cross-maps, intergenic and UTR sequence does not.
# Align the same 50k pairs to mouse GRCm38 head-to-head against human, and
# separately look at what the most abundant reads actually are, so a
# contamination explanation (rRNA, adapter dimer, polyA) can be ruled in or out.

set -euo pipefail

cd /dfs9/ucightf-lab/projects/SotoJ/260909_SotoJ_mRNAseq-limmavoom

SAMPLE=xR106-L5-G3-P041-GATGCGTC-GCAGCCTC
HUMAN=/dfs8/commondata/hisat2/references/grch38_snp_tran/genome_snp_tran
MOUSE=/dfs9/ucightf-lab/projects/REFERENCES/mouse/Hisat2/GRCm38/gencode.vM24.annotation.genome_tran
MM10=/dfs8/commondata/hisat2/references/mm10/genome
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

set +o pipefail
zcat "output/trimmed/${SAMPLE}_trimmed_1P.fq.gz" | head -200000 > "$WORK/sub_1P.fq"
zcat "output/trimmed/${SAMPLE}_trimmed_2P.fq.gz" | head -200000 > "$WORK/sub_2P.fq"
set -o pipefail

module load hisat2/2.2.1

echo "############ HUMAN GRCh38 (grch38_snp_tran) ############"
hisat2 -p 8 -x "$HUMAN" -1 "$WORK/sub_1P.fq" -2 "$WORK/sub_2P.fq" -S /dev/null

echo "############ MOUSE GRCm38 (gencode vM24 tran) ############"
hisat2 -p 8 -x "$MOUSE" -1 "$WORK/sub_1P.fq" -2 "$WORK/sub_2P.fq" -S /dev/null

echo "############ MOUSE mm10 (commondata) ############"
hisat2 -p 8 -x "$MM10" -1 "$WORK/sub_1P.fq" -2 "$WORK/sub_2P.fq" -S /dev/null

module unload hisat2/2.2.1

echo "############ most abundant R1 sequences (contamination check) ############"
# If this is adapter dimer, rRNA or polyA, a handful of sequences dominate.
awk 'NR%4==2' "$WORK/sub_1P.fq" | sort | uniq -c | sort -rn | head -15

echo "############ most abundant 30-mer prefixes ############"
awk 'NR%4==2{print substr($0,1,30)}' "$WORK/sub_1P.fq" | sort | uniq -c | sort -rn | head -15

echo "############ GC content distribution ############"
awk 'NR%4==2{n=length($0); g=gsub(/[GCgc]/,""); printf "%d\n", (g*100)/n}' "$WORK/sub_1P.fq" \
    | sort -n | uniq -c | awk '{printf "GC%%=%-4s %s\n", $2, $1}' | head -40
