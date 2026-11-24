#!/bin/bash
#SBATCH --job-name=rnaseq_testdata
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=make_test_data-%j.out
#SBATCH --error=make_test_data-%j.err
#
# Rebuilds tests/data from public sources. Run from the repo root:
#   sbatch tests/make_test_data.sbatch
# Only needed when the test set itself changes; the result is committed.
#
# Reference: a 10 Mb window of GRCh38 chromosome 22 (Ensembl 84, the release
# of the production human GTF), renamed `chr22win` and shifted to start at 1.
# Reads: PRJNA929683 (HepG2 control vs AT7519, 4 x 2, 2x150 NovaSeq 6000),
# first 4M pairs per run, aligned to full GRCh38, pairs with both mates in the
# window kept and downsampled to READ_PAIRS per sample.
set -euo pipefail

WIN_START=30000000
WIN_END=40000000
READ_PAIRS=25000
FULL_INDEX=/dfs8/commondata/hisat2/references/grch38_snp_tran/genome_snp_tran
FULL_GTF=/dfs8/commondata/hisat2/references/grch38_snp_tran/Homo_sapiens.GRCh38.84.gtf
ENSEMBL_CHR22=https://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz

OUT=tests/data
# ~21 GB of uncompressed reads; keep it off the node-local /tmp.
WORK=${WORK:-$PWD/.testdata_work}
mkdir -p "$OUT/ref" "$OUT/fastq" "$WORK"
THREADS=${SLURM_CPUS_PER_TASK:-4}

module load samtools/1.15.1
module load hisat2/2.2.1

# --- reference window ------------------------------------------------------
curl -fsSL "$ENSEMBL_CHR22" | gunzip -c > "$WORK/chr22.fa"
samtools faidx "$WORK/chr22.fa"
samtools faidx "$WORK/chr22.fa" "22:$((WIN_START + 1))-$WIN_END" \
    | sed "1s/.*/>chr22win/" | gzip -n > "$OUT/ref/genome.fa.gz"

# Genes lying wholly inside the window, with every feature shifted.
awk -v s="$WIN_START" -v e="$WIN_END" -F'\t' 'BEGIN{OFS="\t"}
    $1=="22" && $3=="gene" && $4>s && $5<=e {
        match($9, /gene_id "[^"]+"/); keep[substr($9, RSTART, RLENGTH)] = 1 }
    $1=="22" { match($9, /gene_id "[^"]+"/); id = substr($9, RSTART, RLENGTH)
        if (id in keep) { $1 = "chr22win"; $4 -= s; $5 -= s; print } }' \
    <(awk -F'\t' '$1=="22"' "$FULL_GTF") <(awk -F'\t' '$1=="22"' "$FULL_GTF") \
    | awk '!seen[$0]++' | gzip -n > "$OUT/ref/annotation.gtf.gz"

# --- reads ----------------------------------------------------------------
declare -A RUNS=(
    [HepG2_ctrl_rep1]=SRR23280625 [HepG2_ctrl_rep2]=SRR23280624
    [HepG2_ctrl_rep3]=SRR23280623 [HepG2_ctrl_rep4]=SRR23280622
    [HepG2_AT7519_rep1]=SRR23280621 [HepG2_AT7519_rep2]=SRR23280620
    [HepG2_AT7519_rep3]=SRR23280619 [HepG2_AT7519_rep4]=SRR23280618
)
fetch() {  # fetch SAMPLE RUN MATE: first 4M reads of one FASTQ
    local base="https://ftp.sra.ebi.ac.uk/vol1/fastq/${2:0:6}/0${2: -2}/$2"
    # Closing the pipe early makes curl/gunzip fail; that is expected here.
    { curl -sSL "$base/${2}_$3.fastq.gz" | gunzip -c | head -n 16000000; } \
        > "$WORK/${1}_$3.fq" 2>/dev/null || true
}
for sample in "${!RUNS[@]}"; do
    for mate in 1 2; do fetch "$sample" "${RUNS[$sample]}" "$mate" & done
done
wait

for sample in "${!RUNS[@]}"; do
    [ "$(wc -l < "$WORK/${sample}_1.fq")" -eq 16000000 ] || { echo "short download: $sample" >&2; exit 1; }
    hisat2 -p "$THREADS" --no-unal -x "$FULL_INDEX" \
        -1 "$WORK/${sample}_1.fq" -2 "$WORK/${sample}_2.fq" 2> "$WORK/$sample.hisat2.log" \
        | samtools view -b -f 2 -q 10 - \
        | samtools sort -@ 4 -o "$WORK/$sample.bam" -
    samtools index "$WORK/$sample.bam"
    # Read names of pairs whose both mates fall in the window, then a fixed
    # random subset of them (shuf with a fixed random source is reproducible).
    samtools view "$WORK/$sample.bam" "22:$((WIN_START + 1))-$WIN_END" \
        | awk '$7=="=" && $8>'"$WIN_START"' && $8<='"$WIN_END"' {print $1}' | sort -u \
        | shuf -n "$READ_PAIRS" --random-source=<(yes 42) | sort > "$WORK/$sample.names"
    for mate in 1 2; do
        awk 'NR==FNR {keep["@"$1]=1; next}
             FNR%4==1 {split($1, a, /[ \/]/); p = (a[1] in keep)}
             p' "$WORK/$sample.names" "$WORK/${sample}_$mate.fq" \
            | gzip -n > "$OUT/fastq/${sample}-R$mate.fastq.gz"
    done
    rm -f "$WORK/${sample}"_?.fq "$WORK/$sample.bam"*
    echo "$sample $(grep overall "$WORK/$sample.hisat2.log")"
done

rm -rf "$WORK"
du -sh "$OUT"
