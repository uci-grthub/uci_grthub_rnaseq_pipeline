#!/bin/bash
#SBATCH --job-name=screen_species
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --error=logs/screen_species-%j.err
#SBATCH --output=logs/screen_species-%j.out

# Why: sample P041 aligns 1.7% to human GRCh38 but 98.9% to mouse GRCm38, so the
# run was aligned against the wrong genome. Before refilling the metadata
# species column, confirm every sample individually -- the whole reason
# config.species_references.yaml exists is that a previous run (xR098 P070/P071)
# turned out to be a MIXED human/mouse batch. Do not assume all 27 match P041.

set -euo pipefail

cd /dfs9/ucightf-lab/projects/SotoJ/260909_SotoJ_mRNAseq-limmavoom

HUMAN=/dfs8/commondata/hisat2/references/grch38_snp_tran/genome_snp_tran
MOUSE=/dfs9/ucightf-lab/projects/REFERENCES/mouse/Hisat2/GRCm38/gencode.vM24.annotation.genome_tran
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

module load hisat2/2.2.1

OUT=metadata/species_screen.tsv
printf "sample\thuman_pct\tmouse_pct\tcall\n" > "$OUT"

rate() {  # parse hisat2's stderr summary into a bare percentage
    grep -o '[0-9.]*% overall alignment rate' | grep -o '[0-9.]*'
}

for r1 in output/trimmed/*_trimmed_1P.fq.gz; do
    s=$(basename "$r1" _trimmed_1P.fq.gz)

    # 20k pairs is ample to separate ~2% from ~99%.
    set +o pipefail
    zcat "$r1" | head -80000 > "$WORK/r1.fq"
    zcat "output/trimmed/${s}_trimmed_2P.fq.gz" | head -80000 > "$WORK/r2.fq"
    set -o pipefail

    h=$(hisat2 -p 8 -x "$HUMAN" -1 "$WORK/r1.fq" -2 "$WORK/r2.fq" -S /dev/null 2>&1 | rate)
    m=$(hisat2 -p 8 -x "$MOUSE" -1 "$WORK/r1.fq" -2 "$WORK/r2.fq" -S /dev/null 2>&1 | rate)

    call=$(awk -v h="$h" -v m="$m" 'BEGIN{
        if (m > h + 20)      print "mouse";
        else if (h > m + 20) print "human";
        else                 print "AMBIGUOUS";
    }')
    printf "%s\t%s\t%s\t%s\n" "$s" "$h" "$m" "$call" | tee -a "$OUT"
done

module unload hisat2/2.2.1

echo "############ SUMMARY ############"
awk -F'\t' 'NR>1{c[$4]++} END{for (k in c) printf "%-10s %d\n", k, c[k]}' "$OUT"
