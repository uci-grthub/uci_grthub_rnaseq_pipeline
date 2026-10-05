#!/bin/bash
#SBATCH --job-name=build_ensdb
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --error=logs/build_ensdb-%j.err
#SBATCH --output=logs/build_ensdb-%j.out

##############################################################################
# Build the gene-annotation EnsDb for a reference GTF.
#
# Gene symbols and biotypes for the DE scripts come from an ensembldb EnsDb
# built from the same GTF the counts were made against. org.*.eg.db is keyed on
# Entrez and leaves ~13% of genes unlabelled (all lncRNA/TEC/pseudogene, plus
# some protein-coding); a TxDb cannot help because makeTxDbFromGFF() discards
# the gene_name attribute entirely.
#
# This is a script rather than a workflow rule because the database depends only
# on the GTF -- nothing about any one project -- so it is built once, stored
# beside the reference, and read by every project using that annotation. The
# paths below are the ones config.species_references.yaml points at.
#
# Usage:
#   bash src/build_gene_annotation.sh mouse         # build one species
#   bash src/build_gene_annotation.sh human
#   bash src/build_gene_annotation.sh all           # both
#   sbatch src/build_gene_annotation.sh mouse       # as a batch job
#
# The mouse GTF is 1.87M lines and rtracklayer::import() on it is the expensive
# step, so prefer sbatch over an interactive run. Measured on that GTF: 1m17s
# wall, 1.9 GB peak RSS (job 56794167), which is what the 8G request is sized
# against -- raise it for a substantially larger annotation.
#
# Re-run only when the reference GTF changes. Existing databases are left alone
# unless --force is given.
##############################################################################

set -euo pipefail

REFS=/dfs9/ucightf-lab/kstachel/REFERENCES

# species | gtf | output sqlite | ensembldb organism | genome build | release
# The human GTF lives on read-only commondata, so its database goes in the
# writable reference tree instead of beside the GTF.
read -r -d '' SPECIES_TABLE <<'EOF' || true
mouse|/dfs9/ucightf-lab/kstachel/REFERENCES/mouse/gencode.vM24.annotation.gtf|/dfs9/ucightf-lab/kstachel/REFERENCES/mouse/gencode.vM24.annotation.EnsDb.sqlite|Mus_musculus|GRCm38|100
human|/dfs8/commondata/hisat2/references/grch38_snp_tran/Homo_sapiens.GRCh38.84.gtf|/dfs9/ucightf-lab/kstachel/REFERENCES/human/Homo_sapiens.GRCh38.84.EnsDb.sqlite|Homo_sapiens|GRCh38|84
EOF

# sbatch copies the submitted script into a spool directory, so BASH_SOURCE
# points at /export/spool/... and not into the repo. SLURM_SUBMIT_DIR is where
# sbatch was invoked (the project root); fall back to BASH_SOURCE when running
# interactively.
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    SCRIPT_DIR="$SLURM_SUBMIT_DIR/src"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
if [ ! -r "$SCRIPT_DIR/gene_annotation.R" ]; then
    echo "Cannot find gene_annotation.R under $SCRIPT_DIR" >&2
    echo "Run this from the project root (bash src/build_gene_annotation.sh ...)" >&2
    exit 2
fi
FORCE=0
TARGET="${1:-all}"
[ "${2:-}" = "--force" ] && FORCE=1
[ "$TARGET" = "--force" ] && { FORCE=1; TARGET=all; }

build_one() {
    local species="$1" gtf="$2" out="$3" organism="$4" genome="$5" release="$6"

    echo "=== $species ==="
    if [ ! -r "$gtf" ]; then
        echo "  GTF not readable, skipping: $gtf" >&2
        return 1
    fi
    if [ -s "$out" ] && [ "$FORCE" -eq 0 ]; then
        echo "  already built, leaving alone: $out"
        echo "  (pass --force to rebuild)"
        return 0
    fi
    mkdir -p "$(dirname "$out")"

    echo "  GTF: $gtf"
    echo "  out: $out"
    module load R/4.5.2
    Rscript "$SCRIPT_DIR/gene_annotation.R" "$gtf" "$out" "$organism" "$genome" "$release"
    module unload R/4.5.2
}

status=0
while IFS='|' read -r species gtf out organism genome release; do
    [ -z "$species" ] && continue
    if [ "$TARGET" = "all" ] || [ "$TARGET" = "$species" ]; then
        build_one "$species" "$gtf" "$out" "$organism" "$genome" "$release" || status=1
    fi
done <<< "$SPECIES_TABLE"

if [ "$TARGET" != "all" ] && ! grep -q "^$TARGET|" <<< "$SPECIES_TABLE"; then
    echo "Unknown species '$TARGET'. Known: $(cut -d'|' -f1 <<< "$SPECIES_TABLE" | tr '\n' ' ')" >&2
    exit 2
fi

exit $status
