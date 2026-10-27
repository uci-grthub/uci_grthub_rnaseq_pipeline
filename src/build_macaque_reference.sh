#!/bin/bash
#SBATCH --job-name=build_macaque_ref
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --error=logs/build_macaque_ref-%j.err
#SBATCH --output=logs/build_macaque_ref-%j.out

##############################################################################
# Build the macaque references that config.species_references.yaml points at:
#
#   macaque      Mmul_10 (rheMac10), Ensembl release 116 annotation
#   macaque_t2t  T2T-MMU8v2.0 (GCF_049350105.2), NCBI RefSeq annotation
#                GCF_049350105.2-RS_2025_08
#
# Stages:
#
#   download  fetch genome FASTA + GTF and verify checksums. HPC3 compute
#             nodes reach ftp.ensembl.org and ftp.ncbi.nlm.nih.gov, so this
#             can run under sbatch as well as on a login node.
#   index     HISAT2 index with splice sites and exons (matching the existing
#             human/mouse "genome_tran" indices), a gffread transcript FASTA
#             cut from the same GTF + genome, and a decoy-aware Salmon index.
#   all       download, then index.
#
# Usage:
#   sbatch src/build_macaque_reference.sh all macaque
#   sbatch src/build_macaque_reference.sh all macaque_t2t
#   bash   src/build_macaque_reference.sh download macaque   # just the files
#
# hisat2-build --ss --exon on a ~3 Gb primate genome peaks around 160-200 GB,
# which is what the 240G request is sized against; it lands on the 256G/515G
# standard nodes. Each step is skipped when its output already exists unless
# --force is given, so a failed index job can be resubmitted as is.
#
# After indexing, build the gene annotation database with
#   sbatch src/build_gene_annotation.sh macaque
##############################################################################

set -euo pipefail

REFS=/dfs9/ucightf-lab/kstachel/REFERENCES/macaque

ENSEMBL_RELEASE=116
ENSEMBL_BASE=https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}
NCBI_ACCESSION=GCF_049350105.2_T2T-MMU8v2.0
NCBI_BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/049/350/105/${NCBI_ACCESSION}

# See build_gene_annotation.sh: under sbatch BASH_SOURCE points into the spool.
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
    PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

STAGE="${1:-}"
SPECIES="${2:-}"
FORCE=0
[ "${3:-}" = "--force" ] && FORCE=1

# Per-species layout. Unpacked file names keep the upstream names so the
# provenance of a path in the config is obvious.
case "$SPECIES" in
    macaque)
        OUT="$REFS/Mmul_10"
        GENOME_URL="$ENSEMBL_BASE/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna_sm.toplevel.fa.gz"
        GTF_URL="$ENSEMBL_BASE/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.${ENSEMBL_RELEASE}.gtf.gz"
        ;;
    macaque_t2t)
        OUT="$REFS/T2T-MMU8v2.0"
        GENOME_URL="$NCBI_BASE/${NCBI_ACCESSION}_genomic.fna.gz"
        GTF_URL="$NCBI_BASE/${NCBI_ACCESSION}_genomic.gtf.gz"
        ;;
    *)
        echo "Usage: $0 {download|index|all} {macaque|macaque_t2t} [--force]" >&2
        exit 2
        ;;
esac

GENOME_GZ="$OUT/$(basename "$GENOME_URL")"
GTF_GZ="$OUT/$(basename "$GTF_URL")"
GENOME="${GENOME_GZ%.gz}"
GTF="${GTF_GZ%.gz}"
HISAT2_PREFIX="$OUT/hisat2/genome_tran"
TRANSCRIPTS="$OUT/transcripts.gffread.fa"
SALMON_INDEX="$OUT/salmon_index"

needs_build() {
    [ "$FORCE" -eq 1 ] || [ ! -s "$1" ]
}

# Ensembl publishes BSD `sum` checksums (CHECKSUMS); NCBI publishes md5
# (md5checksums.txt). Either way a mismatch aborts before anything is indexed.
verify_checksum() {
    local file="$1" name
    name="$(basename "$file")"
    if [ "$SPECIES" = "macaque" ]; then
        local sums_url="$(dirname "$2")/CHECKSUMS" expected actual
        expected="$(curl -fsS "$sums_url" | awk -v f="$name" '$3 == f {print $1, $2}')"
        actual="$(sum "$file" | awk '{print $1, $2}')"
        if [ -z "$expected" ] || [ "$expected" != "$actual" ]; then
            echo "  checksum mismatch for $name: expected '$expected', got '$actual'" >&2
            return 1
        fi
    else
        local expected actual
        expected="$(curl -fsS "$NCBI_BASE/md5checksums.txt" | awk -v f="./$name" '$2 == f {print $1}')"
        actual="$(md5sum "$file" | awk '{print $1}')"
        if [ -z "$expected" ] || [ "$expected" != "$actual" ]; then
            echo "  md5 mismatch for $name: expected '$expected', got '$actual'" >&2
            return 1
        fi
    fi
    echo "  checksum ok: $name"
}

fetch() {
    local url="$1" gz="$2"
    if needs_build "${gz%.gz}"; then
        echo "  fetching $url"
        curl -fsSL --retry 3 -o "$gz" "$url"
        verify_checksum "$gz" "$url"
        gunzip -f "$gz"
    else
        echo "  already present, leaving alone: ${gz%.gz}"
    fi
}

do_download() {
    mkdir -p "$OUT"
    echo "=== download $SPECIES -> $OUT ==="
    fetch "$GENOME_URL" "$GENOME_GZ"
    fetch "$GTF_URL" "$GTF_GZ"
    {
        echo "species: $SPECIES"
        echo "genome: $GENOME_URL"
        echo "gtf: $GTF_URL"
        echo "downloaded: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    } > "$OUT/SOURCE.txt"
    (cd "$OUT" && md5sum "$(basename "$GENOME")" "$(basename "$GTF")" > md5sums.txt)
    echo "  wrote $OUT/SOURCE.txt and $OUT/md5sums.txt"
}

do_index() {
    for f in "$GENOME" "$GTF"; do
        if [ ! -s "$f" ]; then
            echo "Missing $f; run the download stage first." >&2
            exit 1
        fi
    done
    local threads="${SLURM_CPUS_PER_TASK:-8}"
    echo "=== index $SPECIES ($threads threads) ==="

    mkdir -p "$(dirname "$HISAT2_PREFIX")"
    if needs_build "$HISAT2_PREFIX.1.ht2"; then
        module load hisat2/2.2.1
        hisat2_extract_splice_sites.py "$GTF" > "$OUT/hisat2/splice_sites.txt"
        hisat2_extract_exons.py "$GTF" > "$OUT/hisat2/exons.txt"
        hisat2-build -p "$threads" \
            --ss "$OUT/hisat2/splice_sites.txt" \
            --exon "$OUT/hisat2/exons.txt" \
            "$GENOME" "$HISAT2_PREFIX"
        module unload hisat2/2.2.1
    else
        echo "  HISAT2 index already built: $HISAT2_PREFIX"
    fi

    # Transcripts are cut from the same GTF + genome so Salmon's transcript IDs
    # are exactly the GTF's transcript_id values that tximport maps to genes.
    if needs_build "$TRANSCRIPTS"; then
        pixi run --manifest-path "$PROJECT_ROOT/pixi.toml" -e refbuild \
            gffread -w "$TRANSCRIPTS" -g "$GENOME" "$GTF"
    else
        echo "  transcript FASTA already built: $TRANSCRIPTS"
    fi

    if needs_build "$SALMON_INDEX/versionInfo.json"; then
        module load salmon/1.8.0
        bash "$PROJECT_ROOT/src/build_salmon_index.sh" \
            -t "$TRANSCRIPTS" -g "$GENOME" -o "$SALMON_INDEX" -p "$threads"
        module unload salmon/1.8.0
    else
        echo "  Salmon index already built: $SALMON_INDEX"
    fi
}

case "$STAGE" in
    download) do_download ;;
    index) do_index ;;
    all) do_download && do_index ;;
    *)
        echo "Usage: $0 {download|index|all} {macaque|macaque_t2t} [--force]" >&2
        exit 2
        ;;
esac
