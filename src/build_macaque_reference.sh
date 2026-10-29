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
#   macaque_refseq
#                Mmul_10 (GCF_003339765.1), NCBI RefSeq Annotation Release
#                103. No HISAT2 index: it is the same assembly as the Ensembl
#                Mmul_10 genome, so the macaque BAMs are counted against this
#                GTF with its sequence names translated to Ensembl's
#                (*.ensembl_names.gtf), and Salmon quantifies RefSeq isoforms,
#                e.g. the three HTR2A transcripts Ensembl lacks.
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
#   sbatch --mem=64G src/build_macaque_reference.sh all macaque_refseq
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
BUILD_HISAT2=1

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
        NCBI_ACCESSION=GCF_049350105.2_T2T-MMU8v2.0
        NCBI_BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/049/350/105/${NCBI_ACCESSION}
        GENOME_URL="$NCBI_BASE/${NCBI_ACCESSION}_genomic.fna.gz"
        GTF_URL="$NCBI_BASE/${NCBI_ACCESSION}_genomic.gtf.gz"
        ;;
    macaque_refseq)
        OUT="$REFS/Mmul_10_RefSeq"
        NCBI_ACCESSION=GCF_003339765.1_Mmul_10
        NCBI_BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/${NCBI_ACCESSION}
        GENOME_URL="$NCBI_BASE/${NCBI_ACCESSION}_genomic.fna.gz"
        GTF_URL="$NCBI_BASE/${NCBI_ACCESSION}_genomic.gtf.gz"
        BUILD_HISAT2=0
        ;;
    *)
        echo "Usage: $0 {download|index|all} {macaque|macaque_t2t|macaque_refseq} [--force]" >&2
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
ENSEMBL_NAMES_GTF="${GTF%.gtf}.ensembl_names.gtf"
ENSEMBL_FAI="$REFS/Mmul_10/Macaca_mulatta.Mmul_10.dna_sm.toplevel.fa.fai"

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

    if [ "$BUILD_HISAT2" -eq 0 ]; then
        echo "  HISAT2 index not built for $SPECIES"
    elif needs_build "$HISAT2_PREFIX.1.ht2"; then
        mkdir -p "$(dirname "$HISAT2_PREFIX")"
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
        # NCBI GTFs give `gene` lines an empty transcript_id, which gffread
        # rejects ("no valid ID found for GFF record"); transcripts are cut
        # from the exon lines, so gene lines are dropped first.
        local gffread_gtf="$GTF"
        if [ "$SPECIES" != "macaque" ]; then
            gffread_gtf="$OUT/gffread_input.gtf"
            awk -F'\t' '$3 != "gene"' "$GTF" > "$gffread_gtf"
        fi
        pixi run --manifest-path "$PROJECT_ROOT/pixi.toml" -e refbuild \
            gffread -w "$TRANSCRIPTS" -g "$GENOME" "$gffread_gtf"
        [ "$gffread_gtf" = "$GTF" ] || rm -f "$gffread_gtf"
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

    if [ "$SPECIES" = "macaque_refseq" ]; then
        if needs_build "$ENSEMBL_NAMES_GTF"; then
            rename_gtf_to_ensembl
        else
            echo "  Ensembl-named GTF already built: $ENSEMBL_NAMES_GTF"
        fi
    fi
}

# RefSeq and Ensembl Mmul_10 are one assembly under two naming schemes
# (NC_041754.1 vs 1, NW_021160096.1 vs QNVO02000261.1). The NCBI assembly
# report maps each RefSeq accession to its chromosome (Assigned-Molecule) or
# scaffold GenBank accession, which are Ensembl's names. Every mapped name must
# exist in the Ensembl genome with the same length, or the BAM coordinates
# would not line up with the annotation.
rename_gtf_to_ensembl() {
    local report="$OUT/${NCBI_ACCESSION}_assembly_report.txt"
    local map="$OUT/refseq_to_ensembl_names.tsv"
    curl -fsSL --retry 3 -o "$report" "$NCBI_BASE/${NCBI_ACCESSION}_assembly_report.txt"
    awk -F'\t' '!/^#/ {
        name = ($2 == "assembled-molecule") ? $3 : $5
        print $7 "\t" name "\t" $9
    }' "$report" > "$map"

    local bad
    bad="$(awk -F'\t' 'NR == FNR { len[$1] = $2; next }
        !($2 in len) || len[$2] != $3 { print $1 " -> " $2 }' "$ENSEMBL_FAI" "$map")"
    if [ -n "$bad" ]; then
        echo "  RefSeq names with no equal-length Ensembl sequence:" >&2
        echo "$bad" | head >&2
        return 1
    fi

    awk -F'\t' -v OFS='\t' 'NR == FNR { name[$1] = $2; next }
        /^#/ { print; next }
        !($1 in name) { print "unmapped seqname " $1 > "/dev/stderr"; exit 1 }
        { $1 = name[$1]; print }' "$map" "$GTF" > "$ENSEMBL_NAMES_GTF.tmp"
    mv "$ENSEMBL_NAMES_GTF.tmp" "$ENSEMBL_NAMES_GTF"
    echo "  wrote $ENSEMBL_NAMES_GTF ($(wc -l < "$map") sequences mapped)"
}

case "$STAGE" in
    download) do_download ;;
    index) do_index ;;
    all) do_download && do_index ;;
    *)
        echo "Usage: $0 {download|index|all} {macaque|macaque_t2t|macaque_refseq} [--force]" >&2
        exit 2
        ;;
esac
