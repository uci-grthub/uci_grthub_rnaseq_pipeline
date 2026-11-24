#!/bin/bash
#SBATCH --job-name=rnaseq_v090_baseline
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=make_expected-%j.out
#SBATCH --error=make_expected-%j.err
#
# Produce tests/expected from the pre-modernization pipeline (tag v0.9.0) on
# the test set, with the tools v0.9.0 used (HPC3 modules + .sif images).
# Run from the repo root: sbatch tests/make_expected_v0.9.0.sbatch
# Note: v0.9.0 mails its hard-coded address on start/finish.
set -euo pipefail
REPO=$PWD
BASE=${BASE:-$REPO/.baseline_v0.9.0}

[ -f tests/data/ref/index/hisat2/genome_tran.1.ht2 ] || \
    pixi run -e tools env THREADS=8 tests/build_mini_indexes.sh

[ -d "$BASE" ] || git worktree add --detach "$BASE" v0.9.0
mkdir -p "$BASE/data" "$BASE/metadata"
ln -sfn "$REPO/tests/data/fastq" "$BASE/data/FASTQ"
cp tests/config/samples.csv "$BASE/metadata/metadata.csv"
IDX=$REPO/tests/data/ref/index
cat > "$BASE/config.local.yaml" <<YAML
default_species: test_mini
params:
  cpus: 8
ncbi:
  tissue: "liver"
  instrument_model: "Illumina NovaSeq 6000"
species_references:
  test_mini:
    organism: "Homo sapiens"
    genome_build: "GRCh38 chr22:30-40Mb window (test)"
    annotation: "Ensembl release 84"
    hisat2_index: "$IDX/hisat2/genome_tran"
    gtf: "$IDX/annotation.gtf"
    salmon_index: "$IDX/salmon/salmon_index"
    ensdb: ""
YAML

cd "$BASE"
pixi run -e snakemake-dfs snakemake --cores 8 --printshellcmds --show-failed-logs \
    all output/rmats/test_mini/.done output/ncbi_submission/test_mini/geo_samples.csv
cd "$REPO"
pixi run python3 tests/compare_outputs.py --write-expected "$BASE" tests/expected
