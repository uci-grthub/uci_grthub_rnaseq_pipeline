#!/bin/bash
#SBATCH --job-name=rnaseq_regression
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=3-00:00:00
#SBATCH --output=regression-%j.out
#
# Run a pipeline version on the full-size regression set and compare it with
# the stored reference outputs. Submit from the repo root:
#   sbatch tests/regression/run_regression.sbatch REF REG_DIR [--write-expected]
# REF is a tag or commit. With LEGACY=1, REF is run in the pre-1.0 layout
# (Snakefile at the repo root, config.local.yaml).
# This job is only the Snakemake controller; rules run as their own Slurm jobs.
set -euo pipefail
REF=${1:?usage: run_regression.sbatch REF REG_DIR [--write-expected]}
REG=${2:?usage: run_regression.sbatch REF REG_DIR [--write-expected]}
MODE=${3:-}
REPO=$PWD
RUN=$REG/runs/$REF
export SLURM_ACCOUNT=${SLURM_ACCOUNT:-${SLURM_JOB_ACCOUNT:-sbsandme_lab}}

[ -d "$RUN/pipeline" ] || git -C "$REPO" worktree add --detach "$RUN/pipeline" "$REF"
pixi install --manifest-path "$RUN/pipeline/pixi.toml" --all

if [ "${LEGACY:-0}" = 1 ]; then
    PROJECT=$RUN/pipeline
    mkdir -p "$PROJECT/data" "$PROJECT/metadata"
    ln -sfn "$REG/data/FASTQ" "$PROJECT/data/FASTQ"
    cp "$REG/samples.csv" "$PROJECT/metadata/metadata.csv"
    printf 'default_species: human\nncbi:\n  tissue: "liver"\n' > "$PROJECT/config.local.yaml"
    (cd "$PROJECT" && pixi run -e snakemake-dfs snakemake --executor slurm --jobs 50 \
        --default-resources slurm_account="$SLURM_ACCOUNT" slurm_partition=standard runtime=12h \
        --latency-wait 60 --rerun-incomplete \
        all output/rmats/human/.done output/ncbi_submission/human/geo_samples.csv)
else
    PROJECT=$RUN/project
    mkdir -p "$PROJECT"
    sed -e "s|../../data/FASTQ|$REG/data/FASTQ|; s|../../samples.csv|$REG/samples.csv|" \
        "$REPO/tests/regression/config.yaml" > "$PROJECT/config.yaml"
    pixi run --manifest-path "$RUN/pipeline/pixi.toml" snakemake \
        -s "$RUN/pipeline/workflow/Snakefile" --directory "$PROJECT" \
        --configfile "$PROJECT/config.yaml" --profile "$RUN/pipeline/profiles/hpc3" \
        all output/rmats/human/.done output/ncbi_submission/human/geo_samples.csv
fi

if [ "$MODE" = --write-expected ]; then
    pixi run --manifest-path "$REPO/pixi.toml" python3 "$REPO/tests/compare_outputs.py" \
        --write-expected "$PROJECT" "$REG/expected"
else
    pixi run --manifest-path "$REPO/pixi.toml" python3 "$REPO/tests/compare_outputs.py" \
        "$PROJECT" "$REG/expected" | tee "$RUN/comparison.txt"
fi
