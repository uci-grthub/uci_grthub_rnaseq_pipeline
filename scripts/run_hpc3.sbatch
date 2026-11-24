#!/bin/bash
#SBATCH --job-name=rnaseq_controller
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=3-00:00:00
#SBATCH --output=snakemake-%j.out
#SBATCH --error=snakemake-%j.err
#
# Snakemake controller for one project. Submit it from the project directory
# (the one holding config.yaml and the sample sheet):
#   cd /path/to/project
#   sbatch -A <account> /path/to/uci_grthub_rnaseq_pipeline-vX.Y.Z/scripts/run_hpc3.sbatch
# This job only schedules; each rule runs as its own Slurm job through
# profiles/hpc3. Extra arguments go to snakemake, e.g.
#   sbatch .../run_hpc3.sbatch --software-deployment-method env-modules
#   sbatch .../run_hpc3.sbatch -n
set -euo pipefail

PROJECT_DIR=${SLURM_SUBMIT_DIR:-$PWD}
# sbatch runs a spooled copy of this script, so find the release through the
# PIPELINE_DIR variable or the path it was submitted with.
PIPELINE_DIR=${PIPELINE_DIR:-$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' | xargs dirname | xargs dirname)}
export SLURM_ACCOUNT=${SLURM_ACCOUNT:-$SLURM_JOB_ACCOUNT}

[ -f "$PROJECT_DIR/config.yaml" ] || { echo "No config.yaml in $PROJECT_DIR" >&2; exit 2; }

pixi run --frozen --manifest-path "$PIPELINE_DIR/pixi.toml" snakemake \
    -s "$PIPELINE_DIR/workflow/Snakefile" \
    --directory "$PROJECT_DIR" \
    --configfile "$PROJECT_DIR/config.yaml" \
    --profile "$PIPELINE_DIR/profiles/hpc3" \
    "$@"
