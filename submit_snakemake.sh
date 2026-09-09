#!/bin/bash
# Launch the workflow with one SLURM job per rule.
#
# Snakemake submits and tracks the jobs itself, so this controller process must
# stay alive for the whole run. Run it under sbatch (as below) rather than in an
# interactive shell, where it would die with the terminal:
#
#     sbatch submit_snakemake.sh                      # default target: rule all
#     sbatch submit_snakemake.sh output/feature_count/human_samples_counts.txt
#
# The old --cluster/--cluster-config invocation this file used to carry was
# removed in Snakemake 8; job resources now come from profiles/slurm/config.yaml.
#SBATCH --job-name=snakemake_rnaseq
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --error=logs/snakemake-%j.err
#SBATCH --output=logs/snakemake-%j.out

set -euo pipefail

cd /dfs9/ucightf-lab/projects/SotoJ/260909_SotoJ_mRNAseq-limmavoom
mkdir -p logs

pixi run -e snakemake-dfs snakemake --profile profiles/slurm "$@"
