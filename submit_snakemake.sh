#!/bin/bash
#SBATCH --job-name=snakemake_rnaseq
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --error=snakemake-%j.err
#SBATCH --output=snakemake-%j.out

# This job only runs the Snakemake controller; each rule is submitted as its
# own Slurm job by snakemake-executor-plugin-slurm. Snakemake, the plugin and
# chelae come from the pixi `snakemake-dfs` environment.
#
# The rules' `partition`/`account` resources are not read by the plugin, so
# the Slurm account and partition are set here as defaults. `runtime` needs a
# unit (a bare number is read as seconds) and applies to every rule.
pixi run -e snakemake-dfs snakemake \
    --executor slurm \
    --jobs 10 \
    --default-resources slurm_account=sbsandme_lab slurm_partition=standard runtime=12h \
    --latency-wait 60 \
    --rerun-incomplete \
    --printshellcmds
