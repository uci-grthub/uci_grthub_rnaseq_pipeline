#!/bin/bash
#SBATCH --job-name=tximport_tx
#SBATCH --account=sbsandme_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/tximport_tpm/mouse_transcript_%j.log

set -euo pipefail
cd /dfs9/ucightf-lab/projects/ChatS/260818_ChatS_RNAseq

module load R/4.5.2
Rscript src/tximport_tpm.R \
    output/salmon_by_species/mouse \
    "$(python3 -c "import yaml;print(yaml.safe_load(open('config.species_references.yaml'))['species_references']['mouse']['gtf'])")" \
    output/tpm/mouse
