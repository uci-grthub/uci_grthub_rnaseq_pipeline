#!/bin/bash
# Run the workflow on the test set and compare against tests/expected.
#   tests/run_test.sh [--cores N] [extra snakemake args...]
# Run from the repo root inside pixi (`pixi run tests/run_test.sh`), or let
# CI do it. Outputs go to tests/project/ (gitignored).
set -euo pipefail
cd "$(dirname "$0")/.."
CORES=2
if [ "${1:-}" = "--cores" ]; then CORES=$2; shift 2; fi

[ -f tests/data/ref/index/hisat2/genome_tran.1.ht2 ] || \
    pixi run -e tools tests/build_mini_indexes.sh

CONFIGS=(tests/config/config.yaml)
# HPC3's glibc is too old for bioconda rustqc; use the cluster image there.
[ -f /dfs9/ucightf-lab/kstachel/containers/rustqc.sif ] && CONFIGS+=(tests/config/hpc3.yaml)

mkdir -p tests/project
snakemake -s workflow/Snakefile --directory tests/project \
    --configfile "${CONFIGS[@]}" \
    --cores "$CORES" --printshellcmds --show-failed-logs \
    "$@" \
    all output/rmats/test_mini/.done output/ncbi_submission/test_mini/geo_samples.csv

python3 tests/compare_outputs.py tests/project tests/expected
