# Full-size regression (HPC3)

Run before every release PR is merged. The dataset is the full public
PRJNA929683 study (HepG2 control vs AT7519, 4 replicates each, 2×150 bp,
about 30–46M pairs per run) on the production human reference (GRCh38,
Ensembl 84), so no consent is involved. The test set is a subsample of the
same study.

One-time setup: download the FASTQs and record reference outputs from the
last release:

```bash
REG=/dfs9/ucightf-lab/$USER/regression/rnaseq_PRJNA929683
sbatch tests/regression/fetch_data.sbatch "$REG"           # ~60 GB
# Reference outputs from the previous release (v0.9.0 for the 1.0.0 cutover;
# legacy in-repo layout, so run it with LEGACY=1):
LEGACY=1 sbatch tests/regression/run_regression.sbatch v0.9.0 "$REG" --write-expected
```

For each release candidate:

```bash
sbatch tests/regression/run_regression.sbatch <tag-or-commit> "$REG"
```

This runs the candidate in `$REG/runs/<tag>/` and then runs
`tests/compare_outputs.py` against `$REG/expected/`, with the same tolerances
as CI. Attach `$REG/runs/<tag>/comparison.txt` to the release PR. After the
release is merged, rerun with `--write-expected` so the next comparison is
against it.
