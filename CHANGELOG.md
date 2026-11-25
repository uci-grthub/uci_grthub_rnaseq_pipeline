# Changelog

All notable changes to this pipeline. From v1.0.0 on, release-please
maintains this file from Conventional Commit titles.

The version number says whether results are comparable:
- **major**: results change for the same input.
- **minor**: new optional features; default results are unchanged.
- **patch**: fixes with no change to results.

## [Unreleased] — v1.0.0 (release candidate)

Packaging and process release. On the same inputs it produces the same
counts, TPMs (r > 0.999), QC metrics and rMATS events as v0.9.0. See the
regression summary on the release PR.

### Changed
- The repository follows the Snakemake workflow catalog layout (`workflow/`,
  `config/`, `profiles/`, `tests/`, `docs/`). Projects now run from their own
  directory with `--directory` and `--configfile`; see docs/usage.md for the
  key renames (`deseq2.metadata` → `samples`, `config.local.yaml` → project
  `config.yaml`).
- Tools come from pinned pixi environments (`default`, `tools`, `r`) instead of
  `module load` inside rules. The versions are the same as the v0.9.0 modules,
  except R, which is 4.5.3 instead of 4.5.2 (same DESeq2/tximport versions),
  and Salmon, which is 1.10.3 instead of 1.8.0 because the bioconda 1.8.0
  build segfaults on HPC3 (provisional; see the PR). MultiQC 1.20 needs
  `rich<14`, `plotly<6` and `python-kaleido 0.2.1`, which are pinned.
  `--software-deployment-method env-modules` still uses the HPC3 modules,
  including salmon/1.8.0.
- MultiQC 1.20 comes from pixi instead of a hard-coded .sif path. RustQC uses
  `tools.rustqc_container` when that file exists (HPC3), otherwise pixi rustqc
  0.2.1. RustQC metrics are not part of the regression comparison.
- `submit_snakemake.sh` is replaced by `profiles/hpc3` and `scripts/run_hpc3.sbatch`.
  The Slurm account comes from `SLURM_ACCOUNT` instead of being set per rule.
- Completion e-mails go to `notify_email`, which is off by default; the
  address is no longer hard-coded.
- `ncbi.tissue` defaults to empty and is set per project.

### Added
- Config, reference and sample-sheet schemas, checked at startup.
- `output/versions.yaml` (pipeline version, commit, pixi lock hash, tool
  versions, references) and `output/pipeline_config.yaml` on every run.
- `params.rmats.read_length` (default 150, as before).
- Test set (public PRJNA929683 reads on a 10 Mb chr22 window),
  `tests/compare_outputs.py`, GitHub Actions CI (lint + test run), release-please,
  MIT license, CITATION.cff, docs/METHODS.md.

### Removed
- Committed run artifacts and placeholders (`data/`, `logs/`, `output/`,
  `results/`, `metadata/`, `proj_src/`, `.specstory/`), `environment_simple.yml`,
  and unused config keys from the unmerged limma-voom work
  (`sex_qc`, `limma_voom`, `hisat2.min_alignment_rate`).

## v0.9.0 — 2026-11-24

Pre-modernization baseline: `main` as of commit 4b00554, tagged so later
versions can be compared against it.
