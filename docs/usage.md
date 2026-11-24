# Running a project

The pipeline is installed once per release and run from a separate project
directory. The release is never edited for a project.

## 1. Install a release (once per version)

```bash
cd /dfs9/ucightf-lab/$USER/pipelines        # anywhere shared and durable
git clone --branch v1.0.0 --depth 1 \
    https://github.com/uci-grthub/uci_grthub_rnaseq_pipeline uci_grthub_rnaseq_pipeline-v1.0.0
cd uci_grthub_rnaseq_pipeline-v1.0.0
pixi install --all                          # default, tools and r environments
```

Production runs use a tagged release only, never `main`. A project started on
one major version finishes on it.

## 2. Set up the project directory

```
my_project/
├── config.yaml        # only the keys that differ from config/config.yaml
├── metadata/metadata.csv   # sample sheet (path set by `samples`)
└── data/FASTQ/        # <sample>-R1.fastq.gz / <sample>-R2.fastq.gz
```

A minimal `config.yaml`:

```yaml
default_species: mouse          # or fill the species column per sample
ncbi:
  tissue: "hippocampus"
notify_email: "you@uci.edu"     # optional
```

Every key and its default is in `config/config.yaml`, and the allowed values
are in `workflow/schemas/config.schema.yaml`. Reference genomes are keyed by
species in `config/references.yaml`. To add a project-only reference, put a
`species_references:` entry in the project config.yaml.

The sample sheet is a CSV with at least a `sample` column, plus `species`,
`condition`, `sex`, the barcode columns `i7barcode`/`i5barcode_NovaSeqV1.5`,
and optional `library` (rows sharing a library are aligned together) and
`sample_id`. If the file is missing, the first run writes one, filled in from
the FASTQ names. Fill in species/condition before running for real.

FASTQs can be fetched with `sbatch <release>/scripts/rsync_fastq.sh SRC_DIR`.

## 3. Run on HPC3

From the project directory:

```bash
export SLURM_ACCOUNT=sbsandme_lab   # your lab account
sbatch -A $SLURM_ACCOUNT <release>/scripts/run_hpc3.sbatch -n   # dry run
sbatch -A $SLURM_ACCOUNT <release>/scripts/run_hpc3.sbatch      # real run
```

The controller job submits each rule as its own Slurm job through
`profiles/hpc3`. The equivalent direct command is:

```bash
pixi run --manifest-path <release>/pixi.toml snakemake \
    -s <release>/workflow/Snakefile --directory . --configfile config.yaml \
    --profile <release>/profiles/hpc3
```

### Software: pixi or HPC3 modules

By default every tool comes from the release's pixi environments. Pass
`--software-deployment-method env-modules` to use the HPC3 modules instead,
for the rules that declare them (FastQC, HISAT2, samtools, featureCounts,
Salmon, rMATS, R). `output/versions.yaml` records which source was used.

RustQC runs from the image at `tools.rustqc_container` when that file exists.
HPC3's glibc is too old for the bioconda build, so on HPC3 the image is used.

## 4. Optional targets

These are not in the default run. Request them by name:

```bash
.../run_hpc3.sbatch all output/rmats/<species>/.done
.../run_hpc3.sbatch all output/ncbi_submission/<species>/geo_samples.csv
.../run_hpc3.sbatch output/deseq2/<species>/deseq2_results.csv   # needs de_comparisons.yaml
```

## Migrating a project that ran inside the repo (pre-1.0)

| Before | Now |
| --- | --- |
| `config.local.yaml` in the repo | `config.yaml` in the project dir, passed with `--configfile` |
| `config.yaml` → `deseq2.metadata` | `samples` |
| `deseq2.infer_metadata_from_fastq` | `infer_samples_from_fastq` |
| `config.species_references.yaml` | `config/references.yaml` (in the release) |
| `src/de_comparisons.yaml` | `de_comparisons.yaml` in the project dir (`deseq2.comparisons_config`) |
| `submit_snakemake.sh` | `scripts/run_hpc3.sbatch` + `profiles/hpc3` |
| hard-coded completion e-mail | `notify_email` (off by default) |

Output paths and the report name did not change.

Because v1.0 changes rule code (tool activation lines), Snakemake reruns
every rule of a project that was started on v0.9 in the same directory. Finish
running projects on the version they started on.
