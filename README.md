# RNAseq_Pipeline

This directory contains a Snakemake workflow for processing bulk RNA-seq data. The pipeline automates quality control, trimming, alignment, quantification, and summarization for multiple samples.

## Workflow Steps

1. **FastQC**: Quality control of raw FASTQ files.
2. **Trimmomatic**: Adapter and quality trimming of reads.
3. **HISAT2**: Alignment of trimmed reads to a reference genome.
4. **Samtools**: Sorting and indexing of BAM files.
5. **featureCounts**: Gene-level quantification for all samples (single matrix output).
6. **Salmon**: Transcript-level quantification.
7. **MultiQC**: Aggregated report of QC and quantification results.

![The workflow](rulegraph.png)

## Directory Structure
- `Snakefile`: Main workflow definition.
- `config.yaml`: Configuration file with paths and parameters.
- `submit_snakemake.sh`: Script to submit the workflow to a cluster.
- `data/`: Raw FASTQ files and related data.
- `fastqc/`: FastQC output files.
- `logs/`: Log files for each step.
- `results/`: Processed data outputs (feature counts, alignments, quantifications).
- `multiqc_data/`: MultiQC intermediate files.
- `multiqc_report.html`: Final MultiQC report.

## Installation

```bash
git clone https://github.com/uci-grthub/uci_grthub_rnaseq_pipeline.git
cd uci_grthub_rnaseq_pipeline
```

## Usage

### 1. Prerequisites

Make sure you have Snakemake installed. You can install it using conda:

```bash
conda install -c conda-forge -c bioconda snakemake
```

Or use the pinned [pixi](https://pixi.sh) environment in `pixi.toml`, which installs
Snakemake + `snakefmt` without touching your system/conda setup:

```bash
pixi install -e snakemake-dfs
pixi run -e snakemake-dfs snakemake --cores 8
```

### 2. Configuration
 Edit `config.yaml` to set paths and parameters for your data and references.
- Sample names
- Input/output paths
- Reference file locations
- Tool parameters

### 3. Running the workflow

#### Option A: Submit to SLURM cluster
```bash
sbatch submit_snakemake.sh
```

#### Option B: Run locally (for testing)
```bash
snakemake --cores 8 --use-conda
```

#### Option C: Dry run (to check workflow)
```bash
snakemake --dry-run
```

### 4. Workflow visualization

Generate a workflow diagram:
```bash
snakemake --dag | dot -Tpng > workflow.png
```
### 5. Output
1. **FastQC**: Quality control reports in `fastqc/`.
3. **Results**: Outputs will be found in the `results/` and other specified directories. The main featureCounts output is `results/feature_count/all_samples_counts.txt`.

## Requirements
- Snakemake
- Modules: fastqc, trimmomatic, hisat2, samtools, subread, salmon, singularity
- Cluster environment (recommended)

## Customization
- Adjust sample detection, references, and tool parameters in `config.yaml`.
- Modify `cluster.yaml` for resource allocation.

### For different library types:
- **Non-stranded libraries**: Change `rna_strandness` to "unstranded" and `library_type` to "IU" in `config.yaml`
- **Different strand orientation**: Modify the strandness parameters accordingly

## Key Differences from Original Script

1. **Modular design**: Each step is a separate rule
2. **Dependency management**: Snakemake automatically handles job dependencies
3. **Parallel execution**: Multiple samples can be processed simultaneously
4. **Configuration-driven**: Easy to modify parameters without editing the main workflow
5. **Resource management**: Better integration with SLURM scheduler
6. **Reproducibility**: Workflow tracks input/output dependencies

## Troubleshooting

1. Check SLURM job status: `squeue -u $USER`
2. View workflow status: `snakemake --summary`
3. Check individual rule logs in the SLURM output files

# TODO
1. retrieve counts from featureCounts and Salmon quantification files, and summarize them in a final report
2. Run DESeq2 or edgeR for differential expression analysis


## Contact
For questions or issues, contact: kstachel@uci.edu
