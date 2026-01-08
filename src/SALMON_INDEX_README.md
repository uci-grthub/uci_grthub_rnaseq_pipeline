# Salmon Decoy-Aware Index Building Guide

## Overview

This script builds a **decoy-aware salmon index** for Drosophila melanogaster RNA-seq quantification. A decoy-aware index improves mapping accuracy by reducing false positive mappings to intronic/genomic regions.

## What is a Decoy-Aware Index?

- **Transcriptome**: Contains all transcript sequences (what we want to quantify)
- **Decoys**: Genomic sequences that are NOT part of the transcriptome
- **Benefit**: Reads that would falsely map to intronic regions map to decoys instead, reducing bias

## Files Required

1. **Transcriptome**: `dmel-all-transcript-r6.53.fasta.gz`
   - Contains all known transcripts for D. melanogaster
   
2. **Genome**: `dmel-all-chromosome-r6.53.fasta`
   - Chromosome sequences used as decoys

Both files should be in the `data/` directory.

## Script Workflow

The `build_salmon_index.sh` script performs these steps:

1. **Validate files**: Checks that transcriptome and genome files exist
2. **Create directories**: Sets up temporary and output directories
3. **Decompress transcriptome**: Extracts the gzipped transcriptome file
4. **Generate decoys file**: Extracts sequence headers from the genome
5. **Create gentrome**: Concatenates transcriptome + genome sequences
6. **Build index**: Runs `salmon index` with decoy information

## Usage

```bash
# Navigate to project directory
cd /dfs9/ucightf-lab/projects/JafaM/251216_JafaM_RNAseq_24plex

# Run the script
bash build_salmon_index.sh
```

### With a Job Scheduler (Slurm Example)

```bash
#!/bin/bash
#SBATCH --job-name=salmon_index
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=02:00:00

cd /dfs9/ucightf-lab/projects/JafaM/251216_JafaM_RNAseq_24plex
bash build_salmon_index.sh
```

## Output

The script generates:

- **`salmon_index/`**: The final salmon index (used for quantification)
- **`salmon_temp/`**: Temporary files including:
  - `dmel-all-transcript-r6.53.fasta`: Decompressed transcriptome
  - `decoys.txt`: List of genomic sequence names
  - `gentrome.fasta`: Combined transcriptome + genome

## Using the Index for Quantification

Once the index is built, quantify RNA-seq samples:

```bash
# Single-end reads
salmon quant -i salmon_index -l A -r sample.fq -o output_quant

# Paired-end reads
salmon quant -i salmon_index -l A -1 sample_R1.fq -2 sample_R2.fq -o output_quant
```

## Customization

To modify the script, edit these variables:

```bash
DATA_DIR="${SCRIPT_DIR}/data"           # Location of reference files
INDEX_DIR="${SCRIPT_DIR}/salmon_index"  # Where to save index
TEMP_DIR="${SCRIPT_DIR}/salmon_temp"    # Temporary files location
```

Also adjust the number of threads in the salmon command:
```bash
--threads 8  # Change to match your available CPUs
```

## Requirements

- **salmon**: Must be installed and in PATH
- **gunzip**: For decompressing files (standard Unix utility)
- **disk space**: ~10-15 GB for index and temporary files
- **memory**: 16-32 GB recommended

## Troubleshooting

**Error: Transcriptome file not found**
- Verify `data/dmel-all-transcript-r6.53.fasta.gz` exists
- Check file permissions

**Error: Genome file not found**
- Verify `data/dmel-all-chromosome-r6.53.fasta` exists
- May need to decompress: `gunzip data/dmel-all-chromosome-r6.53.fasta.gz`

**Index build fails**
- Check available disk space: `df -h`
- Check available memory: `free -h`
- Ensure salmon is properly installed: `salmon --version`

## References

- [Salmon Documentation](https://salmon.readthedocs.io/)
- [Decoy Sequences for Alignment](https://salmon.readthedocs.io/en/latest/salmon.html#using-decoys)
- [FlyBase D. melanogaster Genome](https://flybase.org/)
