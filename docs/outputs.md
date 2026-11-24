# Outputs

All paths are relative to the project directory. `output/` is `paths.output`.

| Deliverable | Location | Default run |
| --- | --- | --- |
| Raw read QC | `output/fastqc/<sample>/` | yes |
| Trimmed reads, chelae metrics | `output/trimmed/` | yes |
| Alignments (one BAM per library), HISAT2 summaries | `output/hisat2_alignment/` | yes |
| Duplicate-marking metrics | `output/hisat2_alignment/<library>_markdup_metrics.txt` | yes |
| Post-alignment RNA QC | `output/rustqc/<library>/` | yes |
| featureCounts table | `output/feature_count/<species>_samples_counts.txt` | yes |
| Gene count matrix (raw counts, CPM, annotation) | `output/counts/<species>/` | yes |
| Correlation, clustering, PCA | `output/sample_qc/<species>/` | yes |
| Salmon transcript quantification | `output/salmon/<library>_salmon_quant/` | yes |
| Gene-level TPM, transcript counts/TPM | `output/tpm/<species>/` | yes |
| Aggregated QC report | `output/multiqc_report.html` | yes |
| Resolved run config | `output/pipeline_config.yaml` | yes |
| Pipeline/tool/reference versions | `output/versions.yaml` | yes |
| Project report | `RNAseq_Project_Report.pdf` | yes |
| Alternative splicing | `output/rmats/<species>/` | on request |
| GEO/SRA submission package | `output/ncbi_submission/<species>/` | on request |
| DESeq2 results | `output/deseq2/<species>/` | on request |

Logs are in `logs/<rule>/` and run-time benchmarks in `benchmarks/<rule>/`.

## versions.yaml

```yaml
pipeline: {name, version, commit}        # commit ends in -dirty for local edits
software:
  source: pixi | env-modules
  pixi_lock_sha256: ...
  tools: {fastqc: ..., hisat2: ..., ...}
references:
  <species>: {genome_build, annotation, hisat2_index, gtf, gtf_sha256, salmon_index, ...}
written: <UTC timestamp>
```

Put the pipeline version and DOI from this file in the delivery notes.

## NCBI submission package

`output/ncbi_submission/<species>/` contains:

- `geo_samples.csv`: paste into the SAMPLES section of the GEO metadata workbook.
- `sra_metadata.csv`
- `md5sums.txt`: checksums for every raw and processed file. Raw FASTQ checksums are reused from `<paths.data>/md5sums.txt` when the sequencing core supplied one.
- `SUBMISSION_README.txt`: what to upload.

Tissue and instrument come from the `ncbi` section of the config.
