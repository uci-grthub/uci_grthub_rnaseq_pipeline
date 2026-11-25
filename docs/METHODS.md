# Methods text

Paste the paragraph below into a methods section, and update it for any
parameters the project changed (`output/pipeline_config.yaml`). Get the exact
versions from the run's `output/versions.yaml`. This file is updated with
every release.

## uci_grthub_rnaseq_pipeline v1.0.0

Bulk RNA-seq data were processed with uci_grthub_rnaseq_pipeline v1.0.0
(UCI Genomics Research and Technology Hub; DOI: _to be added on release_),
a Snakemake (v9.16.3) workflow. Raw read quality was assessed with FastQC
v0.11.9. Adapters and low-quality bases were removed with chelae v0.2.1
(TruSeq adapter preset, 4-base sliding window with a minimum quality of 15,
minimum read length 20), and read pairs were kept or discarded together.
Trimmed reads were aligned to the reference genome with HISAT2 v2.2.1
(`--rna-strandness RF --dta-cufflinks`). Lanes of the same library were
aligned together. Alignments were sorted, mate-fixed and duplicate-marked with
samtools v1.15.1. Gene-level counts were produced with featureCounts (Subread
v2.0.1) on exon features summarised by gene ID, with paired-end fragment
counting and reverse-strand specificity (`-s 2`). Transcript abundance was
estimated with Salmon (v1.8.0 with the HPC3 modules, v1.10.3 with the pixi
environment; `output/versions.yaml` records which) (selective alignment, `--validateMappings
--gcBias`, library type ISR) and summarised to gene level with tximport
v1.38 in R v4.5. Sample-level QC included correlation, hierarchical
clustering and principal component analysis on variance-stabilised counts
(DESeq2 v1.50), and post-alignment RNA QC with RustQC. QC results were
aggregated with MultiQC v1.20.

| Species key | Genome | Annotation |
| --- | --- | --- |
| `human` | GRCh38 | Ensembl 84 |
| `mouse` | GRCm38 (mm10) | GENCODE vM24 |
| `macaque` | Mmul_10 (rheMac10) | Ensembl 116 |
| `macaque_t2t` | T2T-MMU8v2.0 (GCF_049350105.2) | NCBI RefSeq RS_2025_08 |

Optional steps, when they were run:
- Alternative splicing was quantified with rMATS v4.3.0.
- Differential expression was tested with DESeq2 v1.50.

## Tool citations

- Snakemake: Mölder F. et al. Sustainable data analysis with Snakemake. F1000Research 10:33 (2021).
- FastQC: Andrews S. FastQC (Babraham Bioinformatics).
- HISAT2: Kim D. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37:907–915 (2019).
- SAMtools: Danecek P. et al. Twelve years of SAMtools and BCFtools. GigaScience 10:giab008 (2021).
- featureCounts: Liao Y. et al. featureCounts. Bioinformatics 30:923–930 (2014).
- Salmon: Patro R. et al. Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods 14:417–419 (2017).
- tximport: Soneson C. et al. Differential analyses for RNA-seq. F1000Research 4:1521 (2015).
- DESeq2: Love M.I. et al. Genome Biol 15:550 (2014).
- rMATS: Shen S. et al. rMATS. PNAS 111:E5593–E5601 (2014).
- MultiQC: Ewels P. et al. MultiQC. Bioinformatics 32:3047–3048 (2016).
