EMAIL = "kstachel@uci.edu"


onstart:
    shell("mail -s 'STARTED' {EMAIL} < {log}")


onsuccess:
    shell("mail -s 'DONE' {EMAIL} < {log}")


onerror:
    shell("mail -s 'ERROR' {EMAIL} < {log}")


# Snakemake workflow for RNA-seq analysis
# Generalized to process multiple samples from a data directory

import csv
import glob
import os
import re

from snakemake.exceptions import WorkflowError


# Load configuration
configfile: "config.yaml"


configfile: "config.species_references.yaml"


_local_config = "config.local.yaml"
if os.path.exists(_local_config):

    configfile: _local_config


# Extract configuration variables
DATA_PATH = config["paths"]["data"]
OUTPUT_DIR = config["paths"]["output"]

# Get sample names from FASTQ files in data directory
# Support paired-end files named like "*_r1.fq.gz"/"*_r2.fq.gz",
# "*_R1.fq.gz"/"*_R2.fq.gz", "*-R1.fastq.gz"/"*-R2.fastq.gz",
# and Illumina output named like "<sample>-READ1-Sequences.txt.gz"/"<sample>-READ2-Sequences.txt.gz".
suffixes_all = [
    "_r1.fq.gz",
    "_r2.fq.gz",
    "_R1.fq.gz",
    "_R2.fq.gz",
    "-R1.fastq.gz",
    "-R2.fastq.gz",
    "-READ1-Sequences.txt.gz",
    "-READ2-Sequences.txt.gz",
]

# Collect samples. Handle two common layouts:
# 1) data/FASTQ/<sample>/*-READ1-Sequences.txt.gz  (each sample in its own folder)
# 2) data/FASTQ/*_r1.fq.gz (files directly under DATA_PATH)

sample_set = set()
try:
    entries = os.listdir(DATA_PATH)
except Exception:
    entries = []

for entry in entries:
    path = os.path.join(DATA_PATH, entry)
    # If entry is a directory, check whether it contains read files with expected suffixes
    if os.path.isdir(path):
        for suf in suffixes_all:
            matches = glob.glob(os.path.join(path, f"*{suf}"))
            if matches:
                sample_set.add(entry)
                break
    else:
        # entry is a file directly under DATA_PATH
        base = os.path.basename(entry)
        for suf in suffixes_all:
            if base.endswith(suf):
                sample_set.add(base[: -len(suf)])
                break

# If nothing found yet, fall back to a recursive glob search for common R1 patterns
if not sample_set:
    matches = []
    for pat in ["*_r1.fq.gz", "*_R1.fq.gz", "*-R1.fastq.gz"]:
        matches += glob.glob(os.path.join(DATA_PATH, "**", pat), recursive=True)
    for fn in matches:
        base = os.path.basename(fn)
        for suf in suffixes_all:
            if base.endswith(suf):
                sample_set.add(base[: -len(suf)])
                break

SAMPLES = sorted(sample_set)
print(f"Found {len(SAMPLES)} samples: {SAMPLES}", file=sys.stderr)


def _infer_barcodes_from_sample(sample_name):
    # Typical format: xR074-L8-G3-P057-ATGTACCT-TAGGTATG -> last two parts are i7/i5
    parts = sample_name.split("-")
    i7 = parts[-2] if len(parts) >= 2 else ""
    i5 = parts[-1] if len(parts) >= 1 else ""

    def clean(idx):
        return idx if re.fullmatch(r"[ACGTN]+", idx or "") and len(idx) >= 6 else ""

    return clean(i7), clean(i5)


# Create a default metadata file if none exists yet, so that rules depending
# on it (generate_report, deseq2) resolve during DAG building, including
# `snakemake -n`. When infer_metadata_from_fastq is true (default), rows are
# pre-populated with sample names and barcodes inferred from FASTQ filenames;
# sex/condition are left blank for the user to fill in before a real run.
METADATA_PATH = config["deseq2"]["metadata"]
INFER_METADATA_FROM_FASTQ = config["deseq2"].get("infer_metadata_from_fastq", True)
if not os.path.exists(METADATA_PATH):
    os.makedirs(os.path.dirname(METADATA_PATH) or ".", exist_ok=True)
    with open(METADATA_PATH, "w") as fh:
        fh.write("sample,i7barcode,i5barcode_NovaSeqV1.5,sex,condition,species\n")
        if INFER_METADATA_FROM_FASTQ:
            for sample in SAMPLES:
                i7, i5 = _infer_barcodes_from_sample(sample)
                fh.write(f"{sample},{i7},{i5},,,\n")

# Reference paths
ADAPTER_PATH = config["references"]["adapters"]
HISAT2_INDEX = config["references"]["hisat2_index"]
GTF_PATH = config["references"]["gtf"]
SALMON_INDEX = config["references"]["salmon_index"]
TRIMMOMATIC_JAR = config["tools"]["trimmomatic"]
RUSTQC_CONTAINER = "/dfs9/ucightf-lab/kstachel/containers/rustqc.sif"

# Per-sample species -> reference paths, so samples from different organisms
# (e.g. mouse libraries mixed into an otherwise-human run) get aligned and
# quantified against the correct genome instead of all using HISAT2_INDEX above.
#
# A blank species is a hard error, never a silent fall back to human. The
# auto-inferred metadata stub cannot know the organism, so it leaves the column
# empty; defaulting that to human sent all 27 xR106-L5-G3 libraries through the
# human GRCh38 index at a 1.7% alignment rate (they were mouse: 98.9% against
# GRCm38) and still produced a perfectly well-formed counts matrix. The same
# thing happened earlier to xR098-L1-G5-P070/P071. Nothing downstream can
# detect it, so refuse to build the DAG instead -- before any compute is spent.
SPECIES_REFERENCES = config["species_references"]

SAMPLE_SPECIES = {}
_metadata_samples = set()
_species_problems = []
with open(METADATA_PATH, newline="") as fh:
    for row in csv.DictReader(fh):
        sample = (row.get("sample") or "").strip()
        if not sample:
            continue
        _metadata_samples.add(sample)
        species = (row.get("species") or "").strip().lower()
        if not species:
            _species_problems.append(f"  {sample}: 'species' is blank")
        elif species not in SPECIES_REFERENCES:
            _species_problems.append(f"  {sample}: unknown species {species!r}")
        else:
            SAMPLE_SPECIES[sample] = species

for _sample in SAMPLES:
    if _sample not in _metadata_samples:
        _species_problems.append(f"  {_sample}: no row in {METADATA_PATH}")

if _species_problems:
    raise WorkflowError(
        "Cannot resolve a reference genome for every sample:\n"
        + "\n".join(sorted(_species_problems))
        + f"\n\nFill the 'species' column in {METADATA_PATH} using one of: "
        + ", ".join(sorted(SPECIES_REFERENCES))
        + ".\nThere is deliberately no default -- guessing wrong aligns the run "
        "against the wrong genome and still yields a plausible-looking counts "
        "matrix."
    )

SPECIES_LIST = sorted({SAMPLE_SPECIES[sample] for sample in SAMPLES})


def species_ref(sample, key):
    return SPECIES_REFERENCES[SAMPLE_SPECIES[sample]][key]


# Rule all - defines final outputs
rule all:
    input:
        # # FastQC reports
        expand(
            f"{OUTPUT_DIR}/fastqc/{{sample}}/{{sample}}-R1_fastqc.html", sample=SAMPLES
        ),
        expand(
            f"{OUTPUT_DIR}/fastqc/{{sample}}/{{sample}}-R2_fastqc.html", sample=SAMPLES
        ),
        # RustQC markers
        expand(f"{OUTPUT_DIR}/rustqc/{{sample}}/.done", sample=SAMPLES),
        # Trimmed files
        expand(f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_1P.fq.gz", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_2P.fq.gz", sample=SAMPLES),
        # HISAT2 alignment and counting
        expand(
            f"{OUTPUT_DIR}/hisat2_alignment/{{sample}}_align_sorted_markdup.bam",
            sample=SAMPLES,
        ),
        expand(
            f"{OUTPUT_DIR}/hisat2_alignment/{{sample}}_align_sorted_markdup.bam.bai",
            sample=SAMPLES,
        ),
        expand(
            f"{OUTPUT_DIR}/feature_count/{{species}}_samples_counts.txt",
            species=SPECIES_LIST,
        ),
        # Clean gene-level raw count matrix (deliverable form of featureCounts)
        expand(f"{OUTPUT_DIR}/counts/{{species}}/gene_counts.csv", species=SPECIES_LIST),
        # MultiQC report
        f"{OUTPUT_DIR}/multiqc_report.html",
        # limma-voom differential expression (the method used for this project)
        expand(
            f"{OUTPUT_DIR}/limma_voom/{{species}}/limma_voom_results.csv",
            species=SPECIES_LIST,
        ),
        # Project report
        "RNAseq_Project_Report.pdf",
        # DESeq2 results not included by default -- run on demand with e.g.
        # `snakemake output/deseq2/mouse/deseq2_results.csv`


# Rule 0: FastQC on raw FASTQ files
rule fastqc:
    input:
        r1=f"{DATA_PATH}/{{sample}}-R1.fastq.gz",
        r2=f"{DATA_PATH}/{{sample}}-R2.fastq.gz",
    output:
        r1_html=f"{OUTPUT_DIR}/fastqc/{{sample}}/{{sample}}-R1_fastqc.html",
        r1_zip=f"{OUTPUT_DIR}/fastqc/{{sample}}/{{sample}}-R1_fastqc.zip",
        r2_html=f"{OUTPUT_DIR}/fastqc/{{sample}}/{{sample}}-R2_fastqc.html",
        r2_zip=f"{OUTPUT_DIR}/fastqc/{{sample}}/{{sample}}-R2_fastqc.zip",
    threads: 2
    resources:
        mem_mb=4000,
        cpus=2,
        partition="standard",
        account="sbsandme_lab",
    params:
        out_dir=f"{OUTPUT_DIR}/fastqc/{{sample}}",
    log:
        "logs/fastqc/{sample}.log",
    benchmark:
        "benchmarks/fastqc/{sample}.tsv"
    shell:
        """
        exec > {log} 2>&1
        module load fastqc/0.11.9
        rm -rf {params.out_dir}
        mkdir -p {params.out_dir}
        fastqc -o {params.out_dir} -t {threads} {input.r1} {input.r2}
        module unload fastqc/0.11.9
        """


# Rule 0b: RustQC RNA QC on aligned BAM files via singularity
rule rustqc:
    input:
        bam=f"{OUTPUT_DIR}/hisat2_alignment/{{sample}}_align_sorted_markdup.bam",
    output:
        done=f"{OUTPUT_DIR}/rustqc/{{sample}}/.done",
    threads: 8
    resources:
        mem_mb=24000,
        cpus=config["params"]["cpus"],
        partition="standard",
        account="sbsandme_lab",
    params:
        out_dir=f"{OUTPUT_DIR}/rustqc/{{sample}}",
        gtf_path=lambda wildcards: species_ref(wildcards.sample, "gtf"),
    log:
        "logs/rustqc/{sample}.log",
    benchmark:
        "benchmarks/rustqc/{sample}.tsv"
    shell:
        """
        exec > {log} 2>&1
        module load singularity/3.11.3
        mkdir -p {params.out_dir}
        singularity exec {RUSTQC_CONTAINER} rustqc rna \
            --gtf {params.gtf_path} \
            --outdir {params.out_dir} \
            --threads {threads} \
            --paired \
            --skip-dup-check \
            {input.bam}
        touch {output.done}
        module unload singularity/3.11.3
        """


# Rule 1: Trimming with Trimmomatic
rule trimmomatic:
    input:
        r1=f"{DATA_PATH}/{{sample}}-R1.fastq.gz",
        r2=f"{DATA_PATH}/{{sample}}-R2.fastq.gz",
    output:
        r1_paired=f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_1P.fq.gz",
        r1_unpaired=f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_1U.fq.gz",
        r2_paired=f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_2P.fq.gz",
        r2_unpaired=f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_2U.fq.gz",
    threads: 8
    resources:
        mem_mb=4000,
        cpus=config["params"]["cpus"],
        partition="standard",
        account="sbsandme_lab",
    params:
        adapter_path=ADAPTER_PATH,
        trimmed_base=f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed.fq.gz",
    log:
        f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmomatic.log",
    benchmark:
        "benchmarks/trimmomatic/{sample}.tsv"
    shell:
        """
        exec > {log} 2>&1
        java -jar {TRIMMOMATIC_JAR} PE \
            -threads {threads} -phred33 \
            -baseout {params.trimmed_base} \
            {input.r1} {input.r2} \
            ILLUMINACLIP:{params.adapter_path}:{config[params][trimmomatic][illuminaclip]} \
            SLIDINGWINDOW:{config[params][trimmomatic][sliding_window]} \
            MINLEN:{config[params][trimmomatic][min_length]}
        """


# Rule 2: HISAT2 alignment
rule hisat2_align:
    input:
        r1=f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_1P.fq.gz",
        r2=f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_2P.fq.gz",
    output:
        bam=f"{OUTPUT_DIR}/hisat2_alignment/{{sample}}_align.bam",
        summary=f"{OUTPUT_DIR}/hisat2_alignment/alignment_summary/{{sample}}_summary.align",
    threads: 8
    resources:
        mem_mb=24000,
        cpus=config["params"]["cpus"],
        partition="standard",
        account="sbsandme_lab",
    params:
        hisat2_index=lambda wildcards: species_ref(wildcards.sample, "hisat2_index"),
        summary_path=f"{OUTPUT_DIR}/hisat2_alignment/alignment_summary",
        min_alignment_rate=config["params"]["hisat2"]["min_alignment_rate"],
    log:
        "logs/hisat2_align/{sample}.log",
    benchmark:
        "benchmarks/hisat2_align/{sample}.tsv"
    shell:
        """
        exec 2> {log}
        module load hisat2/2.2.1
        module load samtools/1.15.1

        hisat2 -p {threads} -t --qc-filter --rna-strandness {config[params][hisat2][rna_strandness]} \
            --summary-file {output.summary} \
            -x {params.hisat2_index} --dta-cufflinks \
            -1 {input.r1} -2 {input.r2} \
            | samtools sort -n -@ 2 \
            | samtools fixmate -m -@ 2 - {output.bam}

        # Aligning against the wrong genome is not an error condition to hisat2:
        # it exits 0 and emits a valid, sorted, well-formed BAM that simply has
        # almost nothing in it. featureCounts then happily builds a counts
        # matrix out of the residue. Gate on the rate here so the failure
        # surfaces on the first sample instead of in a DE result nobody can
        # interpret. An unparseable summary fails too -- that means hisat2 died
        # partway and the BAM is truncated.
        rate=$(grep -oE '[0-9.]+% overall alignment rate' {output.summary} \
            | grep -oE '[0-9.]+' || true)
        if [ -z "$rate" ]; then
            echo "ERROR: no alignment rate in {output.summary}; hisat2 did not finish" >&2
            exit 1
        fi
        min_rate={params.min_alignment_rate}
        if [ "$(awk -v r="$rate" -v m="$min_rate" 'BEGIN {{ print (r+0 < m+0) ? "low" : "ok" }}')" = "low" ]; then
            echo "ERROR: {wildcards.sample} aligned at ${{rate}}% against {params.hisat2_index}," >&2
            echo "       below the ${{min_rate}}% minimum. This is what the wrong reference genome" >&2
            echo "       looks like -- check the 'species' column in the metadata for this sample." >&2
            exit 1
        fi
        echo "{wildcards.sample}: ${{rate}}% overall alignment rate (minimum ${{min_rate}}%)"

        module unload samtools/1.15.1
        module unload hisat2/2.2.1
        """


# Rule 3: Sort and index BAM file
rule sort_bam:
    input:
        bam=f"{OUTPUT_DIR}/hisat2_alignment/{{sample}}_align.bam",
    output:
        sorted_bam=f"{OUTPUT_DIR}/hisat2_alignment/{{sample}}_align_sorted.bam",
        index=f"{OUTPUT_DIR}/hisat2_alignment/{{sample}}_align_sorted.bam.bai",
    threads: 8
    resources:
        mem_mb=24000,
        cpus=config["params"]["cpus"],
        partition="standard",
        account="sbsandme_lab",
    log:
        "logs/sort_bam/{sample}.log",
    benchmark:
        "benchmarks/sort_bam/{sample}.tsv"
    shell:
        """
        exec > {log} 2>&1
        module load samtools/1.15.1

        samtools sort -@ {threads} -o {output.sorted_bam} {input.bam}
        samtools index -@ {threads} {output.sorted_bam}

        module unload samtools/1.15.1
        """


# Rule 3b: Mark duplicates
rule markdup:
    input:
        sorted_bam=f"{OUTPUT_DIR}/hisat2_alignment/{{sample}}_align_sorted.bam",
    output:
        markdup_bam=f"{OUTPUT_DIR}/hisat2_alignment/{{sample}}_align_sorted_markdup.bam",
        markdup_bai=f"{OUTPUT_DIR}/hisat2_alignment/{{sample}}_align_sorted_markdup.bam.bai",
        metrics=f"{OUTPUT_DIR}/hisat2_alignment/{{sample}}_markdup_metrics.txt",
    threads: 8
    resources:
        mem_mb=24000,
        cpus=config["params"]["cpus"],
        partition="standard",
        account="sbsandme_lab",
    log:
        "logs/markdup/{sample}.log",
    benchmark:
        "benchmarks/markdup/{sample}.tsv"
    shell:
        """
        exec > {log} 2>&1
        module load samtools/1.15.1

        samtools markdup -@ {threads} -f {output.metrics} {input.sorted_bam} {output.markdup_bam}
        samtools index -@ {threads} {output.markdup_bam}

        module unload samtools/1.15.1
        """


# Rule 4: Feature counting (all samples of a given species together --
# featureCounts takes a single -a GTF, so samples from different organisms
# cannot be combined into one run and are split by species instead)
rule feature_counts_all:
    input:
        bam_files=lambda wildcards: expand(
            f"{OUTPUT_DIR}/hisat2_alignment/{{sample}}_align_sorted_markdup.bam",
            sample=[
                sample
                for sample in SAMPLES
                if SAMPLE_SPECIES[sample] == wildcards.species
            ],
        ),
    output:
        counts=f"{OUTPUT_DIR}/feature_count/{{species}}_samples_counts.txt",
    threads: 4
    resources:
        mem_mb=24000,
        cpus=4,
        partition="standard",
        account="sbsandme_lab",
    params:
        gtf_path=lambda wildcards: SPECIES_REFERENCES[wildcards.species]["gtf"],
    log:
        "logs/feature_counts_all/{species}.log",
    benchmark:
        "benchmarks/feature_counts_all/{species}.tsv"
    shell:
        """
        exec > {log} 2>&1
        module load subread/2.0.1
        featureCounts -s {config[params][feature_counts][strandness]} -p -t exon -g gene_id -T {threads} \
            -a {params.gtf_path} \
            -o {output.counts} {input.bam_files}
        module unload subread/2.0.1
        """


# Rule 4a: Clean gene-level count matrix. featureCounts writes a command-line
# comment, six annotation columns and BAM paths as column headers, none of
# which are usable as a delivered count matrix, so split it into plain
# counts/annotation/CPM tables keyed by the metadata sample_id.
rule gene_count_matrix:
    input:
        counts=f"{OUTPUT_DIR}/feature_count/{{species}}_samples_counts.txt",
        metadata=config["deseq2"]["metadata"],
    output:
        counts_csv=f"{OUTPUT_DIR}/counts/{{species}}/gene_counts.csv",
        counts_rds=f"{OUTPUT_DIR}/counts/{{species}}/gene_counts.rds",
        cpm_csv=f"{OUTPUT_DIR}/counts/{{species}}/gene_counts_cpm.csv",
        annotation=f"{OUTPUT_DIR}/counts/{{species}}/gene_annotation.csv",
        metrics=f"{OUTPUT_DIR}/counts/{{species}}/count_matrix_metrics.csv",
    threads: 2
    resources:
        mem_mb=16000,
        cpus=2,
        partition="standard",
        account="sbsandme_lab",
    params:
        out_dir=f"{OUTPUT_DIR}/counts/{{species}}",
    log:
        "logs/gene_count_matrix/{species}.log",
    benchmark:
        "benchmarks/gene_count_matrix/{species}.tsv"
    shell:
        """
        exec > {log} 2>&1
        module load R/4.5.2
        Rscript src/count_matrix.R {input.counts} {input.metadata} {params.out_dir}
        module unload R/4.5.2
        """


# Rule 4c: Sample-level QC -- sequencing/count metrics, sample-sample
# correlation, hierarchical clustering and PCA on variance-stabilised counts.
rule sample_qc:
    input:
        counts_csv=f"{OUTPUT_DIR}/counts/{{species}}/gene_counts.csv",
        metadata=config["deseq2"]["metadata"],
    output:
        metrics=f"{OUTPUT_DIR}/sample_qc/{{species}}/sample_metrics.csv",
        cor_spearman=f"{OUTPUT_DIR}/sample_qc/{{species}}/sample_correlation_spearman.csv",
        cor_pearson=f"{OUTPUT_DIR}/sample_qc/{{species}}/sample_correlation_pearson.csv",
        cor_heatmap=f"{OUTPUT_DIR}/sample_qc/{{species}}/sample_correlation_spearman_heatmap.png",
        distances=f"{OUTPUT_DIR}/sample_qc/{{species}}/sample_distance_euclidean.csv",
        dendrogram=f"{OUTPUT_DIR}/sample_qc/{{species}}/sample_clustering_dendrogram.png",
        pca_csv=f"{OUTPUT_DIR}/sample_qc/{{species}}/pca_coordinates.csv",
        pca_variance=f"{OUTPUT_DIR}/sample_qc/{{species}}/pca_variance_explained.csv",
        pca_plot=f"{OUTPUT_DIR}/sample_qc/{{species}}/pca_plot.png",
        scree_plot=f"{OUTPUT_DIR}/sample_qc/{{species}}/pca_scree_plot.png",
    threads: 2
    resources:
        mem_mb=16000,
        cpus=2,
        partition="standard",
        account="sbsandme_lab",
    params:
        out_dir=f"{OUTPUT_DIR}/sample_qc/{{species}}",
    log:
        "logs/sample_qc/{species}.log",
    benchmark:
        "benchmarks/sample_qc/{species}.tsv"
    shell:
        """
        exec > {log} 2>&1
        module load R/4.5.2
        Rscript src/sample_qc.R {input.counts_csv} {input.metadata} {params.out_dir}
        module unload R/4.5.2
        """


# Rule 4b: rMATS alternative splicing analysis (split by species, same
# reasoning as feature_counts_all -- one GTF per run)
rule rmats:
    input:
        bam_files=lambda wildcards: expand(
            f"{OUTPUT_DIR}/hisat2_alignment/{{sample}}_align_sorted_markdup.bam",
            sample=[
                sample
                for sample in SAMPLES
                if SAMPLE_SPECIES[sample] == wildcards.species
            ],
        ),
    output:
        done=f"{OUTPUT_DIR}/rmats/{{species}}/.done",
    threads: 8
    resources:
        mem_mb=32000,
        cpus=8,
        partition="standard",
        account="sbsandme_lab",
    params:
        bam_list=f"{OUTPUT_DIR}/rmats/{{species}}/bam_files.txt",
        output_dir=f"{OUTPUT_DIR}/rmats/{{species}}",
        gtf_path=lambda wildcards: SPECIES_REFERENCES[wildcards.species]["gtf"],
        read_length=150,
    log:
        "logs/rmats/{species}.log",
    benchmark:
        "benchmarks/rmats/{species}.tsv"
    shell:
        """
        exec > {log} 2>&1
        module load rMATS/4.3.0

        rm -rf {params.output_dir}
        mkdir -p {params.output_dir}/tmp

        # Create BAM file list for rMATS
        echo "{input.bam_files}" | tr ' ' ',' >{params.bam_list}

        # Run rMATS
        rmats.py --b1 {params.bam_list} --gtf {params.gtf_path} \
            -t paired --readLength {params.read_length} \
            --od {params.output_dir} --tmp {params.output_dir}/tmp \
            --nthread {threads}

        touch {output.done}

        module unload rMATS/4.3.0
        """


# Rule 5: Salmon quantification
rule salmon_quant:
    input:
        r1=f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_1P.fq.gz",
        r2=f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_2P.fq.gz",
    output:
        quant=f"{OUTPUT_DIR}/salmon/{{sample}}_salmon_quant/{{sample}}_quant.sf",
    threads: 8
    resources:
        mem_mb=24000,
        cpus=config["params"]["cpus"],
        partition="standard",
        account="sbsandme_lab",
    params:
        salmon_index=lambda wildcards: species_ref(wildcards.sample, "salmon_index"),
        output_dir=f"{OUTPUT_DIR}/salmon/{{sample}}_salmon_quant",
        temp_quant=f"{OUTPUT_DIR}/salmon/{{sample}}_salmon_quant/quant.sf",
    log:
        "logs/salmon_quant/{sample}.log",
    benchmark:
        "benchmarks/salmon_quant/{sample}.tsv"
    shell:
        """
        exec > {log} 2>&1
        module load salmon/1.8.0

        salmon quant -i {params.salmon_index} -l {config[params][salmon][library_type]} \
            -1 {input.r1} -2 {input.r2} \
            -p {threads} --validateMappings --gcBias \
            -o {params.output_dir} \
            --allowDovetail

        # Rename the quant.sf file
        mv {params.temp_quant} {output.quant}

        module unload salmon/1.8.0
        """


# Rule 6: Calculate TPM using tximport from Salmon quantification (per
# species -- tx2gene mapping comes from one GTF, so mouse and human transcript
# IDs can't be imported together)
rule tximport_tpm:
    input:
        quant_files=lambda wildcards: expand(
            f"{OUTPUT_DIR}/salmon/{{sample}}_salmon_quant/{{sample}}_quant.sf",
            sample=[
                sample
                for sample in SAMPLES
                if SAMPLE_SPECIES[sample] == wildcards.species
            ],
        ),
    output:
        tpm_csv=f"{OUTPUT_DIR}/tpm/{{species}}/tpm_salmon.csv",
        tpm_rds=f"{OUTPUT_DIR}/tpm/{{species}}/tpm_salmon.rds",
        txi_rds=f"{OUTPUT_DIR}/tpm/{{species}}/txi_salmon.rds",
        tx_counts_csv=f"{OUTPUT_DIR}/tpm/{{species}}/transcript_counts.csv",
        tx_tpm_csv=f"{OUTPUT_DIR}/tpm/{{species}}/transcript_tpm.csv",
        txi_tx_rds=f"{OUTPUT_DIR}/tpm/{{species}}/txi_transcript_salmon.rds",
    threads: 2
    resources:
        mem_mb=8000,
        cpus=2,
        partition="standard",
        account="sbsandme_lab",
    params:
        staged_salmon_dir=f"{OUTPUT_DIR}/salmon_by_species/{{species}}",
        tpm_dir=f"{OUTPUT_DIR}/tpm/{{species}}",
        gtf_path=lambda wildcards: SPECIES_REFERENCES[wildcards.species]["gtf"],
    log:
        "logs/tximport_tpm/{species}.log",
    benchmark:
        "benchmarks/tximport_tpm/{species}.tsv"
    shell:
        """
        exec > {log} 2>&1
        rm -rf {params.staged_salmon_dir}
        mkdir -p {params.staged_salmon_dir}
        for quant_file in {input.quant_files}; do
            ln -s "$(readlink -f "$(dirname "$quant_file")")" \
                "{params.staged_salmon_dir}/$(basename "$(dirname "$quant_file")")"
        done

        module load R/4.5.2
        Rscript src/tximport_tpm.R {params.staged_salmon_dir} {params.gtf_path} {params.tpm_dir}
        module unload R/4.5.2
        """


# Rule 7: MultiQC report
rule multiqc:
    input:
        expand(f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_1P.fq.gz", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_2P.fq.gz", sample=SAMPLES),
        expand(
            f"{OUTPUT_DIR}/fastqc/{{sample}}/{{sample}}-R1_fastqc.html", sample=SAMPLES
        ),
        expand(
            f"{OUTPUT_DIR}/fastqc/{{sample}}/{{sample}}-R2_fastqc.html", sample=SAMPLES
        ),
        expand(
            f"{OUTPUT_DIR}/hisat2_alignment/{{sample}}_align_sorted_markdup.bam",
            sample=SAMPLES,
        ),
        expand(f"{OUTPUT_DIR}/rustqc/{{sample}}/.done", sample=SAMPLES),
        # expand(f"{OUTPUT_DIR}/salmon/{{sample}}_salmon_quant/{{sample}}_quant.sf", sample=SAMPLES)
    output:
        report=f"{OUTPUT_DIR}/multiqc_report.html",
    threads: 2
    resources:
        mem_mb=4000,
        cpus=2,
        partition="standard",
        account="sbsandme_lab",
    log:
        "logs/multiqc/multiqc.log",
    benchmark:
        "benchmarks/multiqc/multiqc.tsv"
    shell:
        """
        exec > {log} 2>&1
        rm -f {OUTPUT_DIR}/multiqc_report.html {OUTPUT_DIR}/multiqc_report_1.html
        rm -rf {OUTPUT_DIR}/multiqc_data {OUTPUT_DIR}/multiqc_data_1
        module load singularity/3.11.3
        singularity run /dfs9/ucightf-lab/kstachel/TOOLS/multiqc-1.20.sif multiqc {OUTPUT_DIR} -o {OUTPUT_DIR} --force
        module unload singularity/3.11.3
        """


# Rule 8: NCBI GEO/SRA submission package. Checksums for the raw FASTQ files
# are taken from data/FASTQ/md5sums.txt when the sequencing core supplied one,
# so this does not re-hash the full raw data.
rule ncbi_submission:
    input:
        metadata=config["deseq2"]["metadata"],
        counts_csv=f"{OUTPUT_DIR}/counts/{{species}}/gene_counts.csv",
        tpm_csv=f"{OUTPUT_DIR}/tpm/{{species}}/tpm_salmon.csv",
    output:
        geo=f"{OUTPUT_DIR}/ncbi_submission/{{species}}/geo_samples.csv",
        sra=f"{OUTPUT_DIR}/ncbi_submission/{{species}}/sra_metadata.csv",
        md5=f"{OUTPUT_DIR}/ncbi_submission/{{species}}/md5sums.txt",
        readme=f"{OUTPUT_DIR}/ncbi_submission/{{species}}/SUBMISSION_README.txt",
    threads: 1
    resources:
        mem_mb=4000,
        cpus=1,
        partition="standard",
        account="sbsandme_lab",
    params:
        out_dir=f"{OUTPUT_DIR}/ncbi_submission/{{species}}",
        counts_dir=f"{OUTPUT_DIR}/counts/{{species}}",
        tpm_dir=f"{OUTPUT_DIR}/tpm/{{species}}",
        tissue=config["ncbi"]["tissue"],
        instrument_model=config["ncbi"]["instrument_model"],
        genome_build=lambda wildcards: os.path.basename(
            SPECIES_REFERENCES[wildcards.species]["hisat2_index"]
        ),
        annotation=lambda wildcards: os.path.basename(
            SPECIES_REFERENCES[wildcards.species]["gtf"]
        ),
    log:
        "logs/ncbi_submission/{species}.log",
    benchmark:
        "benchmarks/ncbi_submission/{species}.tsv"
    shell:
        """
        exec > {log} 2>&1
        python3 src/prepare_ncbi_submission.py \
            --metadata {input.metadata} \
            --fastq-dir {DATA_PATH} \
            --processed-dir {params.counts_dir} \
            --tpm-dir {params.tpm_dir} \
            --species {wildcards.species} \
            --tissue "{params.tissue}" \
            --instrument-model "{params.instrument_model}" \
            --genome-build "{params.genome_build}" \
            --annotation "{params.annotation}" \
            --output-dir {params.out_dir}
        """


# Rule 7: Generate project report
#
# limma-voom is an input because it is this project's differential-expression
# method and belongs in every report. ncbi_submission and deseq2 deliberately
# are not: generate_report.py reads whatever output of theirs happens to be on
# disk and prints a "not found, run the rule" note when there is none, so
# listing them here would only force those rules to run in every build. Run
# them on demand instead, e.g.
#   snakemake output/ncbi_submission/mouse/geo_samples.csv
# and re-run this rule afterwards to fold their results into the report.
rule generate_report:
    input:
        counts=expand(
            f"{OUTPUT_DIR}/counts/{{species}}/gene_counts.csv", species=SPECIES_LIST
        ),
        limma_voom=expand(
            f"{OUTPUT_DIR}/limma_voom/{{species}}/limma_voom_comparisons_manifest.csv",
            species=SPECIES_LIST,
        ),
        pca=expand(
            f"{OUTPUT_DIR}/limma_voom/{{species}}/pca_all_samples_{{colour_var}}.png",
            species=SPECIES_LIST,
            colour_var=["condition", "age_group"],
        ),
        multiqc=f"{OUTPUT_DIR}/multiqc_report.html",
        metadata=config["deseq2"]["metadata"],
    output:
        report="RNAseq_Project_Report.pdf",
    threads: 1
    resources:
        mem_mb=4000,
        cpus=1,
        partition="standard",
        account="sbsandme_lab",
    params:
        species=",".join(SPECIES_LIST),
    log:
        "logs/generate_report/generate_report.log",
    benchmark:
        "benchmarks/generate_report/generate_report.tsv"
    shell:
        """
        exec > {log} 2>&1
        python3 src/generate_report.py \
            --fastq-dir {DATA_PATH} \
            --metadata {input.metadata} \
            --species {params.species} \
            --output {output.report}
        """


# Rule 9: DESeq2 differential expression analysis (per species -- gene IDs
# and comparisons don't carry across organisms, so this is not a combined run)
rule deseq2:
    input:
        counts=f"{OUTPUT_DIR}/feature_count/{{species}}_samples_counts.txt",
        metadata=config["deseq2"]["metadata"],
        comparisons_config=config["deseq2"]["comparisons_config"],
        script="src/deseq2_analysis.R",
    output:
        results=f"{OUTPUT_DIR}/deseq2/{{species}}/deseq2_results.csv",
        rds=f"{OUTPUT_DIR}/deseq2/{{species}}/dds.rds",
        manifest=f"{OUTPUT_DIR}/deseq2/{{species}}/deseq2_comparisons_manifest.csv",
    threads: 1
    resources:
        mem_mb=8000,
        cpus=1,
        partition="standard",
        account="sbsandme_lab",
    params:
        out_dir=f"{OUTPUT_DIR}/deseq2/{{species}}",
    log:
        "logs/deseq2/{species}.log",
    benchmark:
        "benchmarks/deseq2/{species}.tsv"
    shell:
        """
        exec > {log} 2>&1
        module load R/4.5.2
        Rscript {input.script} {input.counts} {input.metadata} \
            {params.out_dir} {input.comparisons_config}
        module unload R/4.5.2
        """



# Rule 10: limma-voom differential expression analysis (per species). This is
# the project's DE method, so unlike deseq2 it IS in `rule all`. Runs the same
# comparisons config as deseq2 so the two methods can be compared
# contrast-for-contrast.
#
# One model is fitted per species over all that species' samples, on the
# age_group x region cells, blocked on the animal -- every mouse contributed all
# three brain regions, so the samples are not independent. Each comparison in
# the config is a contrast of that single fit. `block_var` and `quality_weights`
# are the two modelling knobs; see their comments in config.yaml.
# Rule 9b: sample-tracking gate. Infers each library's sex from Y-linked genes
# and Xist and refuses to let the DE rule run on material that contradicts
# itself -- a library expressing both marker sets is contaminated or pooled, and
# an animal whose brain regions disagree about which set they express has a
# labelling problem no model can absorb. The `sex` column in metadata.csv is
# empty for this project, so expression is the only sample-tracking evidence
# there is.
#
# This runs BEFORE limma_voom rather than inside it because a swap invalidates
# the blocking that the whole pooled design rests on; discovering it in a
# warning halfway down the DE log is too late. Set sex_qc.strict to false in
# config.yaml to downgrade it to a report and let the pipeline through.
#
# On failure Snakemake removes the output CSV, so the script also prints the
# full per-library table into logs/sex_qc/{species}.log -- read that.
rule sex_qc:
    input:
        counts=f"{OUTPUT_DIR}/feature_count/{{species}}_samples_counts.txt",
        metadata=config["deseq2"]["metadata"],
        script="src/infer_sex.R",
    output:
        inferred_sex=f"{OUTPUT_DIR}/sex_qc/{{species}}/inferred_sex.csv",
    threads: 1
    resources:
        mem_mb=8000,
        cpus=1,
        partition="standard",
        account="sbsandme_lab",
    params:
        out_dir=f"{OUTPUT_DIR}/sex_qc/{{species}}",
        strict=lambda w: "--strict" if config.get("sex_qc", {}).get("strict", True) else "",
    log:
        "logs/sex_qc/{species}.log",
    benchmark:
        "benchmarks/sex_qc/{species}.tsv"
    shell:
        """
        exec > {log} 2>&1
        module load R/4.5.2
        Rscript {input.script} {input.counts} {input.metadata} \
            {params.out_dir} {params.strict}
        module unload R/4.5.2
        """


rule limma_voom:
    input:
        counts=f"{OUTPUT_DIR}/feature_count/{{species}}_samples_counts.txt",
        metadata=config["deseq2"]["metadata"],
        comparisons_config=config["limma_voom"]["comparisons_config"],
        script="src/limma_voom_analysis.R",
        sex_module="src/infer_sex.R",
        sex_qc=f"{OUTPUT_DIR}/sex_qc/{{species}}/inferred_sex.csv",
    output:
        results=f"{OUTPUT_DIR}/limma_voom/{{species}}/limma_voom_results.csv",
        rds=f"{OUTPUT_DIR}/limma_voom/{{species}}/efit.rds",
        manifest=f"{OUTPUT_DIR}/limma_voom/{{species}}/limma_voom_comparisons_manifest.csv",
        pca_condition=f"{OUTPUT_DIR}/limma_voom/{{species}}/pca_all_samples_condition.png",
        pca_age_group=f"{OUTPUT_DIR}/limma_voom/{{species}}/pca_all_samples_age_group.png",
        pca_inferred_sex=f"{OUTPUT_DIR}/limma_voom/{{species}}/pca_all_samples_inferred_sex.png",
        inferred_sex=f"{OUTPUT_DIR}/limma_voom/{{species}}/inferred_sex.csv",
    threads: 1
    resources:
        mem_mb=8000,
        cpus=1,
        partition="standard",
        account="sbsandme_lab",
    params:
        out_dir=f"{OUTPUT_DIR}/limma_voom/{{species}}",
        block_var=config["limma_voom"].get("block_var", "mouse_id"),
        quality_weights=config["limma_voom"].get("quality_weights", "per_sample"),
    log:
        "logs/limma_voom/{species}.log",
    benchmark:
        "benchmarks/limma_voom/{species}.tsv"
    shell:
        """
        exec > {log} 2>&1
        module load R/4.5.2
        Rscript {input.script} {input.counts} {input.metadata} \
            {params.out_dir} {input.comparisons_config} {params.block_var} \
            {params.quality_weights}
        module unload R/4.5.2
        """
