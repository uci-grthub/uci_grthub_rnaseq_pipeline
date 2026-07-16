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
DEFAULT_SPECIES = config.get("default_species", "human")
SPECIES_REFERENCES = config["species_references"]

SAMPLE_SPECIES = {}
with open(METADATA_PATH, newline="") as fh:
    for row in csv.DictReader(fh):
        species = (row.get("species") or "").strip()
        SAMPLE_SPECIES[row["sample"]] = species if species else DEFAULT_SPECIES

SPECIES_LIST = sorted(
    {SAMPLE_SPECIES.get(sample, DEFAULT_SPECIES) for sample in SAMPLES}
)


def species_ref(sample, key):
    species = SAMPLE_SPECIES.get(sample, DEFAULT_SPECIES)
    return SPECIES_REFERENCES[species][key]


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
        expand(f"{OUTPUT_DIR}/rmats/{{species}}/.done", species=SPECIES_LIST),
        # Salmon quantification
        expand(
            f"{OUTPUT_DIR}/salmon/{{sample}}_salmon_quant/{{sample}}_quant.sf",
            sample=SAMPLES,
        ),
        # TPM quantification using tximport
        expand(f"{OUTPUT_DIR}/tpm/{{species}}/tpm_salmon.csv", species=SPECIES_LIST),
        # MultiQC report
        f"{OUTPUT_DIR}/multiqc_report.html",
        # Project report
        "RNAseq_Project_Report.pdf",
        # DESeq2 results not included by default -- run on demand with e.g.
        # `snakemake output/deseq2/human/deseq2_results.csv`
        # # iSEE app2.R file
        # "isee_uci/shiny-server/test_app/app.R"


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
                if SAMPLE_SPECIES.get(sample, DEFAULT_SPECIES) == wildcards.species
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


# Rule 4b: rMATS alternative splicing analysis (split by species, same
# reasoning as feature_counts_all -- one GTF per run)
rule rmats:
    input:
        bam_files=lambda wildcards: expand(
            f"{OUTPUT_DIR}/hisat2_alignment/{{sample}}_align_sorted_markdup.bam",
            sample=[
                sample
                for sample in SAMPLES
                if SAMPLE_SPECIES.get(sample, DEFAULT_SPECIES) == wildcards.species
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
                if SAMPLE_SPECIES.get(sample, DEFAULT_SPECIES) == wildcards.species
            ],
        ),
    output:
        tpm_csv=f"{OUTPUT_DIR}/tpm/{{species}}/tpm_salmon.csv",
        tpm_rds=f"{OUTPUT_DIR}/tpm/{{species}}/tpm_salmon.rds",
        txi_rds=f"{OUTPUT_DIR}/tpm/{{species}}/txi_salmon.rds",
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

        module load R/4.2.2
        Rscript src/tximport_tpm.R {params.staged_salmon_dir} {params.gtf_path} {params.tpm_dir}
        module unload R/4.2.2
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


# Rule 7: Generate project report
rule generate_report:
    input:
        counts=expand(
            f"{OUTPUT_DIR}/feature_count/{{species}}_samples_counts.txt",
            species=SPECIES_LIST,
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
    log:
        "logs/generate_report/generate_report.log",
    benchmark:
        "benchmarks/generate_report/generate_report.tsv"
    shell:
        """
        exec > {log} 2>&1
        python3 proj_src/generate_report.py \
            --fastq-dir {DATA_PATH} \
            --metadata {input.metadata} \
            --output {output.report}
        """


# Rule 9: DESeq2 differential expression analysis (per species -- gene IDs
# and comparisons don't carry across organisms, so this is not a combined run)
rule deseq2:
    input:
        counts=f"{OUTPUT_DIR}/feature_count/{{species}}_samples_counts.txt",
        metadata=config["deseq2"]["metadata"],
        comparisons_config=config["deseq2"]["comparisons_config"],
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
        Rscript proj_src/deseq2_analysis.R {input.counts} {input.metadata} \
            {params.out_dir} {input.comparisons_config}
        module unload R/4.5.2
        if [ "{wildcards.species}" = "{DEFAULT_SPECIES}" ]; then
            cp {output.rds} isee_uci/shiny-server/test_app/dds.rds
        fi
        """


# Rule 10: Generate parametric iSEE app2.R file (built from the default
# species' DESeq2 run -- the iSEE app/comparisons config is human-specific)
rule generate_isee_app:
    input:
        template="templates/app.R.template",
        dds=f"{OUTPUT_DIR}/deseq2/{DEFAULT_SPECIES}/dds.rds",
    output:
        app="isee_uci/shiny-server/test_app/app.R",
    params:
        deseq2_condition=config["isee_app"]["condition"],
        group_a=config["isee_app"]["group_a"],
        group_b=config["isee_app"]["group_b"],
    log:
        "logs/generate_isee_app/generate_isee_app.log",
    benchmark:
        "benchmarks/generate_isee_app/generate_isee_app.tsv"
    shell:
        """
        exec > {log} 2>&1
        # Create the output directory if it doesn't exist
        mkdir -p $(dirname {output.app})

        # Replace placeholders in template with actual parameters
        sed -e 's|{{DESEQ2_CONDITION}}|{params.deseq2_condition}|g' \
            -e 's|{{GROUP_A}}|{params.group_a}|g' \
            -e 's|{{GROUP_B}}|{params.group_b}|g' \
            {input.template} >{output.app}
        """
