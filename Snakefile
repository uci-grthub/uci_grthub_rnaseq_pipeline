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
RUSTQC_CONTAINER = "/dfs9/ucightf-lab/kstachel/containers/rustqc.sif"

# Per-sample species -> reference paths, so samples from different organisms
# (e.g. mouse libraries mixed into an otherwise-human run) get aligned and
# quantified against the correct genome (config.species_references.yaml).
# There is deliberately no built-in default: a macaque library aligned to
# GRCh38 still clears the alignment-rate gate, so a blank species must stop
# the run rather than silently pick a genome. A project can opt in to one by
# setting default_species in config.local.yaml.
DEFAULT_SPECIES = config.get("default_species") or ""
SPECIES_REFERENCES = config["species_references"]

# Keyed by FASTQ sample name. The metadata `sample` column may hold that name
# or an experimental ID (hFB-1); in the latter case the i7/i5 pair ending the
# FASTQ name is the only link, so it is the fallback match (as in count_matrix.R).
def _barcode_key(i7, i5):
    return f"{i7.strip().upper()}-{i5.strip().upper()}" if i7 and i5 else ""


with open(METADATA_PATH, newline="") as fh:
    _metadata_rows = list(csv.DictReader(fh))

_meta_by_name = {row["sample"]: row for row in _metadata_rows}
_meta_by_barcode = {}
for _row in _metadata_rows:
    _i7 = next((v for k, v in _row.items() if k and k.lower().startswith("i7")), "")
    _i5 = next((v for k, v in _row.items() if k and k.lower().startswith("i5")), "")
    _key = _barcode_key(_i7 or "", _i5 or "")
    if _key:
        _meta_by_barcode.setdefault(_key, []).append(_row)

SAMPLE_SPECIES = {}
_unmatched = []
_no_species = []
for sample in SAMPLES:
    row = _meta_by_name.get(sample)
    if row is None:
        candidates = _meta_by_barcode.get(
            _barcode_key(*_infer_barcodes_from_sample(sample)), []
        )
        row = candidates[0] if len(candidates) == 1 else None
    if row is None:
        _unmatched.append(sample)
    species = ((row or {}).get("species") or "").strip() or DEFAULT_SPECIES
    if not species:
        _no_species.append(sample)
    SAMPLE_SPECIES[sample] = species

if _no_species:
    raise ValueError(
        f"No species for {len(_no_species)} sample(s): {_no_species}. Fill the "
        f"`species` column of {METADATA_PATH} (one of "
        f"{sorted(SPECIES_REFERENCES)}), or set default_species in config.local.yaml."
    )

_unknown = sorted(
    {species for species in SAMPLE_SPECIES.values()} - set(SPECIES_REFERENCES)
)
if _unknown:
    raise ValueError(
        f"Unknown species {_unknown} in {METADATA_PATH} or default_species; "
        f"config.species_references.yaml defines {sorted(SPECIES_REFERENCES)}."
    )

if _unmatched:
    print(
        f"WARNING: no metadata row for {len(_unmatched)} sample(s), using "
        f"default_species '{DEFAULT_SPECIES}': {_unmatched}",
        file=sys.stderr,
    )

SPECIES_LIST = sorted(
    {SAMPLE_SPECIES.get(sample, DEFAULT_SPECIES) for sample in SAMPLES}
)


# A library sequenced on several lanes is trimmed per lane and aligned once,
# with HISAT2 given every lane's reads, so duplicate marking, counting and QC
# see one BAM per library. The metadata `library` column names the library; a
# sample without one is its own library.
def _library_of(sample):
    row = _meta_by_name.get(sample)
    if row is None:
        candidates = _meta_by_barcode.get(
            _barcode_key(*_infer_barcodes_from_sample(sample)), []
        )
        row = candidates[0] if len(candidates) == 1 else None
    return ((row or {}).get("library") or "").strip() or sample


LIBRARY_LANES = {}
for sample in SAMPLES:
    LIBRARY_LANES.setdefault(_library_of(sample), []).append(sample)
LIBRARIES = sorted(LIBRARY_LANES)

LIBRARY_SPECIES = {}
for _library, _lanes in LIBRARY_LANES.items():
    _species = {SAMPLE_SPECIES[lane] for lane in _lanes}
    if len(_species) > 1:
        raise ValueError(
            f"Lanes of library {_library} have different species: {sorted(_species)}"
        )
    LIBRARY_SPECIES[_library] = _species.pop()

print(
    f"Aligning {len(SAMPLES)} lane-level samples as {len(LIBRARIES)} libraries",
    file=sys.stderr,
)


def library_ref(library, key):
    return SPECIES_REFERENCES[LIBRARY_SPECIES[library]][key]


def libraries_of_species(species):
    return [lib for lib in LIBRARIES if LIBRARY_SPECIES[lib] == species]


def lane_files(pattern):
    """Input function: `pattern` expanded over the lanes of the library."""
    return lambda wildcards: expand(pattern, sample=LIBRARY_LANES[wildcards.library])


wildcard_constraints:
    sample="|".join(re.escape(s) for s in SAMPLES) or "[^/]+",
    library="|".join(re.escape(lib) for lib in LIBRARIES) or "[^/]+",


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
        expand(f"{OUTPUT_DIR}/rustqc/{{library}}/.done", library=LIBRARIES),
        # Trimmed files
        expand(f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_1P.fq.gz", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_2P.fq.gz", sample=SAMPLES),
        # HISAT2 alignment (one BAM per library) and counting
        expand(
            f"{OUTPUT_DIR}/hisat2_alignment/{{library}}_align_sorted_markdup.bam",
            library=LIBRARIES,
        ),
        expand(
            f"{OUTPUT_DIR}/hisat2_alignment/{{library}}_align_sorted_markdup.bam.bai",
            library=LIBRARIES,
        ),
        expand(
            f"{OUTPUT_DIR}/feature_count/{{species}}_samples_counts.txt",
            species=SPECIES_LIST,
        ),
        # Clean gene-level raw count matrix (deliverable form of featureCounts)
        expand(f"{OUTPUT_DIR}/counts/{{species}}/gene_counts.csv", species=SPECIES_LIST),
        # Sample correlation, clustering and PCA
        expand(f"{OUTPUT_DIR}/sample_qc/{{species}}/pca_plot.png", species=SPECIES_LIST),
        # Salmon quantification (one per library)
        expand(
            f"{OUTPUT_DIR}/salmon/{{library}}_salmon_quant/{{library}}_quant.sf",
            library=LIBRARIES,
        ),
        # TPM quantification using tximport
        expand(f"{OUTPUT_DIR}/tpm/{{species}}/tpm_salmon.csv", species=SPECIES_LIST),
        # Transcript-level count and TPM matrices (tximport txOut=TRUE)
        expand(
            f"{OUTPUT_DIR}/tpm/{{species}}/transcript_counts.csv",
            species=SPECIES_LIST,
        ),
        # MultiQC report
        f"{OUTPUT_DIR}/multiqc_report.html",
        # # GEO/SRA submission sheets and checksums
        # expand(
        #     f"{OUTPUT_DIR}/ncbi_submission/{{species}}/geo_samples.csv",
        #     species=SPECIES_LIST,
        # ),
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
        bam=f"{OUTPUT_DIR}/hisat2_alignment/{{library}}_align_sorted_markdup.bam",
    output:
        done=f"{OUTPUT_DIR}/rustqc/{{library}}/.done",
    threads: 8
    resources:
        mem_mb=24000,
        cpus=config["params"]["cpus"],
        partition="standard",
        account="sbsandme_lab",
    params:
        out_dir=f"{OUTPUT_DIR}/rustqc/{{library}}",
        gtf_path=lambda wildcards: library_ref(wildcards.library, "gtf"),
    log:
        "logs/rustqc/{library}.log",
    benchmark:
        "benchmarks/rustqc/{library}.tsv"
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


# Rule 1: Adapter and quality trimming with chelae. Pairs are kept or dropped
# together, so there are no unpaired outputs.
rule chelae_trim:
    input:
        r1=f"{DATA_PATH}/{{sample}}-R1.fastq.gz",
        r2=f"{DATA_PATH}/{{sample}}-R2.fastq.gz",
    output:
        r1_paired=f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_1P.fq.gz",
        r2_paired=f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_2P.fq.gz",
        # fastp-shape JSON; the name matches MultiQC's fastp module pattern
        json=f"{OUTPUT_DIR}/trimmed/{{sample}}.fastp.json",
        metrics=f"{OUTPUT_DIR}/trimmed/{{sample}}_chelae_metrics.tsv",
    threads: 8
    resources:
        mem_mb=4000,
        cpus=config["params"]["cpus"],
        partition="standard",
        account="sbsandme_lab",
    log:
        f"{OUTPUT_DIR}/trimmed/{{sample}}_chelae.log",
    benchmark:
        "benchmarks/chelae_trim/{sample}.tsv"
    shell:
        """
        exec > {log} 2>&1
        chelae trim \
            -i {input.r1} {input.r2} \
            -o {output.r1_paired} {output.r2_paired} \
            -t {threads} \
            --kit {config[params][chelae][kit]} \
            --quality-trim-5p {config[params][chelae][quality_trim]} \
            -l {config[params][chelae][min_length]} \
            -j {output.json} \
            -m {output.metrics}
        """


# Rule 2: HISAT2 alignment, one run per library over all of its lanes (HISAT2
# takes comma-separated lists for -1/-2, kept in the same lane order)
rule hisat2_align:
    input:
        r1=lane_files(f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_1P.fq.gz"),
        r2=lane_files(f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_2P.fq.gz"),
    output:
        bam=f"{OUTPUT_DIR}/hisat2_alignment/{{library}}_align.bam",
        summary=f"{OUTPUT_DIR}/hisat2_alignment/alignment_summary/{{library}}_summary.align",
    threads: 8
    resources:
        mem_mb=24000,
        cpus=config["params"]["cpus"],
        partition="standard",
        account="sbsandme_lab",
    params:
        hisat2_index=lambda wildcards: library_ref(wildcards.library, "hisat2_index"),
        summary_path=f"{OUTPUT_DIR}/hisat2_alignment/alignment_summary",
        r1_list=lambda wildcards, input: ",".join(input.r1),
        r2_list=lambda wildcards, input: ",".join(input.r2),
    log:
        "logs/hisat2_align/{library}.log",
    benchmark:
        "benchmarks/hisat2_align/{library}.tsv"
    shell:
        """
        exec 2> {log}
        module load hisat2/2.2.1
        module load samtools/1.15.1

        hisat2 -p {threads} -t --qc-filter --rna-strandness {config[params][hisat2][rna_strandness]} \
            --summary-file {output.summary} \
            -x {params.hisat2_index} --dta-cufflinks \
            --rg-id {wildcards.library} --rg SM:{wildcards.library} \
            --rg LB:{wildcards.library} --rg PL:ILLUMINA \
            -1 {params.r1_list} -2 {params.r2_list} \
            | samtools sort -n -@ 2 \
            | samtools fixmate -m -@ 2 - {output.bam}

        module unload samtools/1.15.1
        module unload hisat2/2.2.1
        """


# Rule 3: Sort and index BAM file
rule sort_bam:
    input:
        bam=f"{OUTPUT_DIR}/hisat2_alignment/{{library}}_align.bam",
    output:
        sorted_bam=f"{OUTPUT_DIR}/hisat2_alignment/{{library}}_align_sorted.bam",
        index=f"{OUTPUT_DIR}/hisat2_alignment/{{library}}_align_sorted.bam.bai",
    threads: 8
    resources:
        mem_mb=24000,
        cpus=config["params"]["cpus"],
        partition="standard",
        account="sbsandme_lab",
    log:
        "logs/sort_bam/{library}.log",
    benchmark:
        "benchmarks/sort_bam/{library}.tsv"
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
        sorted_bam=f"{OUTPUT_DIR}/hisat2_alignment/{{library}}_align_sorted.bam",
    output:
        markdup_bam=f"{OUTPUT_DIR}/hisat2_alignment/{{library}}_align_sorted_markdup.bam",
        markdup_bai=f"{OUTPUT_DIR}/hisat2_alignment/{{library}}_align_sorted_markdup.bam.bai",
        metrics=f"{OUTPUT_DIR}/hisat2_alignment/{{library}}_markdup_metrics.txt",
    threads: 8
    resources:
        mem_mb=24000,
        cpus=config["params"]["cpus"],
        partition="standard",
        account="sbsandme_lab",
    log:
        "logs/markdup/{library}.log",
    benchmark:
        "benchmarks/markdup/{library}.tsv"
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
            f"{OUTPUT_DIR}/hisat2_alignment/{{library}}_align_sorted_markdup.bam",
            library=libraries_of_species(wildcards.species),
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
        ensdb=lambda wildcards: SPECIES_REFERENCES[wildcards.species].get("ensdb", ""),
    log:
        "logs/gene_count_matrix/{species}.log",
    benchmark:
        "benchmarks/gene_count_matrix/{species}.tsv"
    shell:
        """
        exec > {log} 2>&1
        module load R/4.5.2
        Rscript src/count_matrix.R {input.counts} {input.metadata} {params.out_dir} "{params.ensdb}"
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
            f"{OUTPUT_DIR}/hisat2_alignment/{{library}}_align_sorted_markdup.bam",
            library=libraries_of_species(wildcards.species),
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


# Rule 5: Salmon quantification, one run per library over all of its lanes
# (salmon takes space-separated lists for -1/-2, kept in the same lane order),
# so transcript matrices have the same library columns as the gene counts
rule salmon_quant:
    input:
        r1=lane_files(f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_1P.fq.gz"),
        r2=lane_files(f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_2P.fq.gz"),
    output:
        quant=f"{OUTPUT_DIR}/salmon/{{library}}_salmon_quant/{{library}}_quant.sf",
    threads: 8
    resources:
        # decoy-aware index is ~17 GB on disk; measured peak RSS ~18.5 GB
        mem_mb=24000,
        cpus=config["params"]["cpus"],
        partition="standard",
        account="sbsandme_lab",
    params:
        salmon_index=lambda wildcards: library_ref(wildcards.library, "salmon_index"),
        output_dir=f"{OUTPUT_DIR}/salmon/{{library}}_salmon_quant",
        temp_quant=f"{OUTPUT_DIR}/salmon/{{library}}_salmon_quant/quant.sf",
    log:
        "logs/salmon_quant/{library}.log",
    benchmark:
        "benchmarks/salmon_quant/{library}.tsv"
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
            f"{OUTPUT_DIR}/salmon/{{library}}_salmon_quant/{{library}}_quant.sf",
            library=libraries_of_species(wildcards.species),
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
        ensdb=lambda wildcards: SPECIES_REFERENCES[wildcards.species].get("ensdb", ""),
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
        Rscript src/tximport_tpm.R {params.staged_salmon_dir} {params.gtf_path} {params.tpm_dir} "{params.ensdb}"
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
            f"{OUTPUT_DIR}/hisat2_alignment/{{library}}_align_sorted_markdup.bam",
            library=LIBRARIES,
        ),
        expand(f"{OUTPUT_DIR}/rustqc/{{library}}/.done", library=LIBRARIES),
        expand(
            f"{OUTPUT_DIR}/salmon/{{library}}_salmon_quant/{{library}}_quant.sf",
            library=LIBRARIES,
        ),
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
rule generate_report:
    input:
        counts=expand(
            f"{OUTPUT_DIR}/counts/{{species}}/gene_counts.csv", species=SPECIES_LIST
        ),
        transcript_counts=expand(
            f"{OUTPUT_DIR}/tpm/{{species}}/transcript_counts.csv", species=SPECIES_LIST
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
    output:
        results=f"{OUTPUT_DIR}/deseq2/{{species}}/deseq2_results.csv",
        rds=f"{OUTPUT_DIR}/deseq2/{{species}}/dds.rds",
        manifest=f"{OUTPUT_DIR}/deseq2/{{species}}/deseq2_comparisons_manifest.csv",
        targets=f"{OUTPUT_DIR}/deseq2/{{species}}/target_gene_results.csv",
        per_region=f"{OUTPUT_DIR}/deseq2/{{species}}/target_per_region_log2fc.csv",
    threads: 1
    resources:
        mem_mb=8000,
        cpus=1,
        partition="standard",
        account="sbsandme_lab",
    params:
        out_dir=f"{OUTPUT_DIR}/deseq2/{{species}}",
        ensdb=lambda wildcards: SPECIES_REFERENCES[wildcards.species].get("ensdb", ""),
        gtf=lambda wildcards: SPECIES_REFERENCES[wildcards.species]["gtf"],
    log:
        "logs/deseq2/{species}.log",
    benchmark:
        "benchmarks/deseq2/{species}.tsv"
    shell:
        """
        exec > {log} 2>&1
        module load R/4.5.2
        Rscript src/deseq2_analysis.R {input.counts} {input.metadata} \
            {params.out_dir} {input.comparisons_config} "{params.ensdb}" "{params.gtf}"
        module unload R/4.5.2
        """

