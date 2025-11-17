EMAIL = "kstachel@uci.edu"

onstart:
   shell("mail -s 'STARTED' {EMAIL} < {log}")

onsuccess:
   shell("mail -s 'DONE' {EMAIL} < {log}")

onerror:
   shell("mail -s 'ERROR' {EMAIL} < {log}")

# Snakemake workflow for RNA-seq analysis
# Generalized to process multiple samples from a data directory

import glob
import os

# Load configuration
configfile: "config.yaml"

# Extract configuration variables
DATA_PATH = config["paths"]["data"]
TRIMMED_PATH = config["paths"]["trimmed"]
HISAT2_PATH = config["paths"]["hisat2"]
ALIGNMENT_SUMMARY_PATH = config["paths"]["alignment_summary"]
FEATURE_COUNT_PATH = config["paths"]["feature_count"]
DESEQ2_PATH = config["paths"]["deseq2"]
SALMON_PATH = config["paths"]["salmon"]

# Get sample names from FASTQ files in data directory
FASTQ_FILES = glob.glob(f"{DATA_PATH}/*_r1.fq.gz")
SAMPLES = [os.path.basename(f).replace("_r1.fq.gz", "") for f in FASTQ_FILES]

print(f"Found {len(SAMPLES)} samples: {SAMPLES}")

# Reference paths
ADAPTER_PATH = config["references"]["adapters"]
HISAT2_INDEX = config["references"]["hisat2_index"]
GTF_PATH = config["references"]["gtf"]
SALMON_INDEX = config["references"]["salmon_index"]
TRIMMOMATIC_JAR = config["tools"]["trimmomatic"]

# Rule all - defines final outputs
rule all:
    input:
        # FastQC reports
        expand(f"fastqc/{{sample}}_r1_fastqc.html", sample=SAMPLES),
        expand(f"fastqc/{{sample}}_r2_fastqc.html", sample=SAMPLES),
        # Trimmed files
        expand(f"{TRIMMED_PATH}/{{sample}}_trimmed_1P.fq.gz", sample=SAMPLES),
        expand(f"{TRIMMED_PATH}/{{sample}}_trimmed_2P.fq.gz", sample=SAMPLES),
        # HISAT2 alignment and counting
        expand(f"{HISAT2_PATH}/{{sample}}_align_sorted.bam", sample=SAMPLES),
        expand(f"{HISAT2_PATH}/{{sample}}_align_sorted.bam.bai", sample=SAMPLES),
        f"{FEATURE_COUNT_PATH}/all_samples_counts.txt",
        # Salmon quantification
        expand(f"{SALMON_PATH}/{{sample}}_salmon_quant/{{sample}}_quant.sf", sample=SAMPLES),
        # MultiQC report
        "multiqc_report.html",
        # DESeq2 results
        f"{DESEQ2_PATH}/deseq2_results.csv",
        # iSEE app2.R file
        "isee_uci/shiny-server/test_app/app.R"

# Rule 0: FastQC on raw FASTQ files
rule fastqc:
    input:
        r1 = f"{DATA_PATH}/{{sample}}_r1.fq.gz",
        r2 = f"{DATA_PATH}/{{sample}}_r2.fq.gz"
    output:
        r1_html = f"fastqc/{{sample}}_r1_fastqc.html",
        r1_zip = f"fastqc/{{sample}}_r1_fastqc.zip",
        r2_html = f"fastqc/{{sample}}_r2_fastqc.html",
        r2_zip = f"fastqc/{{sample}}_r2_fastqc.zip"
    threads: 2
    resources:
        mem_mb = 4000,
        cpus = 2,
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        module load fastqc/0.11.9
        mkdir -p fastqc
        fastqc -o fastqc -t {threads} {input.r1} {input.r2}
        module unload fastqc/0.11.9
        """

# Rule 1: Trimming with Trimmomatic
rule trimmomatic:
    input:
        r1 = f"{DATA_PATH}/{{sample}}_r1.fq.gz",
        r2 = f"{DATA_PATH}/{{sample}}_r2.fq.gz"
    output:
        r1_paired = f"{TRIMMED_PATH}/{{sample}}_trimmed_1P.fq.gz",
        r1_unpaired = f"{TRIMMED_PATH}/{{sample}}_trimmed_1U.fq.gz",
        r2_paired = f"{TRIMMED_PATH}/{{sample}}_trimmed_2P.fq.gz",
        r2_unpaired = f"{TRIMMED_PATH}/{{sample}}_trimmed_2U.fq.gz"
    params:
        adapter_path = ADAPTER_PATH,
        trimmed_base = f"{TRIMMED_PATH}/{{sample}}_trimmed.fq.gz"
    threads: 8
    resources:
        mem_mb = 4000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """

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
        r1 = f"{TRIMMED_PATH}/{{sample}}_trimmed_1P.fq.gz",
        r2 = f"{TRIMMED_PATH}/{{sample}}_trimmed_2P.fq.gz"
    output:
        bam = f"{HISAT2_PATH}/{{sample}}_align.bam",
        summary = f"{ALIGNMENT_SUMMARY_PATH}/{{sample}}_summary.align"
    params:
        hisat2_index = HISAT2_INDEX,
        summary_path = ALIGNMENT_SUMMARY_PATH
    threads: 8
    resources:
        mem_mb = 24000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        module load hisat2/2.2.1
        module load samtools/1.10
        
        hisat2 -p {threads} -t --qc-filter --rna-strandness {config[params][hisat2][rna_strandness]} \
        --summary-file {output.summary} \
        -x {params.hisat2_index} --dta-cufflinks \
        -1 {input.r1} -2 {input.r2} | \
        samtools view -@ {threads} -bS > {output.bam}
        
        module unload samtools/1.10
        module unload hisat2/2.2.1
        """

# Rule 3: Sort and index BAM file
rule sort_bam:
    input:
        bam = f"{HISAT2_PATH}/{{sample}}_align.bam"
    output:
        sorted_bam = f"{HISAT2_PATH}/{{sample}}_align_sorted.bam",
        index = f"{HISAT2_PATH}/{{sample}}_align_sorted.bam.bai"
    threads: 8
    resources:
        mem_mb = 24000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        module load samtools/1.10
        
        samtools sort -@ {threads} -o {output.sorted_bam} {input.bam}
        samtools index -@ {threads} {output.sorted_bam}
        
        module unload samtools/1.10
        """


# Rule 4: Feature counting (all samples together)
rule feature_counts_all:
    input:
        bam_files = expand(f"{HISAT2_PATH}/{{sample}}_align_sorted.bam", sample=SAMPLES)
    output:
        counts = f"{FEATURE_COUNT_PATH}/all_samples_counts.txt"
    params:
        gtf_path = GTF_PATH
    threads: 4
    resources:
        mem_mb = 24000,
        cpus = 4,
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        module load subread/2.0.1
        featureCounts -s {config[params][feature_counts][strandness]} -p -t exon -g gene_id -T {threads} \
        -a {params.gtf_path} \
        -o {output.counts} {input.bam_files}
        module unload subread/2.0.1
        """


# Rule 5: Salmon quantification
rule salmon_quant:
    input:
        r1 = f"{TRIMMED_PATH}/{{sample}}_trimmed_1P.fq.gz",
        r2 = f"{TRIMMED_PATH}/{{sample}}_trimmed_2P.fq.gz"
    output:
        quant = f"{SALMON_PATH}/{{sample}}_salmon_quant/{{sample}}_quant.sf"
    params:
        salmon_index = SALMON_INDEX,
        output_dir = f"{SALMON_PATH}/{{sample}}_salmon_quant",
        temp_quant = f"{SALMON_PATH}/{{sample}}_salmon_quant/quant.sf"
    threads: 8
    resources:
        mem_mb = 24000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        module load salmon/1.8.0
        
        salmon quant -i {params.salmon_index} -l {config[params][salmon][library_type]} \
        -1 {input.r1} -2 {input.r2} \
        -p {threads} --validateMappings --gcBias \
        -o {params.output_dir}
        
        # Rename the quant.sf file
        mv {params.temp_quant} {output.quant}
        
        module unload salmon/1.8.0
        """


# Rule 6: MultiQC report
rule multiqc:
    input:
        expand(f"{TRIMMED_PATH}/{{sample}}_trimmed_1P.fq.gz", sample=SAMPLES),
        expand(f"{TRIMMED_PATH}/{{sample}}_trimmed_2P.fq.gz", sample=SAMPLES),
        expand(f"{HISAT2_PATH}/{{sample}}_align_sorted.bam", sample=SAMPLES),
        expand(f"{FEATURE_COUNT_PATH}/{{sample}}_counts.txt", sample=SAMPLES),
        expand(f"{SALMON_PATH}/{{sample}}_salmon_quant/{{sample}}_quant.sf", sample=SAMPLES)
    output:
        report = "multiqc_report.html"
    threads: 2
    resources:
        mem_mb = 4000,
        cpus = 2,
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        module load singularity/3.11.3
        singularity run /dfs9/ucightf-lab/kstachel/TOOLS/multiqc-1.20.sif multiqc . -o .
        module unload singularity/3.11.3
        """
    
# Rule 7: DESeq2 differential expression analysis
rule deseq2:
    input:
        counts = f"{FEATURE_COUNT_PATH}/all_samples_counts.txt",
        metadata = config["deseq2"]["metadata"]
    output:
        results = f"{DESEQ2_PATH}/deseq2_results.csv",
        rds = f"{DESEQ2_PATH}/dds.rds"
    params:
        out_dir = f"{DESEQ2_PATH}",
        deseq2_condition = config["isee_app"]["condition"],
        group_a = config["isee_app"]["group_a"],
        group_b = config["isee_app"]["group_b"]
    threads: 1
    resources:
        mem_mb = 8000,
        cpus = 1,
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        module load R/4.2.2
        Rscript src/deseq2_analysis.R {input.counts} {input.metadata} \
        {params.out_dir} {params.deseq2_condition} {params.group_a} {params.group_b}
        module unload R/4.2.2
        cp {output.rds} isee_uci/shiny-server/test_app/dds.rds
        """

# Rule 8: Generate parametric iSEE app2.R file
rule generate_isee_app:
    input:
        template = "templates/app.R.template",
        dds = f"{DESEQ2_PATH}/dds.rds"
    output:
        app = "isee_uci/shiny-server/test_app/app.R"
    params:
        deseq2_condition = config["isee_app"]["condition"],
        group_a = config["isee_app"]["group_a"],
        group_b = config["isee_app"]["group_b"]
    shell:
        """
        # Create the output directory if it doesn't exist
        mkdir -p $(dirname {output.app})
        
        # Replace placeholders in template with actual parameters
        sed -e 's|{{DESEQ2_CONDITION}}|{params.deseq2_condition}|g' \
            -e 's|{{GROUP_A}}|{params.group_a}|g' \
            -e 's|{{GROUP_B}}|{params.group_b}|g' \
            {input.template} > {output.app}
        """
