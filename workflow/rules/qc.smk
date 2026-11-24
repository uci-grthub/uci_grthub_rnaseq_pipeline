# Raw-read QC (FastQC), post-alignment RNA QC (RustQC) and the MultiQC report.


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
    params:
        out_dir=f"{OUTPUT_DIR}/fastqc/{{sample}}",
        env=pixi_env("tools", has_modules=True),
    envmodules:
        "fastqc/0.11.9",
    log:
        "logs/fastqc/{sample}.log",
    benchmark:
        "benchmarks/fastqc/{sample}.tsv"
    shell:
        """
        exec > {log} 2>&1
        {params.env}
        rm -rf {params.out_dir}
        mkdir -p {params.out_dir}
        fastqc -o {params.out_dir} -t {threads} {input.r1} {input.r2}
        """


# Rule 0b: RustQC RNA QC on aligned BAM files. Uses the Apptainer/Singularity
# image named by tools.rustqc_container when that file exists (HPC3: the
# bioconda rustqc build needs glibc >= 2.38, newer than the cluster's), and
# otherwise rustqc from the pixi `tools` env (CI, newer hosts).
rule rustqc:
    input:
        bam=f"{OUTPUT_DIR}/hisat2_alignment/{{library}}_align_sorted_markdup.bam",
    output:
        done=f"{OUTPUT_DIR}/rustqc/{{library}}/.done",
    threads: 8
    resources:
        mem_mb=24000,
        cpus=config["params"]["cpus"],
    params:
        out_dir=f"{OUTPUT_DIR}/rustqc/{{library}}",
        gtf_path=lambda wildcards: library_ref(wildcards.library, "gtf"),
        container=config["tools"]["rustqc_container"],
        env=pixi_env("tools"),
    log:
        "logs/rustqc/{library}.log",
    benchmark:
        "benchmarks/rustqc/{library}.tsv"
    shell:
        """
        exec > {log} 2>&1
        if [ -n "{params.container}" ] && [ -f "{params.container}" ]; then
            if command -v module >/dev/null 2>&1; then module load singularity/3.11.3; fi
            RUSTQC="singularity exec {params.container} rustqc"
        else
            {params.env}
            RUSTQC=rustqc
        fi
        mkdir -p {params.out_dir}
        $RUSTQC rna \
            --gtf {params.gtf_path} \
            --outdir {params.out_dir} \
            --threads {threads} \
            --paired \
            --skip-dup-check \
            {input.bam}
        touch {output.done}
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
        config=f"{RESOURCES}/multiqc_config.yaml",
    output:
        report=f"{OUTPUT_DIR}/multiqc_report.html",
    threads: 2
    resources:
        mem_mb=4000,
        cpus=2,
    params:
        env=pixi_env("tools"),
    log:
        "logs/multiqc/multiqc.log",
    benchmark:
        "benchmarks/multiqc/multiqc.tsv"
    shell:
        """
        exec > {log} 2>&1
        {params.env}
        rm -f {OUTPUT_DIR}/multiqc_report.html {OUTPUT_DIR}/multiqc_report_1.html
        rm -rf {OUTPUT_DIR}/multiqc_data {OUTPUT_DIR}/multiqc_data_1
        multiqc {OUTPUT_DIR} -o {OUTPUT_DIR} -c {input.config} --force
        """
