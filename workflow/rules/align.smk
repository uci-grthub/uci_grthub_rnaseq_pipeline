# HISAT2 alignment, coordinate sorting and duplicate marking.


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
    params:
        hisat2_index=lambda wildcards: library_ref(wildcards.library, "hisat2_index"),
        summary_path=f"{OUTPUT_DIR}/hisat2_alignment/alignment_summary",
        r1_list=lambda wildcards, input: ",".join(input.r1),
        r2_list=lambda wildcards, input: ",".join(input.r2),
        env=pixi_env("tools", has_modules=True),
    envmodules:
        "hisat2/2.2.1",
        "samtools/1.15.1",
    log:
        "logs/hisat2_align/{library}.log",
    benchmark:
        "benchmarks/hisat2_align/{library}.tsv"
    shell:
        """
        exec 2> {log}
        {params.env}

        hisat2 -p {threads} -t --qc-filter --rna-strandness {config[params][hisat2][rna_strandness]} \
            --summary-file {output.summary} \
            -x {params.hisat2_index} --dta-cufflinks \
            --rg-id {wildcards.library} --rg SM:{wildcards.library} \
            --rg LB:{wildcards.library} --rg PL:ILLUMINA \
            -1 {params.r1_list} -2 {params.r2_list} \
            | samtools sort -n -@ 2 \
            | samtools fixmate -m -@ 2 - {output.bam}
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
    params:
        env=pixi_env("tools", has_modules=True),
    envmodules:
        "samtools/1.15.1",
    log:
        "logs/sort_bam/{library}.log",
    benchmark:
        "benchmarks/sort_bam/{library}.tsv"
    shell:
        """
        exec > {log} 2>&1
        {params.env}

        samtools sort -@ {threads} -o {output.sorted_bam} {input.bam}
        samtools index -@ {threads} {output.sorted_bam}
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
    params:
        env=pixi_env("tools", has_modules=True),
    envmodules:
        "samtools/1.15.1",
    log:
        "logs/markdup/{library}.log",
    benchmark:
        "benchmarks/markdup/{library}.tsv"
    shell:
        """
        exec > {log} 2>&1
        {params.env}

        samtools markdup -@ {threads} -f {output.metrics} {input.sorted_bam} {output.markdup_bam}
        samtools index -@ {threads} {output.markdup_bam}
        """
