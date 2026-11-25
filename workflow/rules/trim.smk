# Adapter and quality trimming.


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
    log:
        f"{OUTPUT_DIR}/trimmed/{{sample}}_chelae.log",
    benchmark:
        "benchmarks/chelae_trim/{sample}.tsv"
    threads: 8
    resources:
        mem_mb=4000,
        cpus=config["params"]["cpus"],
    params:
        env=pixi_env("tools"),
    shell:
        """
        exec >{log} 2>&1
        {params.env}
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
