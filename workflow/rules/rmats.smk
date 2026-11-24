# Alternative splicing (rMATS). Not part of the default target.


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
    params:
        bam_list=f"{OUTPUT_DIR}/rmats/{{species}}/bam_files.txt",
        output_dir=f"{OUTPUT_DIR}/rmats/{{species}}",
        gtf_path=lambda wildcards: SPECIES_REFERENCES[wildcards.species]["gtf"],
        read_length=config["params"]["rmats"]["read_length"],
        env=pixi_env("tools", has_modules=True),
    envmodules:
        "rMATS/4.3.0",
    log:
        "logs/rmats/{species}.log",
    benchmark:
        "benchmarks/rmats/{species}.tsv"
    shell:
        """
        exec > {log} 2>&1
        {params.env}

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
        """
