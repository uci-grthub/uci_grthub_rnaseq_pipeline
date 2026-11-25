# Transcript quantification (Salmon) and tximport TPM/count matrices.


# Rule 5: Salmon quantification, one run per library over all of its lanes
# (salmon takes space-separated lists for -1/-2, kept in the same lane order),
# so transcript matrices have the same library columns as the gene counts
rule salmon_quant:
    input:
        r1=lane_files(f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_1P.fq.gz"),
        r2=lane_files(f"{OUTPUT_DIR}/trimmed/{{sample}}_trimmed_2P.fq.gz"),
    output:
        quant=f"{OUTPUT_DIR}/salmon/{{library}}_salmon_quant/{{library}}_quant.sf",
    log:
        "logs/salmon_quant/{library}.log",
    benchmark:
        "benchmarks/salmon_quant/{library}.tsv"
    envmodules:
        "salmon/1.8.0",
    threads: 8
    resources:
        # decoy-aware index is ~17 GB on disk; measured peak RSS ~18.5 GB
        mem_mb=24000,
        cpus=config["params"]["cpus"],
    params:
        salmon_index=lambda wildcards: library_ref(wildcards.library, "salmon_index"),
        output_dir=f"{OUTPUT_DIR}/salmon/{{library}}_salmon_quant",
        temp_quant=f"{OUTPUT_DIR}/salmon/{{library}}_salmon_quant/quant.sf",
        env=pixi_env("tools", has_modules=True),
    shell:
        """
        exec >{log} 2>&1
        {params.env}

        salmon quant -i {params.salmon_index} -l {config[params][salmon][library_type]} \
            -1 {input.r1} -2 {input.r2} \
            -p {threads} --validateMappings --gcBias \
            -o {params.output_dir} \
            --allowDovetail

        # Rename the quant.sf file
        mv {params.temp_quant} {output.quant}
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
    log:
        "logs/tximport_tpm/{species}.log",
    benchmark:
        "benchmarks/tximport_tpm/{species}.tsv"
    envmodules:
        "R/4.5.2",
    threads: 2
    resources:
        mem_mb=8000,
        cpus=2,
    params:
        staged_salmon_dir=f"{OUTPUT_DIR}/salmon_by_species/{{species}}",
        tpm_dir=f"{OUTPUT_DIR}/tpm/{{species}}",
        gtf_path=lambda wildcards: SPECIES_REFERENCES[wildcards.species]["gtf"],
        ensdb=lambda wildcards: SPECIES_REFERENCES[wildcards.species].get("ensdb", ""),
        env=pixi_env("r", has_modules=True),
    shell:
        """
        exec >{log} 2>&1
        rm -rf {params.staged_salmon_dir}
        mkdir -p {params.staged_salmon_dir}
        for quant_file in {input.quant_files}; do
            ln -s "$(readlink -f "$(dirname "$quant_file")")" \
                "{params.staged_salmon_dir}/$(basename "$(dirname "$quant_file")")"
        done

        {params.env}
        Rscript {SCRIPTS}/tximport_tpm.R {params.staged_salmon_dir} {params.gtf_path} {params.tpm_dir} "{params.ensdb}"
        """
