# NCBI GEO/SRA submission package. Not part of the default target.


# Rule 8: NCBI GEO/SRA submission package. Checksums for the raw FASTQ files
# are taken from <paths.data>/md5sums.txt when the sequencing core supplied one,
# so this does not re-hash the full raw data.
rule ncbi_submission:
    input:
        metadata=METADATA_PATH,
        counts_csv=f"{OUTPUT_DIR}/counts/{{species}}/gene_counts.csv",
    output:
        geo=f"{OUTPUT_DIR}/ncbi_submission/{{species}}/geo_samples.csv",
        sra=f"{OUTPUT_DIR}/ncbi_submission/{{species}}/sra_metadata.csv",
        md5=f"{OUTPUT_DIR}/ncbi_submission/{{species}}/md5sums.txt",
        readme=f"{OUTPUT_DIR}/ncbi_submission/{{species}}/SUBMISSION_README.txt",
    log:
        "logs/ncbi_submission/{species}.log",
    benchmark:
        "benchmarks/ncbi_submission/{species}.tsv"
    threads: 1
    resources:
        mem_mb=4000,
        cpus=1,
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
        env=pixi_env("default"),
    shell:
        """
        exec >{log} 2>&1
        {params.env}
        python3 {SCRIPTS}/prepare_ncbi_submission.py \
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
