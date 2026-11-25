# Gene-level counting, the delivered count matrices, sample QC and the
# on-demand DESeq2 rule.


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
    log:
        "logs/feature_counts_all/{species}.log",
    benchmark:
        "benchmarks/feature_counts_all/{species}.tsv"
    envmodules:
        "subread/2.0.1",
    threads: 4
    resources:
        mem_mb=24000,
        cpus=4,
    params:
        gtf_path=lambda wildcards: SPECIES_REFERENCES[wildcards.species]["gtf"],
        env=pixi_env("tools", has_modules=True),
    shell:
        """
        exec >{log} 2>&1
        {params.env}
        featureCounts -s {config[params][feature_counts][strandness]} -p -t exon -g gene_id -T {threads} \
            -a {params.gtf_path} \
            -o {output.counts} {input.bam_files}
        """


# Rule 4a: Clean gene-level count matrix. featureCounts writes a command-line
# comment, six annotation columns and BAM paths as column headers, none of
# which are usable as a delivered count matrix, so split it into plain
# counts/annotation/CPM tables keyed by the metadata sample_id.
rule gene_count_matrix:
    input:
        counts=f"{OUTPUT_DIR}/feature_count/{{species}}_samples_counts.txt",
        metadata=METADATA_PATH,
    output:
        counts_csv=f"{OUTPUT_DIR}/counts/{{species}}/gene_counts.csv",
        counts_rds=f"{OUTPUT_DIR}/counts/{{species}}/gene_counts.rds",
        cpm_csv=f"{OUTPUT_DIR}/counts/{{species}}/gene_counts_cpm.csv",
        annotation=f"{OUTPUT_DIR}/counts/{{species}}/gene_annotation.csv",
        metrics=f"{OUTPUT_DIR}/counts/{{species}}/count_matrix_metrics.csv",
    log:
        "logs/gene_count_matrix/{species}.log",
    benchmark:
        "benchmarks/gene_count_matrix/{species}.tsv"
    envmodules:
        "R/4.5.2",
    threads: 2
    resources:
        mem_mb=16000,
        cpus=2,
    params:
        out_dir=f"{OUTPUT_DIR}/counts/{{species}}",
        ensdb=lambda wildcards: SPECIES_REFERENCES[wildcards.species].get("ensdb", ""),
        env=pixi_env("r", has_modules=True),
    shell:
        """
        exec >{log} 2>&1
        {params.env}
        Rscript {SCRIPTS}/count_matrix.R {input.counts} {input.metadata} {params.out_dir} "{params.ensdb}"
        """


# Rule 4c: Sample-level QC -- sequencing/count metrics, sample-sample
# correlation, hierarchical clustering and PCA on variance-stabilised counts.
rule sample_qc:
    input:
        counts_csv=f"{OUTPUT_DIR}/counts/{{species}}/gene_counts.csv",
        metadata=METADATA_PATH,
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
    log:
        "logs/sample_qc/{species}.log",
    benchmark:
        "benchmarks/sample_qc/{species}.tsv"
    envmodules:
        "R/4.5.2",
    threads: 2
    resources:
        mem_mb=16000,
        cpus=2,
    params:
        out_dir=f"{OUTPUT_DIR}/sample_qc/{{species}}",
        env=pixi_env("r", has_modules=True),
    shell:
        """
        exec >{log} 2>&1
        {params.env}
        Rscript {SCRIPTS}/sample_qc.R {input.counts_csv} {input.metadata} {params.out_dir}
        """


# Rule 9: DESeq2 differential expression analysis (per species -- gene IDs
# and comparisons don't carry across organisms, so this is not a combined run).
# Not part of the default target; request it explicitly.
rule deseq2:
    input:
        counts=f"{OUTPUT_DIR}/feature_count/{{species}}_samples_counts.txt",
        metadata=METADATA_PATH,
        comparisons_config=config["deseq2"]["comparisons_config"],
    output:
        results=f"{OUTPUT_DIR}/deseq2/{{species}}/deseq2_results.csv",
        rds=f"{OUTPUT_DIR}/deseq2/{{species}}/dds.rds",
        manifest=f"{OUTPUT_DIR}/deseq2/{{species}}/deseq2_comparisons_manifest.csv",
        targets=f"{OUTPUT_DIR}/deseq2/{{species}}/target_gene_results.csv",
        per_region=f"{OUTPUT_DIR}/deseq2/{{species}}/target_per_region_log2fc.csv",
    log:
        "logs/deseq2/{species}.log",
    benchmark:
        "benchmarks/deseq2/{species}.tsv"
    envmodules:
        "R/4.5.2",
    threads: 1
    resources:
        mem_mb=8000,
        cpus=1,
    params:
        out_dir=f"{OUTPUT_DIR}/deseq2/{{species}}",
        ensdb=lambda wildcards: SPECIES_REFERENCES[wildcards.species].get("ensdb", ""),
        gtf=lambda wildcards: SPECIES_REFERENCES[wildcards.species]["gtf"],
        env=pixi_env("r", has_modules=True),
    shell:
        """
        exec >{log} 2>&1
        {params.env}
        Rscript {SCRIPTS}/deseq2_analysis.R {input.counts} {input.metadata} \
            {params.out_dir} {input.comparisons_config} "{params.ensdb}" "{params.gtf}"
        """
