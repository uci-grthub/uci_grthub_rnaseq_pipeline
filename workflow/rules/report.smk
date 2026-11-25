# Run provenance (resolved config, versions.yaml) and the PDF project report.

import json
import subprocess


def _pipeline_version():
    version_file = REPO_DIR / "version.txt"
    return version_file.read_text().strip() if version_file.exists() else "unknown"


def _git_commit():
    """Commit of the workflow checkout, with -dirty when it has local edits."""
    try:
        commit = subprocess.run(
            ["git", "-C", str(REPO_DIR), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            [
                "git",
                "-C",
                str(REPO_DIR),
                "status",
                "--porcelain",
                "--untracked-files=no",
            ],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return f"{commit}-dirty" if dirty else commit


def _file_sha256(path):
    try:
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    except OSError:
        return "missing"


PIPELINE_VERSION = _pipeline_version()
PIPELINE_COMMIT = _git_commit()
PIXI_LOCK_SHA256 = _file_sha256(REPO_DIR / "pixi.lock")


# The config this run resolved (defaults + references + project config +
# command line), so the report and later readers see exactly what was used.
rule resolved_config:
    output:
        f"{OUTPUT_DIR}/pipeline_config.yaml",
    localrule: True
    params:
        cfg=json.dumps(config, sort_keys=True, default=str),
    run:
        import yaml

        with open(output[0], "w") as fh:
            yaml.safe_dump(json.loads(params.cfg), fh, sort_keys=False)


# Pipeline version, commit, pixi lock hash, tool versions and the references
# of the species analysed. Reruns whenever any of those change.
rule versions:
    input:
        config=f"{OUTPUT_DIR}/pipeline_config.yaml",
    output:
        f"{OUTPUT_DIR}/versions.yaml",
    log:
        "logs/versions/versions.log",
    envmodules:
        "fastqc/0.11.9",
        "hisat2/2.2.1",
        "samtools/1.15.1",
        "subread/2.0.1",
        "salmon/1.8.0",
        "rMATS/4.3.0",
        "R/4.5.2",
    params:
        version=PIPELINE_VERSION,
        commit=PIPELINE_COMMIT,
        lock_sha256=PIXI_LOCK_SHA256,
        software="env-modules" if USE_ENVMODULES else "pixi",
        species=",".join(SPECIES_LIST),
        env_tools=pixi_env("tools", has_modules=True),
        env_r=pixi_env("r", has_modules=True),
        env_default=pixi_env("default"),
    shell:
        """
        exec 2>{log}
        tmp=$(mktemp)
        trap 'rm -f "$tmp"' EXIT
        (
            {params.env_tools}
            echo "fastqc: $(fastqc --version 2>&1 | head -n1)"
            echo "hisat2: $(hisat2 --version 2>&1 | head -n1 | sed 's/.*version //')"
            echo "samtools: $(samtools --version 2>&1 | head -n1)"
            echo "featureCounts: $(featureCounts -v 2>&1 | grep -m1 -o 'v[0-9.]*')"
            echo "salmon: $(salmon --version 2>&1 | head -n1)"
            echo "rmats: $(rmats.py --version 2>&1 | head -n1)"
        ) >>"$tmp"
        (
            {params.env_r}
            echo "R: $(Rscript -e 'cat(R.version.string)' 2>&1 | tail -n1)"
            for pkg in DESeq2 tximport ensembldb; do
                echo "$pkg: $(Rscript -e "cat(as.character(packageVersion('$pkg')))" 2>&1 | tail -n1)"
            done
        ) >>"$tmp"
        (
            {params.env_tools}
            echo "chelae: $(chelae --version 2>&1 | head -n1)"
            echo "multiqc: $(multiqc --version 2>&1 | grep -m1 -o 'version [0-9.]*')"
        ) >>"$tmp"
        {params.env_default}
        python3 {SCRIPTS}/write_versions.py \
            --config {input.config} \
            --tool-versions "$tmp" \
            --version "{params.version}" \
            --commit "{params.commit}" \
            --pixi-lock-sha256 "{params.lock_sha256}" \
            --software "{params.software}" \
            --species "{params.species}" \
            --output {output}
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
        metadata=METADATA_PATH,
        config=f"{OUTPUT_DIR}/pipeline_config.yaml",
    output:
        report="RNAseq_Project_Report.pdf",
    log:
        "logs/generate_report/generate_report.log",
    benchmark:
        "benchmarks/generate_report/generate_report.tsv"
    threads: 1
    resources:
        mem_mb=4000,
        cpus=1,
    params:
        species=",".join(SPECIES_LIST),
        deseq2_flag=(
            "--include-deseq2"
            if config.get("report", {}).get("include_deseq2", True)
            else "--no-include-deseq2"
        ),
        env=pixi_env("default"),
    shell:
        """
        exec >{log} 2>&1
        {params.env}
        python3 {SCRIPTS}/generate_report.py \
            --fastq-dir {DATA_PATH} \
            --metadata {input.metadata} \
            --config {input.config} \
            --species {params.species} \
            {params.deseq2_flag} \
            --output {output.report}
        """
