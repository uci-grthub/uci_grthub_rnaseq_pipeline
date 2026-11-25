# Shared setup: config loading and validation, sample discovery, the sample
# sheet, per-species references, lane -> library grouping, and the software
# environment helper. Sample/species/library logic is unchanged from v0.9.0.

import csv
import glob
import hashlib
import os
import re
import sys
from pathlib import Path

import pandas as pd
from snakemake.utils import validate

REPO_DIR = Path(workflow.basedir).parent
SCRIPTS = str(Path(workflow.basedir) / "scripts")
RESOURCES = str(Path(workflow.basedir) / "resources")
PIXI_MANIFEST = str(REPO_DIR / "pixi.toml")


# Defaults and per-species references ship with the workflow; a project's
# --configfile is applied on top of them by Snakemake.
configfile: str(REPO_DIR / "config" / "config.yaml")
configfile: str(REPO_DIR / "config" / "references.yaml")


validate(config, str(Path(workflow.basedir) / "schemas" / "config.schema.yaml"))
for _species, _refs in config["species_references"].items():
    validate(_refs, str(Path(workflow.basedir) / "schemas" / "references.schema.yaml"))


# --- software environments ---------------------------------------------------
# Rules run their tools from a pixi environment (pixi.toml). Rules that also
# declare `envmodules:` use the HPC3 modules instead when Snakemake is run with
# `--software-deployment-method env-modules`.
try:
    from snakemake.settings.types import DeploymentMethod

    USE_ENVMODULES = (
        DeploymentMethod.ENV_MODULES in workflow.deployment_settings.deployment_method
    )
except (ImportError, AttributeError):
    USE_ENVMODULES = False


def pixi_env(name, has_modules=False):
    """Shell prefix that activates pixi environment `name`.

    Empty when the rule's own `envmodules:` will provide the tools.
    """
    if has_modules and USE_ENVMODULES:
        return ":"
    # Activation scripts of some packages read unset variables, and Snakemake
    # runs shell blocks with `set -u`, so relax it while activating.
    return (
        f'set +u; eval "$(pixi shell-hook --frozen --manifest-path {PIXI_MANIFEST} '
        f'-e {name} --shell bash)"; set -u'
    )


# --- run notifications ---------------------------------------------------------
NOTIFY_EMAIL = config.get("notify_email") or ""

if NOTIFY_EMAIL:

    onstart:
        shell("mail -s 'STARTED' {NOTIFY_EMAIL} < {log}")

    onsuccess:
        shell("mail -s 'DONE' {NOTIFY_EMAIL} < {log}")

    onerror:
        shell("mail -s 'ERROR' {NOTIFY_EMAIL} < {log}")


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
except OSError:
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


# Create a default sample sheet if none exists yet, so that rules depending
# on it (generate_report, deseq2) resolve during DAG building, including
# `snakemake -n`. When infer_samples_from_fastq is true (default), rows are
# pre-populated with sample names and barcodes inferred from FASTQ filenames;
# sex/condition are left blank for the user to fill in before a real run.
METADATA_PATH = config["samples"]
INFER_METADATA_FROM_FASTQ = config.get("infer_samples_from_fastq", True)
if not os.path.exists(METADATA_PATH):
    os.makedirs(os.path.dirname(METADATA_PATH) or ".", exist_ok=True)
    with open(METADATA_PATH, "w") as fh:
        fh.write("sample,i7barcode,i5barcode_NovaSeqV1.5,sex,condition,species\n")
        if INFER_METADATA_FROM_FASTQ:
            for sample in SAMPLES:
                i7, i5 = _infer_barcodes_from_sample(sample)
                fh.write(f"{sample},{i7},{i5},,,\n")

validate(
    pd.read_csv(METADATA_PATH, dtype=str, keep_default_na=False),
    str(Path(workflow.basedir) / "schemas" / "samples.schema.yaml"),
)

# Per-sample species -> reference paths, so samples from different organisms
# (e.g. mouse libraries mixed into an otherwise-human run) get aligned and
# quantified against the correct genome (config/references.yaml).
# There is deliberately no built-in default: a macaque library aligned to
# GRCh38 still clears the alignment-rate gate, so a blank species must stop
# the run rather than silently pick a genome. A project can opt in to one by
# setting default_species in its config.yaml.
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
        f"{sorted(SPECIES_REFERENCES)}), or set default_species in the project config.yaml."
    )

_unknown = sorted(
    {species for species in SAMPLE_SPECIES.values()} - set(SPECIES_REFERENCES)
)
if _unknown:
    raise ValueError(
        f"Unknown species {_unknown} in {METADATA_PATH} or default_species; "
        f"species_references defines {sorted(SPECIES_REFERENCES)}."
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
