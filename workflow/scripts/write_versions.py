#!/usr/bin/env python3
"""Write output/versions.yaml: what produced this run's results.

Called by the `versions` rule; every value is passed in so the script itself
has no knowledge of the workflow layout.
"""

import argparse
import datetime
import hashlib
from pathlib import Path

import yaml

REFERENCE_KEYS = (
    "organism",
    "genome_build",
    "annotation",
    "hisat2_index",
    "gtf",
    "salmon_index",
    "ensdb",
)


def sha256(path: str) -> str:
    try:
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    except OSError:
        return "missing"


def parse_tool_versions(path: str) -> dict[str, str]:
    tools = {}
    for line in Path(path).read_text().splitlines():
        name, sep, value = line.partition(": ")
        if sep:
            tools[name.strip()] = value.strip()
    return tools


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="resolved pipeline config YAML")
    parser.add_argument("--tool-versions", required=True, help="'name: version' lines")
    parser.add_argument("--version", required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--pixi-lock-sha256", required=True)
    parser.add_argument("--software", required=True, help="pixi or env-modules")
    parser.add_argument("--species", required=True, help="comma-separated")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.config) as fh:
        config = yaml.safe_load(fh) or {}

    references = {}
    for species in filter(None, args.species.split(",")):
        refs = config.get("species_references", {}).get(species, {})
        entry = {key: refs[key] for key in REFERENCE_KEYS if refs.get(key)}
        if refs.get("gtf"):
            entry["gtf_sha256"] = sha256(refs["gtf"])
        references[species] = entry

    versions = {
        "pipeline": {
            "name": "uci_grthub_rnaseq_pipeline",
            "version": args.version,
            "commit": args.commit,
        },
        "software": {
            "source": args.software,
            "pixi_lock_sha256": args.pixi_lock_sha256,
            "tools": parse_tool_versions(args.tool_versions),
        },
        "references": references,
        "written": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
    }
    with open(args.output, "w") as fh:
        yaml.safe_dump(versions, fh, sort_keys=False)


if __name__ == "__main__":
    main()
