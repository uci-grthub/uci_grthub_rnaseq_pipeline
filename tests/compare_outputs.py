#!/usr/bin/env python3
"""Compare a pipeline run against reference outputs.

Shared by CI (tests/run_test.sh) and the HPC3 regression run, so both apply
the same tolerances:

  gene counts           exact match (featureCounts is deterministic)
  Salmon gene TPM       Pearson r > 0.999 per sample (EM varies with threads)
  QC metrics            within 0.1 percentage points: HISAT2 overall
                        alignment rate, samtools markdup duplication rate,
                        featureCounts assigned fraction
  rMATS                 same number of events per type. The default rule runs
                        rMATS with --b1 only, so there is no FDR to filter on;
                        event counts (fromGTF and junction-count tables) stand
                        in for "significant events" until a two-group design
                        is configured.
  report / GEO package  files exist; the PDF contains the expected section
                        headings; geo_samples.csv lists the same samples

Usage:
  compare_outputs.py RUN_DIR EXPECTED_DIR            compare, exit 1 on failure
  compare_outputs.py --write-expected RUN_DIR OUT    snapshot a run as reference

RUN_DIR is the project directory (holding output/ and the PDF report).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import sys
from pathlib import Path

import pandas as pd

TPM_MIN_PEARSON = 0.999
QC_TOLERANCE_PP = 0.1
RMATS_TYPES = ("SE", "RI", "MXE", "A5SS", "A3SS")
REPORT_NAME = "RNAseq_Project_Report.pdf"


# --- extraction ----------------------------------------------------------------


def species_dirs(output: Path) -> list[str]:
    return sorted(p.name for p in (output / "counts").iterdir() if p.is_dir())


def hisat2_rates(output: Path) -> dict[str, float]:
    rates = {}
    for path in sorted((output / "hisat2_alignment" / "alignment_summary").glob("*_summary.align")):
        match = re.search(r"([\d.]+)% overall alignment rate", path.read_text())
        if match:
            rates[path.name.removesuffix("_summary.align")] = float(match.group(1))
    return rates


def duplication_rates(output: Path) -> dict[str, float]:
    """Percent of mapped reads marked duplicate, from samtools markdup -f."""
    rates = {}
    for path in sorted((output / "hisat2_alignment").glob("*_markdup_metrics.txt")):
        fields = {}
        for line in path.read_text().splitlines():
            key, sep, value = line.partition(":")
            if sep:
                fields[key.strip()] = value.strip()
        total = float(fields.get("READ", "nan")) - float(fields.get("EXCLUDED", "0"))
        dups = float(fields.get("DUPLICATE TOTAL", "nan"))
        if total > 0:
            rates[path.name.removesuffix("_markdup_metrics.txt")] = 100 * dups / total
    return rates


def assigned_rates(output: Path) -> dict[str, float]:
    """Percent of featureCounts fragments assigned, per BAM column."""
    rates = {}
    for path in sorted((output / "feature_count").glob("*_samples_counts.txt.summary")):
        table = pd.read_csv(path, sep="\t", index_col=0)
        for column in table.columns:
            total = table[column].sum()
            library = Path(column).name.removesuffix("_align_sorted_markdup.bam")
            if total > 0:
                rates[library] = 100 * table.loc["Assigned", column] / total
    return rates


def qc_metrics(output: Path) -> pd.DataFrame:
    frame = pd.DataFrame(
        {
            "alignment_rate": pd.Series(hisat2_rates(output), dtype=float),
            "duplication_rate": pd.Series(duplication_rates(output), dtype=float),
            "assigned_rate": pd.Series(assigned_rates(output), dtype=float),
        }
    )
    frame.index.name = "library"
    return frame.sort_index()


def rmats_events(output: Path) -> dict[str, dict[str, int]]:
    events = {}
    for species_dir in sorted((output / "rmats").glob("*")):
        if not species_dir.is_dir():
            continue
        per_type = {}
        for event in RMATS_TYPES:
            for kind, pattern in (("fromGTF", f"fromGTF.{event}.txt"), ("JC", f"{event}.MATS.JC.txt")):
                path = species_dir / pattern
                if path.exists():
                    with path.open() as fh:
                        per_type[f"{kind}.{event}"] = max(sum(1 for _ in fh) - 1, 0)
        events[species_dir.name] = per_type
    return events


def report_text(run_dir: Path) -> str:
    from pypdf import PdfReader

    reader = PdfReader(str(run_dir / REPORT_NAME))
    return "\n".join(page.extract_text() or "" for page in reader.pages)


def report_sections(run_dir: Path) -> list[str]:
    """Numbered section headings ("3. Methods") in the PDF report."""
    headings = []
    for line in report_text(run_dir).splitlines():
        line = line.strip()
        if re.fullmatch(r"\d+(\.\d+)*\.?\s+[A-Z][\w /&()+-]{2,60}", line) and line not in headings:
            headings.append(line)
    return headings


def geo_samples(output: Path, species: str) -> list[str]:
    path = output / "ncbi_submission" / species / "geo_samples.csv"
    with path.open(newline="") as fh:
        rows = list(csv.reader(fh))
    return sorted(row[0] for row in rows[1:] if row and row[0])


# --- snapshot --------------------------------------------------------------------


def write_expected(run_dir: Path, out: Path) -> None:
    output = run_dir / "output"
    out.mkdir(parents=True, exist_ok=True)
    summary: dict[str, object] = {"species": {}}
    for species in species_dirs(output):
        dest = out / species
        dest.mkdir(exist_ok=True)
        shutil.copy(output / "counts" / species / "gene_counts.csv", dest / "gene_counts.csv")
        shutil.copy(output / "tpm" / species / "tpm_salmon.csv", dest / "tpm_salmon.csv")
        geo = output / "ncbi_submission" / species / "geo_samples.csv"
        summary["species"][species] = {"geo_samples": geo_samples(output, species) if geo.exists() else None}
    qc_metrics(output).to_csv(out / "qc_metrics.tsv", sep="\t", float_format="%.4f")
    summary["rmats_events"] = rmats_events(output)
    summary["report_sections"] = report_sections(run_dir)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"wrote reference outputs to {out}")


# --- comparison ----------------------------------------------------------------


class Checks:
    def __init__(self) -> None:
        self.failures: list[str] = []
        self.passed = 0

    def check(self, ok: bool, label: str, detail: str = "") -> None:
        if ok:
            self.passed += 1
            print(f"PASS  {label}")
        else:
            self.failures.append(label)
            print(f"FAIL  {label}{': ' + detail if detail else ''}")


def compare_counts(checks: Checks, run: Path, expected: Path) -> None:
    got = pd.read_csv(run, index_col=0).sort_index()
    want = pd.read_csv(expected, index_col=0).sort_index()
    same_shape = got.shape == want.shape and list(got.columns) == list(want.columns)
    if not same_shape:
        checks.check(False, f"gene counts {run}", f"shape {got.shape} vs {want.shape}")
        return
    diff = (got != want).to_numpy().sum()
    checks.check(diff == 0 and list(got.index) == list(want.index), f"gene counts exact ({run.parent.name})", f"{diff} cells differ")


def compare_tpm(checks: Checks, run: Path, expected: Path) -> None:
    got = pd.read_csv(run, index_col=0)
    want = pd.read_csv(expected, index_col=0)
    genes = got.index.intersection(want.index)
    checks.check(len(genes) == len(want.index), f"TPM gene set ({run.parent.name})", f"{len(genes)} of {len(want.index)} genes shared")
    for sample in want.columns:
        if sample not in got.columns:
            checks.check(False, f"TPM {sample}", "sample missing")
            continue
        r = got.loc[genes, sample].corr(want.loc[genes, sample], method="pearson")
        ok = not math.isnan(r) and r > TPM_MIN_PEARSON
        checks.check(ok, f"TPM Pearson {sample} r={r:.6f}", f"needs > {TPM_MIN_PEARSON}")


def compare_qc(checks: Checks, got: pd.DataFrame, want: pd.DataFrame) -> None:
    for library in want.index:
        for metric in want.columns:
            expected_value = want.loc[library, metric]
            value = got[metric].get(library, float("nan")) if metric in got else float("nan")
            delta = abs(value - expected_value)
            ok = not math.isnan(delta) and delta <= QC_TOLERANCE_PP
            checks.check(ok, f"{metric} {library} {value:.3f} vs {expected_value:.3f}", f"|diff| {delta:.3f} > {QC_TOLERANCE_PP} pp")


def compare(run_dir: Path, expected_dir: Path) -> int:
    output = run_dir / "output"
    summary = json.loads((expected_dir / "summary.json").read_text())
    checks = Checks()

    for species, info in summary["species"].items():
        compare_counts(checks, output / "counts" / species / "gene_counts.csv", expected_dir / species / "gene_counts.csv")
        compare_tpm(checks, output / "tpm" / species / "tpm_salmon.csv", expected_dir / species / "tpm_salmon.csv")
        if info.get("geo_samples") is not None:
            geo_dir = output / "ncbi_submission" / species
            for name in ("geo_samples.csv", "sra_metadata.csv", "md5sums.txt", "SUBMISSION_README.txt"):
                checks.check((geo_dir / name).is_file(), f"GEO package file {species}/{name}")
            if (geo_dir / "geo_samples.csv").is_file():
                checks.check(geo_samples(output, species) == info["geo_samples"], f"GEO samples {species}")

    want_qc = pd.read_csv(expected_dir / "qc_metrics.tsv", sep="\t", index_col=0)
    compare_qc(checks, qc_metrics(output), want_qc)

    got_events = rmats_events(output)
    for species, per_type in summary.get("rmats_events", {}).items():
        for key, count in per_type.items():
            value = got_events.get(species, {}).get(key)
            checks.check(value == count, f"rMATS {species} {key} {value} vs {count}")

    report = run_dir / REPORT_NAME
    checks.check(report.is_file(), "PDF report generated")
    if report.is_file():
        text = report_text(run_dir)
        for heading in summary.get("report_sections", []):
            checks.check(heading in text, f"report section '{heading}'")
    checks.check((output / "versions.yaml").is_file(), "versions.yaml written")

    print(f"\n{checks.passed} passed, {len(checks.failures)} failed")
    return 1 if checks.failures else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--write-expected", action="store_true", help="snapshot RUN_DIR into the second path")
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("expected_dir", type=Path)
    args = parser.parse_args()
    if args.write_expected:
        write_expected(args.run_dir, args.expected_dir)
        return 0
    return compare(args.run_dir, args.expected_dir)


if __name__ == "__main__":
    sys.exit(main())
