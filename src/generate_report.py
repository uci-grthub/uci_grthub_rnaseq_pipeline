#!/usr/bin/env python3
"""
Generate a project report for bulk RNA-seq analysis driven by Snakemake
with optional summaries from produced outputs (MultiQC, featureCounts, DESeq2).

USAGE:
    python generate_report.py [OPTIONS]

OPTIONS:
    --fastq-dir DIR           Path to FASTQ directory (sample subfolders supported)
                              Default: data/FASTQ

    --output FILE             Output PDF path
                              Default: RNAseq_Project_Report.pdf

    --author NAME             Report author name
                              Default: Kevin Stachelek

    --padj-threshold FLOAT    Padj threshold for DE gene counts
                              Default: 0.05

    --fast                    Fast mode: skip heavy scans (e.g., DESeq2 CSV padj counting)
                              Default: off

    --species LIST            Comma-separated species analysed; selects the per-species
                              references, count matrix and QC outputs to report on
                              Default: inferred from the metadata species column

OUTPUTS:
    - PDF report with project information, the reference genome and annotation used,
      pipeline details (FastQC, Trimmomatic, HISAT2, featureCounts, Salmon),
      per-sample alignment statistics from MultiQC, count matrix summary, sample
      correlation/clustering/PCA figures, DESeq2 contrasts with significant gene
      counts, a deliverable file index, and the NCBI submission package status.

REQUIREMENTS:
    - reportlab: for PDF generation
    - Python 3.8+
"""

import os
import glob
import gzip
import re
import json
import csv
from datetime import datetime
from pathlib import Path
from collections import defaultdict
from reportlab.lib.pagesizes import letter, A4, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (
    BaseDocTemplate,
    PageTemplate,
    Frame,
    NextPageTemplate,
    Table,
    TableStyle,
    Paragraph,
    Spacer,
    PageBreak,
    Image,
)
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT


class MetadataSummary:
    """Load and summarize project metadata from a CSV file.

    Column names vary between the sequencing-core sheet and the pipeline's own
    metadata.csv, so lookups accept several spellings (case/spacing-insensitive):
      - Sample Name / sample
      - i7 Barcode Sequence / i7 / i7barcode
      - i5 Barcode Sequence / i5 / i5barcode_NovaSeqV1.5
      - Organism / species (optional)
      - Sample ID, Condition, Sex, Strain, Description, Comments (optional)
    """

    def __init__(self, path: str | None):
        self.path = path
        self.rows = []  # normalized dicts
        self.index_to_sample = {}
        self.name_to_sample = {}
        self._load()

    @staticmethod
    def _first(norm: dict, *keys):
        """Return the first non-empty value among several possible column names."""
        for key in keys:
            val = norm.get(key)
            if val:
                return val
        return None

    @staticmethod
    def _norm(s: str) -> str:
        return re.sub(r"\s+", " ", s.strip())

    def _normalize_key(self, k: str) -> str:
        kk = k.strip().lower()
        kk = kk.replace(" ", "_")
        kk = kk.replace("(", "").replace(")", "")
        return kk

    def _pair(self, i7: str | None, i5: str | None) -> str | None:
        if not i7 or not i5:
            return None
        i7u = self._norm(i7).upper()
        i5u = self._norm(i5).upper()
        if not re.fullmatch(r"[ACGTN]+", i7u) or not re.fullmatch(r"[ACGTN]+", i5u):
            return None
        return f"{i7u}-{i5u}"

    def _load(self):
        if not self.path or not os.path.isfile(self.path):
            return
        try:
            with open(self.path, 'r', newline='') as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    norm = {self._normalize_key(k): (v.strip() if isinstance(v, str) else v) for k, v in row.items() if k}
                    sample = self._first(norm, 'sample_name', 'sample', 'name')
                    # metadata.csv spells these i7barcode / i5barcode_NovaSeqV1.5
                    i7 = self._first(norm, 'i7_barcode_sequence', 'i7barcode', 'i7')
                    i5 = self._first(norm, 'i5_barcode_sequence', 'i5barcode', 'i5')
                    if not i5:
                        i5 = next(
                            (v for k, v in norm.items() if k.startswith('i5barcode') and v),
                            None,
                        )
                    organism = self._first(norm, 'organism', 'species')
                    description = self._first(norm, 'description')
                    comments = self._first(norm, 'comments')
                    idx = self._pair(i7, i5)
                    rec = {
                        'sample_name': sample or '',
                        'sample_id': self._first(norm, 'sample_id') or '',
                        'i7': (i7 or '').upper(),
                        'i5': (i5 or '').upper(),
                        'organism': organism or '',
                        'condition': self._first(norm, 'condition', 'group') or '',
                        'sex': self._first(norm, 'sex') or '',
                        'strain': self._first(norm, 'strain', 'mouse_strain') or '',
                        'description': description or '',
                        'comments': comments or '',
                        'index_pair': idx or '',
                    }
                    self.rows.append(rec)
                    if idx:
                        self.index_to_sample[idx] = rec
                    if sample:
                        self.name_to_sample[sample] = rec
        except Exception as e:
            print(f"Warning: Failed to read metadata CSV {self.path}: {e}")


class FASTQMetadataExtractor:
    """Extract metadata from FASTQ files in layouts used by the Snakefile.

    The read-pair suffixes mirror the Snakefile's own list, so any layout the
    workflow can run is also one the report can describe. Files may sit
    directly under the FASTQ directory or in a per-sample subfolder.
    """

    # (R1 suffix, R2 suffix), checked in order
    READ_SUFFIXES = (
        ("-READ1-Sequences.txt.gz", "-READ2-Sequences.txt.gz"),
        ("-R1.fastq.gz", "-R2.fastq.gz"),
        ("_R1.fastq.gz", "_R2.fastq.gz"),
        ("_r1.fq.gz", "_r2.fq.gz"),
        ("_R1.fq.gz", "_R2.fq.gz"),
        ("-R1.fq.gz", "-R2.fq.gz"),
    )

    def __init__(self, fastq_dir: str, size_only: bool = True):
        self.fastq_dir = fastq_dir
        self.samples = {}
        self.size_only = size_only
        self.parse_samples()

    def parse_samples(self):
        """Parse sample names and file sizes from FASTQ directory."""
        r1_candidates = []
        for r1_suffix, _ in self.READ_SUFFIXES:
            r1_candidates += glob.glob(
                os.path.join(self.fastq_dir, "**", f"*{r1_suffix}"), recursive=True
            )
        r1_candidates = sorted(set(r1_candidates))

        for r1_file in r1_candidates:
            r1_base = os.path.basename(r1_file)
            r1_dir = os.path.dirname(r1_file)

            for r1_suffix, r2_suffix in self.READ_SUFFIXES:
                if r1_base.endswith(r1_suffix):
                    sample = r1_base[: -len(r1_suffix)]
                    r2_file = os.path.join(r1_dir, f"{sample}{r2_suffix}")
                    break
            else:
                continue

            r1_size = self._safe_size_gb(r1_file)
            r2_size = self._safe_size_gb(r2_file) if os.path.exists(r2_file) else 0.0

            i7_index, i5_index = self._infer_indices_from_sample(sample)
            # Size-only estimate by default (avoid decompressing gz on network FS)
            read_count = self._estimate_reads_by_size(r1_file)

            self.samples[sample] = {
                'r1_path': r1_file,
                'r2_path': r2_file if os.path.exists(r2_file) else None,
                'r1_size_gb': r1_size,
                'r2_size_gb': r2_size,
                'read_count': read_count,
                'total_size_gb': (r1_size or 0) + (r2_size or 0),
                'i7_index': i7_index,
                'i5_index': i5_index,
            }

    @staticmethod
    def _safe_size_gb(path: str) -> float:
        try:
            return os.path.getsize(path) / (1024 ** 3)
        except Exception:
            return 0.0

    @staticmethod
    def _infer_indices_from_sample(sample_name: str):
        # Typical format: xR074-L8-G3-P057-ATGTACCT-TAGGTATG -> last two parts are i7/i5
        parts = sample_name.split('-')
        i7 = parts[-2] if len(parts) >= 2 else 'N/A'
        i5 = parts[-1] if len(parts) >= 1 else 'N/A'
        # sanity: keep only plausible index sequences (A/C/G/T and length>=6)
        def clean(idx):
            if re.fullmatch(r"[ACGTN]+", idx or "") and len(idx) >= 6:
                return idx
            return 'N/A'
        return clean(i7), clean(i5)

    def _count_reads_sample(self, fastq_file, sample_lines=1000):
        """Estimate read count by sampling header lines from gzip FASTQ file.
        This is approximate and meant for quick reporting only.
        """
        # Retained for optional future use; not used by default.
        try:
            heads = 0
            with gzip.open(fastq_file, 'rt') as f:
                for _ in range(sample_lines):
                    h = f.readline()
                    if not h:
                        break
                    if h.startswith('@') and not h.startswith('@+'):
                        heads += 1

            total_size = os.path.getsize(fastq_file)
            # Fallback multiplier if file is very small or read failed
            if heads == 0:
                return 0
            # Assume ~100 bytes per read record chunk in gzip on average; heuristic
            approx_reads = int(total_size / 100)
            return max(approx_reads, heads)
        except Exception as e:
            print(f"Warning: Could not estimate reads in {fastq_file}: {e}")
            return 0

    def _estimate_reads_by_size(self, fastq_file: str) -> int:
        """Estimate reads using compressed file size only (very fast heuristic).
        Uses ~100 bytes per read record in gzip as a coarse average.
        """
        try:
            total_size = os.path.getsize(fastq_file)  # bytes (compressed)
            if total_size <= 0:
                return 0
            return int(total_size / 100)
        except Exception:
            return 0

    def get_summary(self):
        total_samples = len(self.samples)
        total_size = sum((s.get('total_size_gb') or 0) for s in self.samples.values())
        avg_size = (total_size / total_samples) if total_samples > 0 else 0

        return {
            'total_samples': total_samples,
            'total_size_gb': total_size,
            'avg_size_gb': avg_size,
            'generation_date': datetime.now().strftime("%B %d, %Y"),
        }


class MultiQCSummary:
    """Per-sample metrics pulled from multiqc_data/multiqc_data.json.

    MultiQC stores its general-statistics table as `report_general_stats_data`,
    a list of per-module dicts keyed by that module's own sample name. Those
    names carry tool-specific suffixes (`-R1`, `_summary`, `_salmon_quant`,
    `_align_sorted_markdup`), so each is normalised back to the pipeline sample
    name before the modules are merged.
    """

    # Suffixes each tool appends to the sample name, longest first so that
    # e.g. '_align_sorted_markdup' is stripped before '_align'.
    SAMPLE_SUFFIXES = (
        '_align_sorted_markdup',
        '_salmon_quant',
        '_summary',
        '-R1',
        '-R2',
        '_R1',
        '_R2',
    )

    def __init__(self, base_dir: str):
        self.base_dir = base_dir
        self.multiqc_json = os.path.join(base_dir, 'multiqc_data', 'multiqc_data.json')
        self.fastqc = {}
        self.hisat2 = {}
        self.trimmomatic = {}
        self.featurecounts = {}
        self.salmon = {}
        self._parse()

    @property
    def available(self) -> bool:
        return bool(self.fastqc or self.hisat2)

    @classmethod
    def _base_sample(cls, name: str) -> str:
        """Strip tool-specific suffixes to recover the pipeline sample name."""
        cleaned = name
        changed = True
        while changed:
            changed = False
            for suffix in cls.SAMPLE_SUFFIXES:
                if cleaned.endswith(suffix):
                    cleaned = cleaned[: -len(suffix)]
                    changed = True
        # Trimmomatic reports as '<sample>_<sample>-R1'; collapse the repeat
        half = len(cleaned) // 2
        if len(cleaned) % 2 == 1 and cleaned[:half] == cleaned[half + 1:]:
            cleaned = cleaned[:half]
        return cleaned

    def _parse(self):
        if not os.path.isfile(self.multiqc_json):
            return
        try:
            with open(self.multiqc_json, 'r') as fh:
                data = json.load(fh)
        except (OSError, json.JSONDecodeError) as e:
            print(f"Warning: Failed parsing {self.multiqc_json}: {e}")
            return

        blocks = data.get('report_general_stats_data') or []
        if isinstance(blocks, dict):  # older MultiQC layout
            blocks = [blocks]

        for block in blocks:
            if not isinstance(block, dict):
                continue
            for raw_sample, metrics in block.items():
                if not isinstance(metrics, dict):
                    continue
                sample = self._base_sample(raw_sample)
                # FastQC reports each read separately; keep R1 as the
                # representative record so counts are per read pair.
                if 'total_sequences' in metrics:
                    if raw_sample.endswith('-R2') or raw_sample.endswith('_R2'):
                        continue
                    self.fastqc[sample] = {
                        'percent_gc': metrics.get('percent_gc'),
                        'total_sequences': metrics.get('total_sequences'),
                        'avg_sequence_length': metrics.get('avg_sequence_length'),
                        'percent_duplicates': metrics.get('percent_duplicates'),
                    }
                if 'overall_alignment_rate' in metrics:
                    self.hisat2[sample] = {
                        'aligned': metrics.get('overall_alignment_rate'),
                        'total_reads': metrics.get('total_reads'),
                        'concordant': metrics.get('paired_aligned_one'),
                    }
                if 'surviving_pct' in metrics:
                    self.trimmomatic[sample] = {
                        'surviving_pct': metrics.get('surviving_pct'),
                        'input_read_pairs': metrics.get('input_read_pairs'),
                        'dropped_pct': metrics.get('dropped_pct'),
                    }
                if 'percent_assigned' in metrics:
                    self.featurecounts[sample] = {
                        'percent_assigned': metrics.get('percent_assigned'),
                        'assigned': metrics.get('Assigned'),
                    }
                if 'percent_mapped' in metrics:
                    self.salmon[sample] = {
                        'percent_mapped': metrics.get('percent_mapped'),
                    }


class FeatureCountsSummary:
    """Summarize featureCounts matrix shape and basic stats."""

    def __init__(self, counts_path: str):
        self.counts_path = counts_path
        self.genes = 0
        self.samples = 0
        self.header_samples = []
        self._scan()

    def _scan(self):
        if not os.path.isfile(self.counts_path):
            return
        try:
            with open(self.counts_path, 'r') as fh:
                # featureCounts prefixes the table with its own command line as
                # a '#' comment; skipping it keeps that line out of the header
                # and out of the gene count.
                reader = csv.reader(
                    (line for line in fh if not line.startswith('#')), delimiter='\t'
                )
                header = next(reader, None)
                if header and len(header) >= 7:
                    # featureCounts: first 6 columns are annotation/meta, counts from 7th
                    self.header_samples = header[6:]
                    self.samples = len(self.header_samples)
                for _ in reader:
                    self.genes += 1
        except Exception as e:
            print(f"Warning: Failed to read counts from {self.counts_path}: {e}")


class DESeq2ResultsSummary:
    """Collect DESeq2 result CSVs and summarize significant gene counts."""

    def __init__(self, deseq_dir: str, padj_thresh: float = 0.05, fast: bool = False):
        self.deseq_dir = deseq_dir
        self.padj_thresh = padj_thresh
        self.fast = fast
        self.contrasts = []  # list of dicts with name, sex, n_sig
        self.pca_pdfs = []
        self._scan()

    def _scan(self):
        if not os.path.isdir(self.deseq_dir):
            return
        # Find result CSVs recursively
        csv_paths = glob.glob(os.path.join(self.deseq_dir, '**', 'results_*.csv'), recursive=True)
        for csv_path in sorted(csv_paths):
            nm = os.path.basename(csv_path).replace('results_', '').replace('.csv', '')
            sex = 'unknown'
            # infer sex from parent folder name like sex_M
            parts = Path(csv_path).parts
            for p in parts:
                if p.startswith('sex_'):
                    sex = p.replace('sex_', '', 1) or 'unknown'
            if self.fast:
                n_sig = 'skipped (fast)'
            else:
                n_sig = self._count_sig(csv_path)
            parsed = self._parse_contrast_name(nm)
            self.contrasts.append({
                'name': nm,
                'sex': sex,
                'n_sig': n_sig,
                'path': csv_path,
                'comparison_type': parsed['comparison_type'],
                'comparison_text': parsed['comparison_text'],
            })

        # PCA PDFs are stored in results/pca_plots_*.pdf
        results_dir = os.path.join(os.path.dirname(self.deseq_dir), '..', 'results')
        # Normalize path
        results_dir = str(Path(results_dir).resolve())
        pca_candidates = glob.glob(os.path.join(results_dir, 'pca_plots*.pdf'))
        self.pca_pdfs = sorted(pca_candidates)

    @staticmethod
    def _parse_contrast_name(name: str) -> dict:
        """Convert contrast filename tokens into readable comparison descriptions."""
        patterns = [
            (
                r"^sex_([^_]+)_main_condition_(.+?)_vs_(.+)$",
                lambda m: {
                    'comparison_type': 'main_condition',
                    'comparison_text': f"Condition main effect: {m.group(2)} vs {m.group(3)}"
                },
            ),
            (
                r"^sex_([^_]+)_main_age_(.+?)_vs_(.+)$",
                lambda m: {
                    'comparison_type': 'main_age',
                    'comparison_text': f"Age main effect: {m.group(2)} vs {m.group(3)}"
                },
            ),
            (
                r"^sex_([^_]+)_age_(.+?)_condition_(.+?)_vs_(.+)$",
                lambda m: {
                    'comparison_type': 'condition_within_age',
                    'comparison_text': f"Condition within age {m.group(2)}: {m.group(3)} vs {m.group(4)}"
                },
            ),
            (
                r"^sex_([^_]+)_condition_(.+?)_age_(.+?)_vs_(.+)$",
                lambda m: {
                    'comparison_type': 'age_within_condition',
                    'comparison_text': f"Age within condition {m.group(2)}: {m.group(3)} vs {m.group(4)}"
                },
            ),
            (
                r"^sex_([^_]+)_interaction_(.+)$",
                lambda m: {
                    'comparison_type': 'interaction',
                    'comparison_text': f"Interaction term: {m.group(2)}"
                },
            ),
        ]

        for pat, fn in patterns:
            mt = re.match(pat, name)
            if mt:
                return fn(mt)

        return {
            'comparison_type': 'other',
            'comparison_text': name,
        }

    def _count_sig(self, csv_path: str) -> int:
        n = 0
        try:
            with open(csv_path, 'r') as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    padj = row.get('padj')
                    if padj is None or padj == '' or padj == 'NA':
                        continue
                    try:
                        if float(padj) < self.padj_thresh:
                            n += 1
                    except ValueError:
                        continue
        except Exception as e:
            print(f"Warning: failed reading {csv_path}: {e}")
        return n


class DESeq2ComparisonsSummary:
    """Read the exported DESeq2 comparison manifest from project root."""

    def __init__(self, csv_path: str):
        self.csv_path = csv_path
        self.rows = []
        self._scan()

    def _scan(self):
        if not self.csv_path or not os.path.isfile(self.csv_path):
            return
        try:
            with open(self.csv_path, 'r', newline='') as fh:
                reader = csv.DictReader(fh)
                self.rows = [row for row in reader]
        except Exception as e:
            print(f"Warning: failed reading DESeq2 comparisons CSV {self.csv_path}: {e}")


class CountMatrixSummary:
    """Summarize the cleaned gene-level count matrix written by count_matrix.R."""

    def __init__(self, counts_dir: str):
        self.counts_dir = counts_dir
        self.counts_csv = os.path.join(counts_dir, 'gene_counts.csv')
        self.metrics_csv = os.path.join(counts_dir, 'count_matrix_metrics.csv')
        self.genes = 0
        self.samples = 0
        self.sample_names = []
        self.metrics = []
        self._scan()

    @property
    def available(self) -> bool:
        return self.genes > 0

    def _scan(self):
        if os.path.isfile(self.counts_csv):
            try:
                with open(self.counts_csv, newline='') as fh:
                    reader = csv.reader(fh)
                    header = next(reader, None)
                    if header:
                        self.sample_names = header[1:]
                        self.samples = len(self.sample_names)
                    for _ in reader:
                        self.genes += 1
            except OSError as e:
                print(f"Warning: Failed to read {self.counts_csv}: {e}")

        if os.path.isfile(self.metrics_csv):
            try:
                with open(self.metrics_csv, newline='') as fh:
                    self.metrics = list(csv.DictReader(fh))
            except OSError as e:
                print(f"Warning: Failed to read {self.metrics_csv}: {e}")


class SampleQCSummary:
    """Collect the correlation / clustering / PCA outputs written by sample_qc.R."""

    FIGURES = [
        ('pca_plot.png', 'Principal component analysis (PC1 vs PC2)'),
        ('pca_scree_plot.png', 'Variance explained by each principal component'),
        ('sample_correlation_spearman_heatmap.png', 'Sample-sample Spearman correlation'),
        ('sample_clustering_dendrogram.png', 'Hierarchical clustering of samples'),
    ]

    def __init__(self, qc_dir: str):
        self.qc_dir = qc_dir
        self.metrics = []
        self.pca_variance = []
        self.correlation_range = None
        self.transformation = ''
        self.figures = []
        self._scan()

    @property
    def available(self) -> bool:
        return bool(self.metrics or self.figures)

    def _scan(self):
        if not os.path.isdir(self.qc_dir):
            return

        metrics_path = os.path.join(self.qc_dir, 'sample_metrics.csv')
        if os.path.isfile(metrics_path):
            try:
                with open(metrics_path, newline='') as fh:
                    self.metrics = list(csv.DictReader(fh))
            except OSError as e:
                print(f"Warning: Failed to read {metrics_path}: {e}")

        variance_path = os.path.join(self.qc_dir, 'pca_variance_explained.csv')
        if os.path.isfile(variance_path):
            try:
                with open(variance_path, newline='') as fh:
                    self.pca_variance = list(csv.DictReader(fh))
            except OSError as e:
                print(f"Warning: Failed to read {variance_path}: {e}")

        transform_path = os.path.join(self.qc_dir, 'transformation.txt')
        if os.path.isfile(transform_path):
            try:
                with open(transform_path) as fh:
                    self.transformation = fh.read().strip()
            except OSError:
                pass

        self.correlation_range = self._correlation_range(
            os.path.join(self.qc_dir, 'sample_correlation_spearman.csv')
        )

        for name, caption in self.FIGURES:
            path = os.path.join(self.qc_dir, name)
            if os.path.isfile(path):
                self.figures.append((path, caption))

    @staticmethod
    def _correlation_range(path: str):
        """Min/max off-diagonal correlation, i.e. how tightly the samples agree."""
        if not os.path.isfile(path):
            return None
        try:
            with open(path, newline='') as fh:
                reader = csv.reader(fh)
                next(reader, None)  # header
                values = []
                for row_idx, row in enumerate(reader):
                    for col_idx, cell in enumerate(row[1:]):
                        if col_idx == row_idx:
                            continue  # self-correlation is always 1
                        try:
                            values.append(float(cell))
                        except ValueError:
                            continue
        except OSError as e:
            print(f"Warning: Failed to read {path}: {e}")
            return None
        if not values:
            return None
        return min(values), max(values)


class NCBISubmissionSummary:
    """Report on the GEO/SRA submission package, if it has been generated."""

    def __init__(self, submission_dir: str):
        self.submission_dir = submission_dir
        self.geo_csv = os.path.join(submission_dir, 'geo_samples.csv')
        self.sra_csv = os.path.join(submission_dir, 'sra_metadata.csv')
        self.md5_txt = os.path.join(submission_dir, 'md5sums.txt')
        self.n_samples = 0
        self.n_checksums = 0
        self._scan()

    @property
    def available(self) -> bool:
        return os.path.isfile(self.geo_csv)

    def _scan(self):
        if os.path.isfile(self.geo_csv):
            try:
                with open(self.geo_csv, newline='') as fh:
                    self.n_samples = max(0, sum(1 for _ in fh) - 1)
            except OSError as e:
                print(f"Warning: Failed to read {self.geo_csv}: {e}")
        if os.path.isfile(self.md5_txt):
            try:
                with open(self.md5_txt) as fh:
                    self.n_checksums = sum(1 for line in fh if line.strip())
            except OSError as e:
                print(f"Warning: Failed to read {self.md5_txt}: {e}")


class ReportGenerator:
    """Generate PDF report summarizing pipeline inputs and outputs."""

    def __init__(self, output_path, author, fastq_dir, padj_thresh=0.05, workdir='.', fast: bool = False, metadata_path: str | None = None, comparisons_csv: str | None = None, species: list[str] | None = None):
        self.output_path = output_path
        self.author = author
        self.fastq_dir = fastq_dir
        self.padj_thresh = padj_thresh
        self.workdir = workdir
        self.fast = fast
        self.cfg = self._load_config()
        # Metadata is loaded before the species list so the species column can
        # inform it when --species was not supplied.
        self.metadata = MetadataSummary(metadata_path)
        self.species_list = species or self._infer_species()
        self.primary_species = self.species_list[0] if self.species_list else 'unknown'
        self.output_dir = self._cfg_get(['paths', 'output'], 'output') or 'output'

        self.extractor = FASTQMetadataExtractor(fastq_dir, size_only=True)
        self.summary = self.extractor.get_summary()
        self.mqc = MultiQCSummary(os.path.join(workdir, self.output_dir))
        self.counts = CountMatrixSummary(
            os.path.join(workdir, self.output_dir, 'counts', self.primary_species)
        )
        self.fc = FeatureCountsSummary(
            os.path.join(
                workdir, self.output_dir, 'feature_count',
                f'{self.primary_species}_samples_counts.txt',
            )
        )
        self.sample_qc = SampleQCSummary(
            os.path.join(workdir, self.output_dir, 'sample_qc', self.primary_species)
        )
        self.ncbi = NCBISubmissionSummary(
            os.path.join(workdir, self.output_dir, 'ncbi_submission', self.primary_species)
        )
        self.deseq = DESeq2ResultsSummary(os.path.join(workdir, self.output_dir, 'deseq2'), padj_thresh=padj_thresh, fast=fast)
        comparisons_path = comparisons_csv or 'deseq2_comparisons.csv'
        if not os.path.isabs(comparisons_path):
            comparisons_path = os.path.join(workdir, comparisons_path)
        self.deseq_manifest = DESeq2ComparisonsSummary(comparisons_path)

    def _infer_species(self):
        """Species actually analysed: metadata species column, else default_species."""
        default_species = self._cfg_get(['default_species'], 'human')
        found = {
            (row.get('organism') or '').strip()
            for row in self.metadata.rows
            if (row.get('organism') or '').strip()
        }
        return sorted(found) if found else [default_species]

    def generate(self):
        """Generate the PDF report."""
        # Set up a document with both portrait and landscape page templates
        left_margin = 0.75*inch
        right_margin = 0.75*inch
        top_margin = 0.75*inch
        bottom_margin = 0.75*inch

        # Frames for portrait and landscape
        pw, ph = letter
        lw, lh = landscape(letter)

        portrait_frame = Frame(
            left_margin,
            bottom_margin,
            pw - left_margin - right_margin,
            ph - top_margin - bottom_margin,
            id='portrait_frame'
        )
        landscape_frame = Frame(
            left_margin,
            bottom_margin,
            lw - left_margin - right_margin,
            lh - top_margin - bottom_margin,
            id='landscape_frame'
        )

        def on_portrait(canvas, doc):
            canvas.setPageSize(letter)

        def on_landscape(canvas, doc):
            canvas.setPageSize(landscape(letter))

        portrait_template = PageTemplate(id='Portrait', frames=[portrait_frame], onPage=on_portrait)
        landscape_template = PageTemplate(id='Landscape', frames=[landscape_frame], onPage=on_landscape)

        doc = BaseDocTemplate(
            self.output_path,
            pagesize=letter,
            rightMargin=right_margin,
            leftMargin=left_margin,
            topMargin=top_margin,
            bottomMargin=bottom_margin,
            pageTemplates=[portrait_template, landscape_template],
        )
        
        elements = []
        styles = getSampleStyleSheet()
        
        # Custom styles
        title_style = ParagraphStyle(
            'CustomTitle',
            parent=styles['Heading1'],
            fontSize=24,
            textColor=colors.HexColor('#1f4788'),
            spaceAfter=12,
            alignment=TA_CENTER,
            fontName='Helvetica-Bold'
        )
        cell_style_small = ParagraphStyle(
            'CellSmall',
            parent=styles['BodyText'],
            fontSize=8,
            leading=10,
            wordWrap='CJK'
        )
        
        heading_style = ParagraphStyle(
            'CustomHeading',
            parent=styles['Heading2'],
            fontSize=14,
            textColor=colors.HexColor('#1f4788'),
            spaceAfter=12,
            spaceBefore=12,
            fontName='Helvetica-Bold'
        )
        
        # Title
        elements.append(Paragraph("RNA-seq Project Report", title_style))
        elements.append(Spacer(1, 0.2*inch))
        
        # Project Information
        elements.append(Paragraph("Project Information", heading_style))
        project_info = [
            ['Generation Date:', self.summary['generation_date']],
            ['Author:', self.author]
        ]
        cpus_cfg = self._cfg_get(['params', 'cpus'])
        if cpus_cfg is not None:
            project_info.append(['Configured CPUs:', str(cpus_cfg)])
        
        info_table = Table(project_info, colWidths=[2.7*inch, 3.3*inch])
        info_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (0, -1), colors.HexColor('#E8EEF7')),
            ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
        ]))
        elements.append(info_table)
        elements.append(Spacer(1, 0.3*inch))
        
        # Summary Statistics
        elements.append(Paragraph("Summary Statistics", heading_style))
        summary_data = [
            ['Metric', 'Value'],
            ['Total Samples', str(self.summary['total_samples'])],
            ['Total Data Size', f"{self.summary['total_size_gb']:.2f} GB"],
            ['Average Size per Sample', f"{self.summary['avg_size_gb']:.2f} GB"],
        ]
        summary_data.append(['Organism(s)', ', '.join(self.species_list)])
        condition_counts = defaultdict(int)
        for row in self.metadata.rows:
            condition = (row.get('condition') or '').strip()
            if condition:
                condition_counts[condition] += 1
        if condition_counts:
            summary_data.append(['Experimental Groups', str(len(condition_counts))])
            replicates = sorted(set(condition_counts.values()))
            replicate_text = (
                str(replicates[0]) if len(replicates) == 1
                else f"{min(replicates)}-{max(replicates)}"
            )
            summary_data.append(['Replicates per Group', replicate_text])
        if self.counts.available:
            summary_data.append(['Genes Quantified', f"{self.counts.genes:,}"])

        summary_table = Table(summary_data, colWidths=[2.5*inch, 2.5*inch])
        summary_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1f4788')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F0F0F0')]),
        ]))
        elements.append(summary_table)
        elements.append(Spacer(1, 0.3*inch))

        body_style = ParagraphStyle(
            'CustomBody',
            parent=styles['BodyText'],
            fontSize=10,
            alignment=TA_LEFT,
            spaceAfter=10,
            leading=14
        )
        path_style = ParagraphStyle(
            'PathCell',
            parent=styles['BodyText'],
            fontSize=7,
            leading=9,
            fontName='Courier',
            wordWrap='CJK'
        )

        # Reference Genome and Annotation -- reported per species, since samples
        # are aligned against the reference matching their metadata species
        # rather than a single project-wide genome.
        elements.append(Paragraph("Reference Genome and Annotation", heading_style))
        for species in self.species_list:
            refs = self._cfg_get(['species_references', species]) or {}
            organism = refs.get('organism', species)
            ref_rows = [['Item', 'Value']]
            ref_rows.append(['Organism', str(organism)])
            if refs.get('genome_build'):
                ref_rows.append(['Genome build', str(refs['genome_build'])])
            if refs.get('annotation'):
                ref_rows.append(['Annotation release', str(refs['annotation'])])
            for label, key in (
                ('HISAT2 index', 'hisat2_index'),
                ('GTF annotation', 'gtf'),
                ('Salmon index', 'salmon_index'),
            ):
                value = refs.get(key) or self._cfg_get(['references', key])
                if value:
                    ref_rows.append([label, Paragraph(str(value), path_style)])
            adapters = self._cfg_get(['references', 'adapters'])
            if adapters:
                ref_rows.append(['Adapter sequences', Paragraph(str(adapters), path_style)])

            if len(self.species_list) > 1:
                elements.append(Paragraph(f"<b>{species}</b>", styles['Heading3']))
            ref_table = Table(ref_rows, colWidths=[1.7*inch, 5.3*inch], repeatRows=1)
            ref_table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1f4788')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTNAME', (0, 1), (0, -1), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 9),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
                ('TOPPADDING', (0, 0), (-1, -1), 5),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
                ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F5F5F5')]),
            ]))
            elements.append(ref_table)
            elements.append(Spacer(1, 0.2*inch))

        # Pipeline Overview
        elements.append(Paragraph("Pipeline Overview", heading_style))

        # FastQC
        elements.append(Paragraph("<b>FastQC v0.11.9</b>", styles['Heading3']))
        elements.append(Paragraph(
            "Quality control was performed on raw FASTQ files using FastQC. Per-sample HTML reports are aggregated by MultiQC.",
            body_style
        ))
        # Trimmomatic
        elements.append(Paragraph("<b>Trimmomatic v0.39</b>", styles['Heading3']))
        tr_adapters = self._display_file(self._cfg_get(['references', 'adapters']))
        tr_illum = self._cfg_get(['params', 'trimmomatic', 'illuminaclip'])
        tr_window = self._cfg_get(['params', 'trimmomatic', 'sliding_window'])
        tr_minlen = self._cfg_get(['params', 'trimmomatic', 'min_length'])
        if tr_adapters and tr_illum and tr_window is not None and tr_minlen is not None:
            tr_text = (
                f"Adapters were removed and reads trimmed in PE mode with "
                f"ILLUMINACLIP:{tr_adapters}:{tr_illum}, SLIDINGWINDOW:{tr_window}, MINLEN:{tr_minlen}."
            )
        else:
            tr_text = (
                "Adapters were removed and reads trimmed with Trimmomatic (PE mode) using project parameters."
            )
        elements.append(Paragraph(tr_text, body_style))
        # HISAT2
        elements.append(Paragraph("<b>HISAT2 v2.2.1 + SAMtools v1.15.1</b>", styles['Heading3']))
        hs_strand = self._cfg_get(['params', 'hisat2', 'rna_strandness'])
        hs_index = self._display_file(self._species_ref('hisat2_index'))
        if hs_strand or hs_index:
            hs_bits = []
            if hs_strand:
                hs_bits.append(f"strandness: {hs_strand}")
            if hs_index:
                hs_bits.append(f"index: {hs_index}")
            hs_text = (
                f"Trimmed reads were aligned with HISAT2 ({', '.join(hs_bits)}). "
                f"BAMs were sorted and indexed with SAMtools."
            )
        else:
            hs_text = (
                "Trimmed reads were aligned to the reference using HISAT2. BAMs were sorted and indexed with SAMtools."
            )
        elements.append(Paragraph(hs_text, body_style))
        # featureCounts
        elements.append(Paragraph("<b>featureCounts v2.0.1</b>", styles['Heading3']))
        fc_strand = self._cfg_get(['params', 'feature_counts', 'strandness'])
        fc_gtf = self._display_file(self._species_ref('gtf'))
        if fc_gtf or fc_strand is not None:
            fc_bits = ["-t exon", "-g gene_id", "-p"]
            if fc_strand is not None:
                fc_bits.append(f"-s {fc_strand}")
            if fc_gtf:
                fc_bits.append(f"-a {fc_gtf}")
            fc_text = (
                "Gene-level counts were generated across all samples with featureCounts ("
                + ", ".join(fc_bits)
                + ")."
            )
        else:
            fc_text = (
                "Gene-level counts were generated in a single run across all samples (exon features, gene_id attribute)."
            )
        elements.append(Paragraph(fc_text, body_style))
        # Salmon (optional)
        elements.append(Paragraph("<b>Salmon v1.8.0</b>", styles['Heading3']))
        sm_lib = self._cfg_get(['params', 'salmon', 'library_type'])
        sm_index = self._display_file(self._species_ref('salmon_index'))
        if sm_lib or sm_index:
            sm_bits = ["--validateMappings", "--gcBias"]
            if sm_lib:
                sm_bits.append(f"-l {sm_lib}")
            if sm_index:
                sm_bits.append(f"-i {sm_index}")
            sm_text = "Transcript-level quantification with Salmon (" + ", ".join(sm_bits) + ")."
        else:
            sm_text = (
                "Transcript-level quantification with Salmon was configured (validateMappings, gcBias); outputs are per-sample if enabled."
            )
        elements.append(Paragraph(sm_text, body_style))
        # MultiQC
        elements.append(Paragraph("<b>MultiQC v1.20</b>", styles['Heading3']))
        elements.append(Paragraph(
            "QC summaries were aggregated into a single HTML report. See multiqc_report.html for details.",
            body_style
        ))
        # Sample QC
        elements.append(Paragraph("<b>Sample QC (DESeq2 vst, R)</b>", styles['Heading3']))
        elements.append(Paragraph(
            "Raw gene-level counts were filtered to expressed genes and variance-stabilised, "
            "then used for sample-sample correlation, hierarchical clustering, and principal "
            "component analysis. See the sample QC section below.",
            body_style
        ))
        # DESeq2
        elements.append(Paragraph("<b>DESeq2 (R)</b>", styles['Heading3']))
        elements.append(Paragraph(
            "Differential expression is available on demand from the comparison definitions in "
            "the project comparisons config; it is not part of the default primary analysis "
            "target. Comparisons that have been run are listed below.",
            body_style
        ))

        # Alignment and mapping statistics, per sample, from the MultiQC summary
        elements.append(Paragraph("Alignment and Mapping Statistics", heading_style))
        alignment_rows = self._get_alignment_rows()
        if alignment_rows:
            align_data = [[
                'Sample',
                'Read Pairs',
                'GC %',
                'Surviving\nTrimming',
                'HISAT2\nAlignment',
                'Assigned to\nGenes',
                'Salmon\nMapping',
            ]]
            for row in alignment_rows:
                align_data.append([
                    Paragraph(row['sample'], cell_style_small),
                    row['total_sequences'],
                    row['percent_gc'],
                    row['surviving'],
                    row['aligned'],
                    row['assigned'],
                    row['salmon'],
                ])
            align_table = Table(
                align_data,
                colWidths=[1.0*inch, 1.1*inch, 0.6*inch, 1.0*inch, 1.1*inch, 1.1*inch, 1.1*inch],
                repeatRows=1,
            )
            align_table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1f4788')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (1, 0), (-1, -1), 'RIGHT'),
                ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 9),
                ('FONTSIZE', (0, 1), (-1, -1), 8),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
                ('TOPPADDING', (0, 0), (-1, -1), 4),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
                ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F5F5F5')]),
            ]))
            elements.append(align_table)
            elements.append(Paragraph(
                "Per-sample FastQC, trimming and alignment metrics in full are in "
                f"{self._relpath(os.path.join(self.output_dir, 'multiqc_report.html'))}.",
                body_style
            ))
        else:
            elements.append(Paragraph(
                "MultiQC data were not found, so per-sample alignment statistics could not be "
                f"tabulated here. See {self._relpath(os.path.join(self.output_dir, 'multiqc_report.html'))}.",
                body_style
            ))

        # Gene-level count matrix
        elements.append(Paragraph("Gene-Level Count Matrix", heading_style))
        if self.counts.available:
            counts_rel = self._relpath(self.counts.counts_csv)
            elements.append(Paragraph(
                f"Raw gene-level counts for {self.counts.samples} samples across "
                f"{self.counts.genes:,} annotated genes were produced by featureCounts and "
                f"written as a plain matrix to {counts_rel}. A CPM-normalised matrix and the "
                "gene annotation (chromosome, span, strand, union-exon length) accompany it in "
                f"{self._relpath(self.counts.counts_dir)}.",
                body_style
            ))
            if self.counts.metrics:
                count_data = [['Sample', 'Assigned Reads', '% Assigned', 'Genes Detected']]
                for row in self.counts.metrics:
                    count_data.append([
                        Paragraph(str(row.get('sample', '')), cell_style_small),
                        self._format_number(row.get('assigned_reads')),
                        self._format_number(row.get('percent_assigned')),
                        self._format_number(row.get('genes_detected')),
                    ])
                counts_table = Table(
                    count_data,
                    colWidths=[2.2*inch, 1.6*inch, 1.4*inch, 1.6*inch],
                    repeatRows=1,
                )
                counts_table.setStyle(TableStyle([
                    ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1f4788')),
                    ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                    ('ALIGN', (1, 0), (-1, -1), 'RIGHT'),
                    ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                    ('FONTSIZE', (0, 0), (-1, 0), 9),
                    ('FONTSIZE', (0, 1), (-1, -1), 8),
                    ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
                    ('TOPPADDING', (0, 0), (-1, -1), 4),
                    ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
                    ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F5F5F5')]),
                ]))
                elements.append(counts_table)
        elif self.fc.genes:
            elements.append(Paragraph(
                f"featureCounts produced {self.fc.genes:,} genes across {self.fc.samples} samples "
                f"in {self._relpath(self.fc.counts_path)}. The cleaned matrix has not been "
                "generated yet; run the gene_count_matrix rule.",
                body_style
            ))
        else:
            elements.append(Paragraph(
                "No count matrix was found. Run the feature_counts_all and gene_count_matrix rules.",
                body_style
            ))
        elements.append(Spacer(1, 0.15*inch))

        # Sample correlation, clustering and PCA
        elements.append(PageBreak())
        elements.append(Paragraph("Sample Correlation, Clustering and PCA", heading_style))
        if self.sample_qc.available:
            qc_text = (
                "Counts were filtered to expressed genes and transformed with the "
                f"{self.sample_qc.transformation or 'variance-stabilising transformation'} "
                "before computing sample-sample correlations, Euclidean distances and PCA. "
                "PCA uses the 500 most variable genes."
            )
            if self.sample_qc.correlation_range:
                low, high = self.sample_qc.correlation_range
                qc_text += (
                    f" Off-diagonal Spearman correlations range from {low:.3f} to {high:.3f}."
                )
            elements.append(Paragraph(qc_text, body_style))

            if self.sample_qc.pca_variance:
                var_data = [['Component', 'Variance Explained (%)']]
                for row in self.sample_qc.pca_variance[:5]:
                    var_data.append([
                        str(row.get('component', '')),
                        self._format_number(row.get('percent_variance')),
                    ])
                var_table = Table(var_data, colWidths=[2.0*inch, 2.5*inch])
                var_table.setStyle(TableStyle([
                    ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1f4788')),
                    ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                    ('ALIGN', (1, 0), (-1, -1), 'RIGHT'),
                    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                    ('FONTSIZE', (0, 0), (-1, -1), 9),
                    ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
                    ('TOPPADDING', (0, 0), (-1, -1), 4),
                    ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
                    ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F5F5F5')]),
                ]))
                elements.append(var_table)
                elements.append(Spacer(1, 0.2*inch))

            caption_style = ParagraphStyle(
                'Caption',
                parent=styles['BodyText'],
                fontSize=8,
                alignment=TA_CENTER,
                textColor=colors.HexColor('#444444'),
                spaceAfter=12,
            )
            for path, caption in self.sample_qc.figures:
                figure = self._scaled_image(path, max_width=6.0*inch, max_height=4.2*inch)
                if figure is None:
                    continue
                elements.append(figure)
                elements.append(Paragraph(caption, caption_style))
            elements.append(Paragraph(
                "Correlation matrices, distance matrices, PCA coordinates and vector (PDF) "
                f"versions of these figures are in {self._relpath(self.sample_qc.qc_dir)}.",
                body_style
            ))
        else:
            elements.append(Paragraph(
                "Sample QC outputs were not found. Run the sample_qc rule to generate the "
                "correlation, clustering and PCA results.",
                body_style
            ))

        # DESeq2 comparisons undertaken
        elements.append(Paragraph("DESeq2 Comparisons Undertaken", heading_style))
        comparison_rows = self._get_deseq_comparison_rows()
        if comparison_rows:
            de_rows = [[
                'Sex',
                'Comparison Type',
                'Comparison',
                f"Significant Genes (padj < {self.padj_thresh})",
            ]]

            for item in comparison_rows:
                de_rows.append([
                    str(item.get('sex', 'unknown')),
                    str(item.get('comparison_type', 'other')),
                    Paragraph(str(item.get('comparison_text', item.get('contrast_name', ''))), cell_style_small),
                    str(item.get('n_sig', 'N/A')),
                ])

            de_table = Table(
                de_rows,
                colWidths=[0.7*inch, 1.3*inch, 3.6*inch, 1.0*inch],
                repeatRows=1,
            )
            de_table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1f4788')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('ALIGN', (2, 1), (2, -1), 'LEFT'),
                ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 9),
                ('FONTSIZE', (0, 1), (-1, -1), 8),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
                ('TOPPADDING', (0, 0), (-1, -1), 6),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
                ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F5F5F5')]),
            ]))
            elements.append(de_table)
        else:
            elements.append(Paragraph(
                "No DESeq2 comparison manifest or result CSV files were found, so no undertaken comparisons could be listed.",
                body_style,
            ))
        elements.append(Spacer(1, 0.15*inch))
        
        # Start a new landscape page for Sample Details
        elements.append(NextPageTemplate('Landscape'))
        elements.append(PageBreak())
        elements.append(Paragraph("Sample Details", heading_style))

        # If metadata available, merge on i7/i5 indexes into a unified table
        merge_with_metadata = bool(self.metadata and self.metadata.rows)

        if merge_with_metadata:
            sample_data = [['Sample Name', 'Sample ID', 'Group', 'Strain', 'i7 Index', 'i5 Index', 'Read Pairs', 'Organism']]
        else:
            sample_data = [['Sample Name', 'i7 Index', 'i5 Index', 'Read Pairs']]

        for sample_name in sorted(self.extractor.samples.keys()):
            sample_info = self.extractor.samples[sample_name]
            i7_idx = sample_info['i7_index']
            i5_idx = sample_info['i5_index']
            # Prefer FastQC's exact count over the compressed-size heuristic
            exact_reads = self.mqc.fastqc.get(sample_name, {}).get('total_sequences')
            if exact_reads:
                est_reads = self._format_number(exact_reads)
            elif sample_info['read_count'] > 0:
                est_reads = f"{sample_info['read_count']:,.0f} (est.)"
            else:
                est_reads = 'N/A'

            if merge_with_metadata:
                # Match on the sample name first; the index pair is the fallback
                # for sheets that name samples differently from the FASTQ files.
                meta = self.metadata.name_to_sample.get(sample_name)
                if meta is None and i7_idx != 'N/A' and i5_idx != 'N/A':
                    key = f"{str(i7_idx).upper()}-{str(i5_idx).upper()}"
                    meta = self.metadata.index_to_sample.get(key)
                meta = meta or {}
                row = [
                    Paragraph(sample_name, cell_style_small),
                    Paragraph(meta.get('sample_id', ''), cell_style_small),
                    Paragraph(meta.get('condition', ''), cell_style_small),
                    Paragraph(meta.get('strain', ''), cell_style_small),
                    i7_idx,
                    i5_idx,
                    est_reads,
                    Paragraph(meta.get('organism', ''), cell_style_small),
                ]
            else:
                row = [
                    Paragraph(sample_name, cell_style_small),
                    i7_idx,
                    i5_idx,
                    est_reads,
                ]
            sample_data.append(row)

        if merge_with_metadata:
            # Columns sum to 9.5 inches, the landscape frame width
            samples_table = Table(
                sample_data,
                colWidths=[2.4*inch, 0.8*inch, 1.9*inch, 0.9*inch, 0.9*inch, 0.9*inch, 0.9*inch, 0.8*inch],
                repeatRows=1,
            )
        else:
            samples_table = Table(
                sample_data,
                colWidths=[2.2*inch, 1.25*inch, 1.25*inch, 1.3*inch],
                repeatRows=1,
            )
        samples_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1f4788')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('ALIGN', (0, 0), (0, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('WORDWRAP', (0, 0), (-1, -1), 'CJK'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 9),
            ('FONTSIZE', (0, 1), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
            ('TOPPADDING', (0, 0), (-1, -1), 6),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F5F5F5')]),
        ]))
        elements.append(samples_table)
        elements.append(Spacer(1, 0.3*inch))
        # Switch back to portrait for the remaining content
        elements.append(NextPageTemplate('Portrait'))

        # Processed data files and NCBI submission package
        elements.append(PageBreak())
        elements.append(Paragraph("Processed Data Files and NCBI Submission", heading_style))

        deliverables = [['Deliverable', 'Location']]
        for label, path in self._deliverable_paths():
            deliverables.append([label, Paragraph(path, path_style)])
        deliverable_table = Table(deliverables, colWidths=[2.6*inch, 4.4*inch], repeatRows=1)
        deliverable_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1f4788')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F5F5F5')]),
        ]))
        elements.append(deliverable_table)
        elements.append(Spacer(1, 0.2*inch))

        if self.ncbi.available:
            elements.append(Paragraph(
                f"A GEO/SRA submission package for {self.ncbi.n_samples} samples has been "
                f"assembled in {self._relpath(self.ncbi.submission_dir)}. It contains "
                "geo_samples.csv (paste into the SAMPLES section of the GEO metadata workbook), "
                "sra_metadata.csv (the SRA workbook), and md5sums.txt covering "
                f"{self.ncbi.n_checksums} raw and processed files. Raw FASTQ checksums are "
                "carried over from the checksum file supplied with the sequencing data rather "
                "than recomputed. SUBMISSION_README.txt in the same directory lists which files "
                "to upload.",
                body_style
            ))
        else:
            elements.append(Paragraph(
                "The NCBI submission package has not been generated yet. Run the "
                "ncbi_submission rule to produce the GEO and SRA metadata sheets and checksums.",
                body_style
            ))

        # References
        elements.append(PageBreak())
        elements.append(Paragraph("References", heading_style))
        
        reference_style = ParagraphStyle(
            'Reference',
            parent=styles['BodyText'],
            fontSize=9,
            leftIndent=0.2*inch,
            spaceAfter=8,
            leading=11,
            textColor=colors.black
        )
        
        references = [
            "Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. http://www.bioinformatics.babraham.ac.uk/projects/fastqc",
            "Bolger et al. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics.",
            "Kim et al. (2015). HISAT: a fast spliced aligner with low memory requirements. Nature Methods.",
            "Liao et al. (2014). featureCounts: assigning sequence reads to genomic features. Bioinformatics.",
            "Patro et al. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods.",
            "Love et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology.",
            "Ewels et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics.",
        ]
        
        for i, ref in enumerate(references, 1):
            elements.append(Paragraph(f"<b>{i}.</b> {ref}", reference_style))
        
        # Build PDF
        doc.build(elements)
        print(f"Report generated successfully: {self.output_path}")

    @staticmethod
    def _format_number(value):
        """Render a CSV field as a thousands-separated number, or pass it through."""
        if value is None or value == '' or value == 'NA':
            return 'N/A'
        try:
            number = float(value)
        except (TypeError, ValueError):
            return str(value)
        if number.is_integer():
            return f"{int(number):,}"
        return f"{number:,.2f}"

    def _sample_label(self, sample: str) -> str:
        """Short experimental label (GAC1) when metadata provides one."""
        meta = self.metadata.name_to_sample.get(sample) or {}
        return meta.get('sample_id') or sample

    def _get_alignment_rows(self):
        """Per-sample QC, trimming, alignment and assignment rates from MultiQC."""
        if not self.mqc.available:
            return []

        def pct(value):
            return 'N/A' if value is None else f"{self._format_number(value)}%"

        samples = sorted(set(self.mqc.fastqc) | set(self.mqc.hisat2))
        rows = []
        for sample in samples:
            fastqc = self.mqc.fastqc.get(sample, {})
            hisat2 = self.mqc.hisat2.get(sample, {})
            trim = self.mqc.trimmomatic.get(sample, {})
            featurecounts = self.mqc.featurecounts.get(sample, {})
            salmon = self.mqc.salmon.get(sample, {})
            rows.append({
                'sample': self._sample_label(sample),
                'total_sequences': self._format_number(fastqc.get('total_sequences')),
                'percent_gc': self._format_number(fastqc.get('percent_gc')),
                'surviving': pct(trim.get('surviving_pct')),
                'aligned': pct(hisat2.get('aligned')),
                'assigned': pct(featurecounts.get('percent_assigned')),
                'salmon': pct(salmon.get('percent_mapped')),
            })
        return rows

    def _scaled_image(self, path: str, max_width: float, max_height: float):
        """Load a figure scaled to fit the frame while preserving aspect ratio."""
        try:
            from reportlab.lib.utils import ImageReader
            width, height = ImageReader(path).getSize()
        except Exception as e:
            print(f"Warning: Could not read figure {path}: {e}")
            return None
        if not width or not height:
            return None
        scale = min(max_width / width, max_height / height)
        return Image(path, width=width * scale, height=height * scale)

    def _deliverable_paths(self):
        """Locations of the primary-analysis deliverables, for the report's file index."""
        species = self.primary_species
        out = self.output_dir
        entries = [
            ('Raw read QC (FastQC)', os.path.join(out, 'fastqc')),
            ('Aggregated QC report', os.path.join(out, 'multiqc_report.html')),
            ('Trimmed reads', os.path.join(out, 'trimmed')),
            ('Alignments (BAM + index)', os.path.join(out, 'hisat2_alignment')),
            ('Alignment summaries', os.path.join(out, 'hisat2_alignment', 'alignment_summary')),
            ('Post-alignment RNA QC', os.path.join(out, 'rustqc')),
            ('featureCounts output', os.path.join(out, 'feature_count', f'{species}_samples_counts.txt')),
            ('Raw gene count matrix', os.path.join(out, 'counts', species, 'gene_counts.csv')),
            ('CPM matrix', os.path.join(out, 'counts', species, 'gene_counts_cpm.csv')),
            ('Gene annotation', os.path.join(out, 'counts', species, 'gene_annotation.csv')),
            ('Transcript quantification (Salmon)', os.path.join(out, 'salmon')),
            ('Gene-level TPM matrix', os.path.join(out, 'tpm', species, 'tpm_salmon.csv')),
            ('Sample QC / correlation / PCA', os.path.join(out, 'sample_qc', species)),
            ('Alternative splicing (rMATS)', os.path.join(out, 'rmats', species)),
            ('NCBI submission package', os.path.join(out, 'ncbi_submission', species)),
            ('Sample metadata', self.metadata.path or 'metadata/metadata.csv'),
        ]
        return [(label, self._relpath(path)) for label, path in entries]

    def _relpath(self, path: str) -> str:
        try:
            return os.path.relpath(path, start=self.workdir)
        except Exception:
            return path

    def _display_file(self, path_val):
        """Return a display-friendly file identifier without absolute directories.
        - If a path-like string, return its basename.
        - Otherwise, return the original value.
        """
        if not path_val or not isinstance(path_val, str):
            return path_val
        # Normalize and strip trailing slash
        pv = path_val.rstrip('/').strip()
        # If it looks like a path (contains '/'), show only basename
        if '/' in pv:
            base = os.path.basename(pv)
            return base or pv
        return pv

    @staticmethod
    def _deep_merge(base: dict, overlay: dict) -> dict:
        merged = dict(base)
        for key, value in overlay.items():
            if isinstance(value, dict) and isinstance(merged.get(key), dict):
                merged[key] = ReportGenerator._deep_merge(merged[key], value)
            else:
                merged[key] = value
        return merged

    def _species_ref(self, key, default=None):
        """Reference path for the primary species, falling back to config references.

        Samples are aligned against per-species references, so reporting
        config['references'] verbatim would name the human genome for a mouse
        run whenever default_species has been overridden.
        """
        refs = self._cfg_get(['species_references', self.primary_species])
        if isinstance(refs, dict) and refs.get(key):
            return refs[key]
        return self._cfg_get(['references', key], default)

    def _cfg_get(self, keys, default=None):
        d = self.cfg or {}
        for k in keys:
            if isinstance(d, dict) and k in d:
                d = d[k]
            else:
                return default
        return d

    def _get_deseq_comparison_rows(self):
        n_sig_by_name = {
            str(item.get('name')): item.get('n_sig', 'N/A')
            for item in self.deseq.contrasts
        }

        if self.deseq_manifest.rows:
            rows = []
            for row in self.deseq_manifest.rows:
                contrast_name = row.get('contrast_name') or row.get('name') or ''
                rows.append({
                    'sex': row.get('sex', 'unknown'),
                    'comparison_type': row.get('comparison_type', 'other'),
                    'comparison_text': row.get('comparison_text', contrast_name),
                    'contrast_name': contrast_name,
                    'n_sig': n_sig_by_name.get(contrast_name, 'N/A'),
                })
            return sorted(
                rows,
                key=lambda x: (str(x.get('sex', '')), str(x.get('comparison_type', '')), str(x.get('contrast_name', '')))
            )

        return sorted(
            [
                {
                    'sex': item.get('sex', 'unknown'),
                    'comparison_type': item.get('comparison_type', 'other'),
                    'comparison_text': item.get('comparison_text', item.get('name', '')),
                    'contrast_name': item.get('name', ''),
                    'n_sig': item.get('n_sig', 'N/A'),
                }
                for item in self.deseq.contrasts
            ],
            key=lambda x: (str(x.get('sex', '')), str(x.get('comparison_type', '')), str(x.get('contrast_name', '')))
        )

    def _load_config(self):
        """Merge the same config files the Snakefile loads, in the same order.

        config.local.yaml overrides config.species_references.yaml overrides
        config.yaml. Loading only one of them (as this used to) drops every
        reference and parameter whenever a local override file exists, because
        that file holds a single key.
        """
        layers = [
            os.path.join(self.workdir, 'config.yaml'),
            os.path.join(self.workdir, 'config.species_references.yaml'),
            os.path.join(self.workdir, 'config.local.yaml'),
        ]
        present = [p for p in layers if os.path.isfile(p)]
        if not present:
            return None

        # Try PyYAML first
        try:
            import yaml  # type: ignore
            merged = {}
            for path in present:
                with open(path, 'r') as fh:
                    loaded = yaml.safe_load(fh) or {}
                merged = self._deep_merge(merged, loaded)
            return merged
        except ImportError:
            pass

        cfg_path = present[0]

        # Fallback: naive extraction of needed keys
        try:
            with open(cfg_path, 'r') as fh:
                text = fh.read()
        except Exception:
            return None

        def find_val(key, cast=None):
            m = re.search(rf"^\s*{re.escape(key)}\s*:\s*\"?([^#\n\"]+)\"?", text, re.MULTILINE)
            if not m:
                return None
            val = m.group(1).strip()
            if cast:
                try:
                    return cast(val)
                except Exception:
                    return val
            return val

        cfg = {
            'references': {
                'adapters': find_val('adapters'),
                'hisat2_index': find_val('hisat2_index'),
                'gtf': find_val('gtf'),
                'salmon_index': find_val('salmon_index'),
            },
            'params': {
                'trimmomatic': {
                    'illuminaclip': find_val('illuminaclip'),
                    'sliding_window': find_val('sliding_window'),
                    'min_length': find_val('min_length', int),
                },
                'hisat2': {
                    'rna_strandness': find_val('rna_strandness'),
                },
                'salmon': {
                    'library_type': find_val('library_type'),
                },
                'feature_counts': {
                    'strandness': find_val('strandness', int),
                },
                'cpus': find_val('cpus', int),
            }
        }
        return cfg


def main():
    """Main function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='Generate an RNA-seq project report (Snakemake + DESeq2)'
    )
    parser.add_argument(
        '--fastq-dir',
        default='data/FASTQ',
        help='Path to FASTQ directory (default: data/FASTQ)'
    )
    parser.add_argument(
        '--output',
        default='RNAseq_Project_Report.pdf',
        help='Output PDF path (default: RNAseq_Project_Report.pdf)'
    )
    parser.add_argument(
        '--author',
        default='Kevin Stachelek',
        help='Report author name'
    )
    parser.add_argument(
        '--padj-threshold',
        type=float,
        default=0.05,
        help='Padj threshold for counting significant genes (default: 0.05)'
    )
    parser.add_argument(
        '--fast',
        action='store_true',
        help='Fast mode: skip heavy scans (e.g., DESeq2 CSV padj counting)'
    )
    parser.add_argument(
        '--metadata',
        default='metadata/metadata.csv',
        help='Path to metadata CSV (default: metadata/metadata.csv)'
    )
    parser.add_argument(
        '--comparisons-csv',
        default='deseq2_comparisons.csv',
        help='Path to DESeq2 comparisons CSV (default: deseq2_comparisons.csv)'
    )
    parser.add_argument(
        '--species',
        default=None,
        help='Comma-separated species analysed (default: inferred from metadata/config)'
    )

    args = parser.parse_args()

    # Validate FASTQ directory
    if not os.path.isdir(args.fastq_dir):
        print(f"Error: FASTQ directory not found: {args.fastq_dir}")
        return 1

    # Generate report
    generator = ReportGenerator(
        output_path=args.output,
        author=args.author,
        fastq_dir=args.fastq_dir,
        padj_thresh=args.padj_threshold,
        workdir='.',
        fast=args.fast,
        metadata_path=args.metadata,
        comparisons_csv=args.comparisons_csv,
        species=[s.strip() for s in args.species.split(',') if s.strip()] if args.species else None
    )

    generator.generate()
    return 0


if __name__ == '__main__':
    exit(main())
