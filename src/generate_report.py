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

OUTPUTS:
    - PDF report with project information, pipeline details (FastQC, Trimmomatic,
      HISAT2, featureCounts, Salmon), optional QC/alignment summaries from
      MultiQC, counts matrix summary, DESeq2 contrasts with significant gene counts,
      and references to generated files.

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

    Expected columns (case/spacing-insensitive):
      - Sample Name
      - i7 Barcode Sequence
      - i5 Barcode Sequence
      - Organism (optional)
      - Description (optional)
      - Comments (optional)
    """

    def __init__(self, path: str | None):
        self.path = path
        self.rows = []  # normalized dicts
        self.index_to_sample = {}
        self._load()

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
                    norm = {self._normalize_key(k): (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
                    sample = norm.get('sample_name') or norm.get('sample') or norm.get('name')
                    i7 = norm.get('i7_barcode_sequence') or norm.get('i7')
                    i5 = norm.get('i5_barcode_sequence') or norm.get('i5')
                    organism = norm.get('organism')
                    description = norm.get('description')
                    comments = norm.get('comments')
                    idx = self._pair(i7, i5)
                    rec = {
                        'sample_name': sample or '',
                        'i7': (i7 or '').upper(),
                        'i5': (i5 or '').upper(),
                        'organism': organism or '',
                        'description': description or '',
                        'comments': comments or '',
                        'index_pair': idx or '',
                    }
                    self.rows.append(rec)
                    if idx:
                        self.index_to_sample[idx] = rec
        except Exception as e:
            print(f"Warning: Failed to read metadata CSV {self.path}: {e}")


class FASTQMetadataExtractor:
    """Extract metadata from FASTQ files in layouts used by the Snakefile.

    Supported patterns:
      - data/FASTQ/<sample>/<sample>-READ1-Sequences.txt.gz (and READ2)
      - data/FASTQ/**/*_r1.fq.gz (and _r2) [fallback]
    """

    def __init__(self, fastq_dir: str, size_only: bool = True):
        self.fastq_dir = fastq_dir
        self.samples = {}
        self.size_only = size_only
        self.parse_samples()

    def parse_samples(self):
        """Parse sample names and file sizes from FASTQ directory."""
        # Primary layout: subfolders per sample with READ1/READ2 naming
        r1_candidates = sorted(
            glob.glob(os.path.join(self.fastq_dir, "*", "*-READ1-Sequences.txt.gz"))
        )

        # Fallback: recursive search for common R1 patterns
        if not r1_candidates:
            r1_candidates = sorted(
                glob.glob(os.path.join(self.fastq_dir, "**", "*_r1.fq.gz"), recursive=True)
                + glob.glob(os.path.join(self.fastq_dir, "**", "*_R1.fq.gz"), recursive=True)
            )

        for r1_file in r1_candidates:
            r1_base = os.path.basename(r1_file)
            r1_dir = os.path.dirname(r1_file)

            if r1_base.endswith("-READ1-Sequences.txt.gz"):
                sample = r1_base.replace("-READ1-Sequences.txt.gz", "")
                r2_file = os.path.join(r1_dir, f"{sample}-READ2-Sequences.txt.gz")
            elif r1_base.endswith("_r1.fq.gz"):
                sample = r1_base[:-len("_r1.fq.gz")]
                r2_file = os.path.join(r1_dir, f"{sample}_r2.fq.gz")
            elif r1_base.endswith("_R1.fq.gz"):
                sample = r1_base[:-len("_R1.fq.gz")]
                r2_file = os.path.join(r1_dir, f"{sample}_R2.fq.gz")
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
    """Parse minimal stats from multiqc_data/multiqc_data.json if present."""

    def __init__(self, base_dir: str):
        self.base_dir = base_dir
        self.multiqc_json = os.path.join(base_dir, 'multiqc_data', 'multiqc_data.json')
        self.fastqc = {}
        self.hisat2 = {}
        self.trimmomatic = {}
        self._parse()

    def _parse(self):
        if not os.path.isfile(self.multiqc_json):
            return
        try:
            with open(self.multiqc_json, 'r') as fh:
                data = json.load(fh)
        except Exception as e:
            print(f"Warning: Failed parsing {self.multiqc_json}: {e}")
            return

        # FastQC module summary
        for rec in data.get('report_data', {}).get('fastqc', {}).get('general_stats', []):
            sample = rec.get('Sample') or rec.get('sample_name') or rec.get('Name')
            if not sample:
                continue
            self.fastqc[sample] = {
                'percent_gc': rec.get('percent_gc'),
                'total_sequences': rec.get('total_sequences'),
                'avg_sequence_length': rec.get('avg_sequence_length'),
            }

        # HISAT2 alignment; MultiQC often stores under 'hisat2'
        for rec in data.get('report_data', {}).get('hisat2', {}).get('general_stats', []):
            sample = rec.get('Sample') or rec.get('sample_name') or rec.get('Name')
            if not sample:
                continue
            self.hisat2[sample] = {
                'aligned': rec.get('hisat2_aligned', rec.get('alignment_rate')),
                'concordant': rec.get('hisat2_concordant_pairs'),
            }

        # Trimmomatic
        for rec in data.get('report_data', {}).get('trimmomatic', {}).get('general_stats', []):
            sample = rec.get('Sample') or rec.get('sample_name') or rec.get('Name')
            if not sample:
                continue
            self.trimmomatic[sample] = {
                'surviving_pct': rec.get('trimmomatic_surviving') or rec.get('percent_surviving'),
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
                reader = csv.reader(fh, delimiter='\t')
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
                if p.startswith('sex_') and len(p) == 5:
                    sex = p[-1]
            if self.fast:
                n_sig = 'skipped (fast)'
            else:
                n_sig = self._count_sig(csv_path)
            self.contrasts.append({'name': nm, 'sex': sex, 'n_sig': n_sig, 'path': csv_path})

        # PCA PDFs are stored in results/pca_plots_*.pdf
        results_dir = os.path.join(os.path.dirname(self.deseq_dir), '..', 'results')
        # Normalize path
        results_dir = str(Path(results_dir).resolve())
        pca_candidates = glob.glob(os.path.join(results_dir, 'pca_plots*.pdf'))
        self.pca_pdfs = sorted(pca_candidates)

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


class ReportGenerator:
    """Generate PDF report summarizing pipeline inputs and outputs."""

    def __init__(self, output_path, author, fastq_dir, padj_thresh=0.05, workdir='.', fast: bool = False, metadata_path: str | None = None):
        self.output_path = output_path
        self.author = author
        self.fastq_dir = fastq_dir
        self.padj_thresh = padj_thresh
        self.workdir = workdir
        self.fast = fast
        self.cfg = self._load_config()
        self.extractor = FASTQMetadataExtractor(fastq_dir, size_only=True)
        self.summary = self.extractor.get_summary()
        self.mqc = MultiQCSummary(workdir)
        self.fc = FeatureCountsSummary(os.path.join(workdir, 'output', 'feature_count', 'all_samples_counts.txt'))
        self.deseq = DESeq2ResultsSummary(os.path.join(workdir, 'output', 'deseq2'), padj_thresh=padj_thresh, fast=fast)
        # Metadata
        self.metadata = MetadataSummary(metadata_path)

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
        
        # Pipeline Overview
        elements.append(Paragraph("Pipeline Overview", heading_style))
        
        body_style = ParagraphStyle(
            'CustomBody',
            parent=styles['BodyText'],
            fontSize=10,
            alignment=TA_LEFT,
            spaceAfter=10,
            leading=14
        )
        
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
        elements.append(Paragraph("<b>HISAT2 v2.2.1 + SAMtools v1.10</b>", styles['Heading3']))
        hs_strand = self._cfg_get(['params', 'hisat2', 'rna_strandness'])
        hs_index = self._display_file(self._cfg_get(['references', 'hisat2_index']))
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
        fc_gtf = self._display_file(self._cfg_get(['references', 'gtf']))
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
        sm_index = self._display_file(self._cfg_get(['references', 'salmon_index']))
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
        # DESeq2
        elements.append(Paragraph("<b>DESeq2 (R)</b>", styles['Heading3']))
        elements.append(Paragraph(
            "Differential expression was performed with a design including condition, age, and their interaction. Analyses were run separately by sex. PCA plots were generated for QC.",
            body_style
        ))
        elements.append(Spacer(1, 0.15*inch))
        
        # Start a new landscape page for Sample Details
        elements.append(NextPageTemplate('Landscape'))
        elements.append(PageBreak())
        elements.append(Paragraph("Sample Details", heading_style))

        # If metadata available, merge on i7/i5 indexes into a unified table
        merge_with_metadata = bool(self.metadata and self.metadata.rows)

        if merge_with_metadata:
            sample_data = [['Sample Name', 'i7 Index', 'i5 Index', 'Est. Reads', 'Organism', 'Description']]
        else:
            sample_data = [['Sample Name', 'i7 Index', 'i5 Index', 'Est. Reads']]

        for sample_name in sorted(self.extractor.samples.keys()):
            sample_info = self.extractor.samples[sample_name]
            i7_idx = sample_info['i7_index']
            i5_idx = sample_info['i5_index']
            est_reads = f"{sample_info['read_count']:,.0f}" if sample_info['read_count'] > 0 else 'N/A'

            if merge_with_metadata and i7_idx != 'N/A' and i5_idx != 'N/A':
                key = f"{str(i7_idx).upper()}-{str(i5_idx).upper()}"
                meta = self.metadata.index_to_sample.get(key)
                organism = meta.get('organism', '') if meta else ''
                desc = meta.get('description', '') if meta else ''
                row = [
                    Paragraph(sample_name, cell_style_small),
                    i7_idx,
                    i5_idx,
                    est_reads,
                    Paragraph(organism, cell_style_small),
                    Paragraph(desc, cell_style_small),
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
            # Columns sum to 7 inches, widen Description to improve wrapping
            samples_table = Table(
                sample_data,
                colWidths=[1.8*inch, 0.9*inch, 0.9*inch, 1.0*inch, 1.0*inch, 1.4*inch],
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

    def _cfg_get(self, keys, default=None):
        d = self.cfg or {}
        for k in keys:
            if isinstance(d, dict) and k in d:
                d = d[k]
            else:
                return default
        return d

    def _load_config(self):
        """Load config.local.yaml if present, else config.yaml. Avoid hard dependency on PyYAML."""
        candidates = [
            os.path.join(self.workdir, 'config.local.yaml'),
            os.path.join(self.workdir, 'config.yaml'),
        ]
        cfg_path = None
        for p in candidates:
            if os.path.isfile(p):
                cfg_path = p
                break
        if not cfg_path:
            return None

        # Try PyYAML first
        try:
            import yaml  # type: ignore
            with open(cfg_path, 'r') as fh:
                return yaml.safe_load(fh)
        except Exception:
            pass

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
        metadata_path=args.metadata
    )

    generator.generate()
    return 0


if __name__ == '__main__':
    exit(main())
