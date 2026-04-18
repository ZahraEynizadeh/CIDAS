#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Path and tool configuration for the CDS → Pfam pipeline.

Configure via environment variables (see README.md). Optional
PIPELINE_EVENT_TABLE_STRIP_PREFIX controls how contrast names are derived from
event table filenames.
"""
import os
import sys
from pathlib import Path


def _require_root():
    raw = os.environ.get("PIPELINE_PROJECT_ROOT", "").strip()
    if not raw:
        print(
            "ERROR: PIPELINE_PROJECT_ROOT is not set.\n"
            "  Set it to the directory that contains your inputs and results, then re-run.\n"
            "See README.md.",
            file=sys.stderr,
        )
        sys.exit(2)
    p = Path(raw).expanduser().resolve()
    if not p.is_dir():
        print("ERROR: PIPELINE_PROJECT_ROOT is not a directory: %s" % p, file=sys.stderr)
        sys.exit(2)
    return p


PROJECT_ROOT = _require_root()


def _under(key, relative_parts):
    v = os.environ.get(key, "").strip()
    if v:
        return Path(v).expanduser().resolve()
    return PROJECT_ROOT.joinpath(*relative_parts)


def _tool(env_key, default_cmd):
    return (os.environ.get(env_key, default_cmd) or default_cmd).strip()


# --- Reference inputs ---
GTF_PATH = _under("PIPELINE_GTF", ("inputs", "refs", "annotation.gtf"))
GENOME_FA = _under("PIPELINE_GENOME_FA", ("inputs", "refs", "genome.fa"))

# --- Significant event TSV directory ---
EVENT_TABLES_DIR = _under(
    "PIPELINE_EVENT_TABLES_DIR", ("results", "publication_tables")
)
EVENT_TABLE_GLOB = os.environ.get(
    "PIPELINE_EVENT_TABLE_GLOB", "unique_only_significant_events_*.tsv"
)

TRANSCRIPT_FA = _under("PIPELINE_TRANSCRIPT_FA", ("tmp", "transcripts_unified.fa"))

OUT_SUMMARY = _under(
    "PIPELINE_CDS_SUMMARY_OUT",
    ("results", "tables", "das_overlap_predicted_cds_summary.tsv"),
)
OUT_EVENTS = _under(
    "PIPELINE_CDS_EVENTS_OUT",
    ("results", "tables", "das_events_predicted_cds_overlap.tsv"),
)

PUB_TABLES_DIR = _under(
    "PIPELINE_PUB_TABLES_DIR", ("results", "publication_tables")
)
OUT_DIR = _under(
    "PIPELINE_PFAM_OUT_DIR", ("results", "tables", "pfam_domain_change")
)
STAGE1_PKL = OUT_DIR / "stage1_records_and_headers.pkl"

GFFREAD = _tool("PIPELINE_GFFREAD", "gffread")
HMMSCAN = _tool("PIPELINE_HMMSCAN", "hmmscan")

_pfam_env = os.environ.get("PIPELINE_PFAM_HMM", "").strip()
PFAM_HMM = Path(_pfam_env).expanduser().resolve() if _pfam_env else None

DOM_IEVAL = float(os.environ.get("PIPELINE_PFAM_IEVAL", "1e-5"))

ROOT = PROJECT_ROOT

# If event tables are named like PREFIX_conditionA_vs_conditionB.tsv, set this prefix so
# the contrast label becomes conditionA_vs_conditionB. If unset, uses PREFIX below; set to
# empty string in the environment to use the full file stem as the contrast label.
_raw_strip = os.environ.get("PIPELINE_EVENT_TABLE_STRIP_PREFIX")
if _raw_strip is None:
    EVENT_TABLE_STRIP_PREFIX = "unique_only_significant_events_"
else:
    EVENT_TABLE_STRIP_PREFIX = _raw_strip


def contrast_label_from_table(path: Path) -> str:
    """Contrast name from event table filename (one table per contrast recommended)."""
    stem = path.stem
    p = EVENT_TABLE_STRIP_PREFIX
    if p and stem.startswith(p):
        return stem[len(p) :]
    return stem


def event_tables():
    if not EVENT_TABLES_DIR.is_dir():
        print("ERROR: EVENT_TABLES_DIR missing: %s" % EVENT_TABLES_DIR, file=sys.stderr)
        sys.exit(2)
    paths = sorted(EVENT_TABLES_DIR.glob(EVENT_TABLE_GLOB))
    if not paths:
        print(
            "ERROR: No files matching %s under %s"
            % (EVENT_TABLE_GLOB, EVENT_TABLES_DIR),
            file=sys.stderr,
        )
        sys.exit(2)
    return paths


def require_pfam_hmm():
    if PFAM_HMM is None or not PFAM_HMM.is_file():
        print(
            "ERROR: Set PIPELINE_PFAM_HMM to your Pfam-A.hmm file (see README).",
            file=sys.stderr,
        )
        sys.exit(2)
    return PFAM_HMM
