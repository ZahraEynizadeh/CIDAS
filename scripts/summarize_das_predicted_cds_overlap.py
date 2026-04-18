#!/usr/bin/env python3
# Python 3.6 compatible (no postponed evaluation of annotations).
"""
rMATS events vs longest-ORF predicted CDS (gffread transcripts → ORFs → genomic intervals).

Emits per-contrast summaries plus a pooled deduped row. Not the same as domain-level
disruption; use hmmscan / InterPro on isoform proteins for that.
"""
import csv
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

_PKG = Path(__file__).resolve().parent.parent
if str(_PKG) not in sys.path:
    sys.path.insert(0, str(_PKG))
import settings

GFFREAD = settings.GFFREAD
GTF_PATH = settings.GTF_PATH
GENOME_FA = settings.GENOME_FA
OUT_SUMMARY = settings.OUT_SUMMARY
OUT_EVENTS = settings.OUT_EVENTS
NOTE_SHORT = "predicted_CDS=longest_ORF_merged; not_Pfam_domain_change"
TRANSCRIPT_FA = settings.TRANSCRIPT_FA

MIN_ORF_NT = 225  # 75 codons

STOPS = {b"TAA", b"TAG", b"TGA"}

re_gene = re.compile(r'gene_id "([^"]+)"')
re_tx = re.compile(r'transcript_id "([^"]+)"')


def parse_gtf_exons(path: Path) -> Tuple[Dict[str, str], Dict[str, List[Tuple[str, int, int, str]]]]:
    """transcript_id -> strand; transcript_id -> exons (chr, start, end inclusive 1-based, strand)."""
    strand_by_tx: Dict[str, str] = {}
    raw: Dict[str, List[Tuple[str, int, int]]] = defaultdict(list)
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            chrom, _src, _feat, start_s, end_s, _sc, strand, _ph, attrs = cols
            mg = re_gene.search(attrs)
            mt = re_tx.search(attrs)
            if not mt:
                continue
            tx = mt.group(1)
            strand_by_tx[tx] = strand
            gs, ge = int(start_s), int(end_s)
            raw[tx].append((chrom, gs, ge))
    ordered: Dict[str, List[Tuple[str, int, int, str]]] = {}
    for tx, exons in raw.items():
        st = strand_by_tx.get(tx, "+")
        if st == "+":
            exons.sort(key=lambda x: x[1])
        else:
            exons.sort(key=lambda x: x[1], reverse=True)
        ordered[tx] = [(c, a, b, st) for c, a, b in exons]
    return strand_by_tx, ordered


def tx_gene_map(path: Path) -> Dict[str, str]:
    """transcript_id -> gene_id (MSTRG)."""
    gmap: Dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            attrs = cols[8]
            mg = re_gene.search(attrs)
            mt = re_tx.search(attrs)
            if mg and mt:
                gmap[mt.group(1)] = mg.group(1)
    return gmap


def longest_orf(seq: bytes) -> Optional[Tuple[int, int]]:
    """Return (start, end) 0-based half-open on transcript for longest ORF, or None."""
    best: Tuple[int, int] = (0, 0)
    slen = len(seq)
    for frame in (0, 1, 2):
        i = frame
        while i + 3 <= slen:
            if seq[i : i + 3] != b"ATG":
                i += 3
                continue
            start = i
            j = i + 3
            while j + 3 <= slen:
                codon = seq[j : j + 3]
                if codon in STOPS:
                    end = j + 3
                    if end - start > best[1] - best[0]:
                        best = (start, end)
                    break
                j += 3
            else:
                end = slen
                if end - start > best[1] - best[0]:
                    best = (start, end)
            i += 3
    if best[1] - best[0] < MIN_ORF_NT:
        return None
    return best


def orf_to_genomic_ivs(
    exons: Sequence[Tuple[str, int, int, str]], t0: int, t1: int
) -> List[Tuple[str, int, int]]:
    """Map transcript half-open [t0,t1) to list of (chr, g_start, g_end) inclusive 1-based."""
    out: List[Tuple[str, int, int]] = []
    pos = 0
    for chrom, gs, ge, st in exons:
        exlen = ge - gs + 1
        exon_t0 = pos
        exon_t1 = pos + exlen
        lo = max(t0, exon_t0)
        hi = min(t1, exon_t1)
        if lo < hi:
            if st == "+":
                g_start = gs + (lo - exon_t0)
                g_end = gs + (hi - exon_t0) - 1
            else:
                g_start = ge - (hi - exon_t0) + 1
                g_end = ge - (lo - exon_t0)
            out.append((chrom, g_start, g_end))
        pos += exlen
    return out


def merge_ivs(ivs: List[Tuple[str, int, int]]) -> List[Tuple[str, int, int]]:
    by_chr: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for c, a, b in ivs:
        by_chr[c].append((a, b))
    merged: List[Tuple[str, int, int]] = []
    for c, lst in by_chr.items():
        lst.sort()
        cur_s, cur_e = lst[0]
        for s, e in lst[1:]:
            if s <= cur_e + 1:
                cur_e = max(cur_e, e)
            else:
                merged.append((c, cur_s, cur_e))
                cur_s, cur_e = s, e
        merged.append((c, cur_s, cur_e))
    return merged


def inclusive_overlap(
    a0: int, a1: int, b0: int, b1: int
) -> bool:
    """Inclusive 1-based intervals [a0,a1] and [b0,b1]."""
    return not (a1 < b0 or b1 < a0)


def iv_overlap_chr(
    chr_a: str, ivs_a: List[Tuple[int, int]], chr_b: str, b0: int, b1: int
) -> bool:
    if chr_a != chr_b:
        return False
    for s, e in ivs_a:
        if inclusive_overlap(s, e, b0, b1):
            return True
    return False


def fnum(x: str) -> Optional[float]:
    x = (x or "").strip()
    if not x or x == ".":
        return None
    try:
        return float(x)
    except ValueError:
        return None


def fint(x: str) -> Optional[int]:
    v = fnum(x)
    if v is None:
        return None
    return int(v)


def event_genomic_intervals(row: dict) -> List[Tuple[str, int, int]]:
    """Return inclusive 1-based intervals (chr, start, end) for the variable region."""
    chrom = row.get("chr", "").strip()
    et = row.get("event_type", "").strip()
    ivs: List[Tuple[str, int, int]] = []

    def add_pair(s0: Optional[int], e1: Optional[int]) -> None:
        if s0 is None or e1 is None:
            return
        # rMATS: 0-based left, 1-based right (inclusive last base)
        g0 = s0 + 1
        g1 = e1
        if g0 <= g1:
            ivs.append((chrom, g0, g1))

    if et == "SE":
        add_pair(fint(row.get("exonStart_0base", "")), fint(row.get("exonEnd", "")))
    elif et == "RI":
        add_pair(fint(row.get("riExonStart_0base", "")), fint(row.get("riExonEnd", "")))
    elif et in ("A3SS", "A5SS"):
        add_pair(fint(row.get("longExonStart_0base", "")), fint(row.get("longExonEnd", "")))
        add_pair(fint(row.get("shortES", "")), fint(row.get("shortEE", "")))
    elif et == "MXE":
        add_pair(fint(row.get("1stExonStart_0base", "")), fint(row.get("1stExonEnd", "")))
        add_pair(fint(row.get("2ndExonStart_0base", "")), fint(row.get("2ndExonEnd", "")))
    return ivs


def row_passes_filters(row: dict) -> bool:
    sig = (row.get("is_significant") or "").strip().lower()
    if sig and sig not in ("true", "1", "yes"):
        return False
    fdr = fnum(row.get("FDR", ""))
    adpsi = fnum(row.get("abs_dPSI", ""))
    if fdr is not None and fdr > 0.05:
        return False
    if adpsi is not None and adpsi < 0.10:
        return False
    return True


def row_dedup_key(row: dict) -> Tuple[str, ...]:
    g = (row.get("GeneID") or "").strip().replace('"', "")
    ct = (row.get("chr") or "").strip()
    et = (row.get("event_type") or "").strip()
    ev_ivs = event_genomic_intervals(row)
    parts = [g, ct, et] + ["%s-%s" % (a, b) for _, a, b in ev_ivs]
    return tuple(parts)


def load_deduped_rows(tab: Path) -> List[dict]:
    seen: Set[Tuple[str, ...]] = set()
    out: List[dict] = []
    with tab.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not row_passes_filters(row):
                continue
            key = row_dedup_key(row)
            if key in seen:
                continue
            seen.add(key)
            out.append(row)
    return out


def evaluate_cds_overlap(
    rows: List[dict],
    gene_merged: Dict[str, List[Tuple[str, int, int]]],
    contrast_label: str = "",
) -> Tuple[Set[str], Set[str], int, List[dict]]:
    genes_all: Set[str] = set()
    genes_hit: Set[str] = set()
    event_hits = 0
    out_rows: List[dict] = []
    for row in rows:
        gid = (row.get("GeneID") or "").strip().replace('"', "")
        genes_all.add(gid)
        cds = gene_merged.get(gid, [])
        ev_ivs = event_genomic_intervals(row)
        hit = False
        for chrom, a0, a1 in ev_ivs:
            for c2, s, e in cds:
                if iv_overlap_chr(c2, [(s, e)], chrom, a0, a1):
                    hit = True
                    break
            if hit:
                break
        if hit:
            genes_hit.add(gid)
            event_hits += 1
        base = {
            k: row.get(k, "")
            for k in (
                "GeneID",
                "Former_GeneID",
                "event_type",
                "chr",
                "strand",
            )
        }
        ct = contrast_label if contrast_label else (row.get("contrast") or "").strip()
        out_rows.append(
            {
                **base,
                "contrast": ct,
                "overlaps_predicted_cds": "true" if hit else "false",
            }
        )
    return genes_all, genes_hit, event_hits, out_rows


def load_pooled_cross_contrast_deduped(event_tables: List[Path]) -> List[dict]:
    """Same genomic event in multiple contrasts counts once (for study-wide totals)."""
    seen: Set[Tuple[str, ...]] = set()
    out: List[dict] = []
    for tab in event_tables:
        with tab.open(newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if not row_passes_filters(row):
                    continue
                key = row_dedup_key(row)
                if key in seen:
                    continue
                seen.add(key)
                out.append(row)
    return out


def ensure_transcript_fa() -> None:
    TRANSCRIPT_FA.parent.mkdir(parents=True, exist_ok=True)
    if TRANSCRIPT_FA.is_file() and TRANSCRIPT_FA.stat().st_size > 1_000_000:
        return
    cmd = [
        str(GFFREAD),
        "-w",
        str(TRANSCRIPT_FA),
        "-g",
        str(GENOME_FA),
        str(GTF_PATH),
    ]
    subprocess.run(cmd, check=True)


def read_fasta_records(path: Path) -> Iterable[Tuple[str, bytes]]:
    name: Optional[str] = None
    chunks: List[bytes] = []
    with path.open("rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(chunks)
                name = line[1:].strip().split()[0].decode()
                chunks = []
            else:
                chunks.append(line.strip().upper())
        if name is not None:
            yield name, b"".join(chunks)


def main() -> None:
    event_tables = settings.event_tables()
    if not event_tables:
        print("No publication event tables found.", file=sys.stderr)
        sys.exit(1)
    gff = str(GFFREAD)
    if not Path(gff).is_file() and not shutil.which(gff):
        print("Missing gffread (not on PATH and not a file): %s" % gff, file=sys.stderr)
        sys.exit(1)

    print("Parsing GTF exon structure…", flush=True)
    _strands, exons_by_tx = parse_gtf_exons(GTF_PATH)
    tx_to_gene = tx_gene_map(GTF_PATH)

    print("Extracting transcript sequences with gffread (cached)…", flush=True)
    ensure_transcript_fa()

    print("Computing longest ORF per transcript and merging per gene…", flush=True)
    gene_cds_ivs: Dict[str, List[Tuple[str, int, int]]] = defaultdict(list)
    n_tx = 0
    for txid, seq in read_fasta_records(TRANSCRIPT_FA):
        n_tx += 1
        if n_tx % 10000 == 0:
            print(f"  …{n_tx} transcripts", flush=True)
        gid = tx_to_gene.get(txid)
        if not gid:
            continue
        exons = exons_by_tx.get(txid)
        if not exons:
            continue
        orf = longest_orf(seq)
        if not orf:
            continue
        t0, t1 = orf
        for triple in orf_to_genomic_ivs(exons, t0, t1):
            gene_cds_ivs[gid].append(triple)

    gene_merged: Dict[str, List[Tuple[str, int, int]]] = {
        g: merge_ivs(ivs) for g, ivs in gene_cds_ivs.items()
    }

    print("Loading rMATS event tables (per contrast + pooled)…", flush=True)
    summary_table: List[List[str]] = []
    all_event_out_rows: List[dict] = []

    for tab in event_tables:
        ctab = settings.contrast_label_from_table(tab)
        rows_ct = load_deduped_rows(tab)
        ga, gh, _eh, ore = evaluate_cds_overlap(rows_ct, gene_merged, ctab)
        ng, nh = len(ga), len(gh)
        pct = (100.0 * nh / ng) if ng else 0.0
        summary_table.append(
            [
                ctab,
                str(len(rows_ct)),
                str(ng),
                str(nh),
                "%.2f" % pct,
                NOTE_SHORT,
            ]
        )
        all_event_out_rows.extend(ore)
        print(
            "  %s: %d events, %d genes, %d overlap predicted CDS (%.1f%%)"
            % (ctab, len(rows_ct), ng, nh, pct),
            flush=True,
        )

    pooled = load_pooled_cross_contrast_deduped(event_tables)
    ga_p, gh_p, _eh_p, _ = evaluate_cds_overlap(pooled, gene_merged, "")
    ngp, nhp = len(ga_p), len(gh_p)
    pct_p = (100.0 * nhp / ngp) if ngp else 0.0
    summary_table.append(
        [
            "pooled_cross_contrast_event_deduped",
            str(len(pooled)),
            str(ngp),
            str(nhp),
            "%.2f" % pct_p,
            NOTE_SHORT + "; same_genomic_event_counted_once",
        ]
    )
    print(
        "  pooled (cross-contrast deduped): %d events, %d genes, %d overlap (%.1f%%)"
        % (len(pooled), ngp, nhp, pct_p),
        flush=True,
    )

    OUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with OUT_SUMMARY.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "contrast",
                "n_significant_events_deduped",
                "n_unique_genes",
                "genes_with_event_overlapping_predicted_cds",
                "percent_genes",
                "note",
            ]
        )
        for row in summary_table:
            w.writerow(row)

    with OUT_EVENTS.open("w", newline="") as fh:
        fields = list(all_event_out_rows[0].keys()) if all_event_out_rows else []
        if fields:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(all_event_out_rows)

    print()
    print("Wrote:", OUT_SUMMARY)
    print("Wrote:", OUT_EVENTS)


if __name__ == "__main__":
    main()
