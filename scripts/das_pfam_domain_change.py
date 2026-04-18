#!/usr/bin/env python3
# Python 3.6 compatible.
"""
Per event: inclusion vs exclusion DNA → longest ORF → single hmmscan batch; flag
unordered Pfam accession-set differences (in silico; not proteomics).
"""
import argparse
import csv
import pickle
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

_PKG = Path(__file__).resolve().parent.parent
if str(_PKG) not in sys.path:
    sys.path.insert(0, str(_PKG))
import settings

GTF_PATH = settings.GTF_PATH
GENOME_FA = settings.GENOME_FA
TRANSCRIPT_FA = settings.TRANSCRIPT_FA
GFFREAD = settings.GFFREAD
HMMSCAN = settings.HMMSCAN
PUB_TABLES_DIR = settings.PUB_TABLES_DIR
OUT_DIR = settings.OUT_DIR
STAGE1_PKL = settings.STAGE1_PKL

MIN_ORF_NT = 225
MIN_PROT_AA = 50
DOM_IEVAL = settings.DOM_IEVAL
STOPS = {b"TAA", b"TAG", b"TGA"}

re_gene = re.compile(r'gene_id "([^"]+)"')
re_tx = re.compile(r'transcript_id "([^"]+)"')

CODON_TABLE = {
    b"TTT": b"F", b"TTC": b"F", b"TTA": b"L", b"TTG": b"L",
    b"TCT": b"S", b"TCC": b"S", b"TCA": b"S", b"TCG": b"S",
    b"TAT": b"Y", b"TAC": b"Y", b"TAA": b"*", b"TAG": b"*",
    b"TGT": b"C", b"TGC": b"C", b"TGA": b"*", b"TGG": b"W",
    b"CTT": b"L", b"CTC": b"L", b"CTA": b"L", b"CTG": b"L",
    b"CCT": b"P", b"CCC": b"P", b"CCA": b"P", b"CCG": b"P",
    b"CAT": b"H", b"CAC": b"H", b"CAA": b"Q", b"CAG": b"Q",
    b"CGT": b"R", b"CGC": b"R", b"CGA": b"R", b"CGG": b"R",
    b"ATT": b"I", b"ATC": b"I", b"ATA": b"I", b"ATG": b"M",
    b"ACT": b"T", b"ACC": b"T", b"ACA": b"T", b"ACG": b"T",
    b"AAT": b"N", b"AAC": b"N", b"AAA": b"K", b"AAG": b"K",
    b"AGT": b"S", b"AGC": b"S", b"AGA": b"R", b"AGG": b"R",
    b"GTT": b"V", b"GTC": b"V", b"GTA": b"V", b"GTG": b"V",
    b"GCT": b"A", b"GCC": b"A", b"GCA": b"A", b"GCG": b"A",
    b"GAT": b"D", b"GAC": b"D", b"GAA": b"E", b"GAG": b"E",
    b"GGT": b"G", b"GGC": b"G", b"GGA": b"G", b"GGG": b"G",
}


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


def revcomp(b: bytes) -> bytes:
    t = b.translate(b"".maketrans(b"ACGTacgt", b"TGCAtgca"))
    return t[::-1]


def load_genome(path: Path) -> Dict[str, bytes]:
    out: Dict[str, bytes] = {}
    name: Optional[str] = None
    chunks: List[bytes] = []
    with path.open("rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    out[name] = b"".join(chunks)
                name = line[1:].strip().split()[0].decode()
                chunks = []
            else:
                chunks.append(line.strip().upper())
        if name is not None:
            out[name] = b"".join(chunks)
    return out


def fetch_inclusive(genome: Dict[str, bytes], chrom: str, g0: int, g1: int) -> bytes:
    s = genome[chrom][g0 - 1 : g1]
    return s


def parse_gtf_exons(path: Path) -> Dict[str, List[Tuple[str, int, int, str]]]:
    strand_by_tx: Dict[str, str] = {}
    raw: Dict[str, List[Tuple[str, int, int]]] = defaultdict(list)
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            chrom, _s, _f, start_s, end_s, _sc, strand, _ph, attrs = cols
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
    return ordered


def tx_gene_map(path: Path) -> Dict[str, str]:
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


def gene_transcripts(tx_to_gene: Dict[str, str], gene: str) -> List[str]:
    return [t for t, g in tx_to_gene.items() if g == gene]


def longest_orf(seq: bytes) -> Optional[Tuple[int, int]]:
    best = (0, 0)
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
                if seq[j : j + 3] in STOPS:
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


def translate_orf(seq: bytes, orf: Tuple[int, int]) -> bytes:
    s, e = orf
    pep = bytearray()
    for i in range(s, e, 3):
        codon = seq[i : i + 3]
        if len(codon) < 3:
            break
        aa = CODON_TABLE.get(codon, b"X")
        if aa == b"*":
            break
        pep.extend(aa)
    return bytes(pep)


def load_fasta_dict(path: Path) -> Dict[str, bytes]:
    out: Dict[str, bytes] = {}
    name: Optional[str] = None
    buf: List[bytes] = []
    with path.open("rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    out[name] = b"".join(buf)
                name = line[1:].strip().split()[0].decode()
                buf = []
            else:
                buf.append(line.strip().upper())
        if name is not None:
            out[name] = b"".join(buf)
    return out


def inclusive_genomic_overlap(a0: int, a1: int, b0: int, b1: int) -> bool:
    return not (a1 < b0 or b1 < a0)


def map_genomic_interval_inside_exon(
    exons: List[Tuple[str, int, int, str]], chrom: str, g0: int, g1: int
) -> Optional[Tuple[int, int]]:
    pos = 0
    for c, gs, ge, st in exons:
        if c != chrom:
            return None
        if gs <= g0 and ge >= g1:
            if st == "+":
                t0 = pos + (g0 - gs)
                t1 = pos + (g1 - gs + 1)
            else:
                t0 = pos + (ge - g1)
                t1 = pos + (ge - g0 + 1)
            return t0, t1
        pos += ge - gs + 1
    return None


def find_tx_covering_interval(
    gene: str,
    tx_to_gene: Dict[str, str],
    exons_by_tx: Dict[str, List[Tuple[str, int, int, str]]],
    chrom: str,
    g0: int,
    g1: int,
) -> Optional[str]:
    for tid in gene_transcripts(tx_to_gene, gene):
        exons = exons_by_tx.get(tid)
        if not exons:
            continue
        if map_genomic_interval_inside_exon(exons, chrom, g0, g1) is not None:
            return tid
    return None


def se_inc_exc_seq(
    gene: str,
    row: dict,
    tx_to_gene: Dict[str, str],
    exons_by_tx: Dict[str, List[Tuple[str, int, int, str]]],
    fa: Dict[str, bytes],
) -> Tuple[Optional[bytes], Optional[bytes], str]:
    chrom = row["chr"].strip()
    c0 = fint(row.get("exonStart_0base", ""))
    c1 = fint(row.get("exonEnd", ""))
    if c0 is None or c1 is None:
        return None, None, "se_missing_coords"
    g0, g1 = c0 + 1, c1
    tid = find_tx_covering_interval(gene, tx_to_gene, exons_by_tx, chrom, g0, g1)
    if not tid or tid not in fa:
        return None, None, "se_no_transcript"
    exons = exons_by_tx[tid]
    mm = map_genomic_interval_inside_exon(exons, chrom, g0, g1)
    if mm is None:
        return None, None, "se_map_fail"
    t0, t1 = mm
    seq = fa[tid]
    inc = seq
    exc = seq[:t0] + seq[t1:]
    return inc, exc, "ok"


def ri_inc_exc_seq(row: dict, genome: Dict[str, bytes]) -> Tuple[Optional[bytes], Optional[bytes], str]:
    chrom = row["chr"].strip()
    st = row["strand"].strip()
    u0 = fint(row.get("upstreamES", ""))
    u1 = fint(row.get("upstreamEE", ""))
    d0 = fint(row.get("downstreamES", ""))
    d1 = fint(row.get("downstreamEE", ""))
    if None in (u0, u1, d0, d1):
        return None, None, "ri_missing_coords"
    u0, u1, d0, d1 = u0 + 1, u1, d0 + 1, d1
    if chrom not in genome:
        return None, None, "ri_bad_chrom"
    if st == "+":
        exc = fetch_inclusive(genome, chrom, u0, u1) + fetch_inclusive(genome, chrom, d0, d1)
        if d0 > u1 + 1:
            intron = fetch_inclusive(genome, chrom, u1 + 1, d0 - 1)
        else:
            intron = b""
        inc = fetch_inclusive(genome, chrom, u0, u1) + intron + fetch_inclusive(genome, chrom, d0, d1)
        return inc, exc, "ok"
    # minus strand: 5' exon has higher genomic coordinates (typically d block right of u)
    def blk(gs: int, ge: int) -> bytes:
        return revcomp(fetch_inclusive(genome, chrom, gs, ge))

    if d0 > u0:
        if d0 > u1 + 1:
            intron_raw = fetch_inclusive(genome, chrom, u1 + 1, d0 - 1)
        else:
            intron_raw = b""
        exc = blk(d0, d1) + blk(u0, u1)
        inc = blk(d0, d1) + revcomp(intron_raw) + blk(u0, u1)
    else:
        if u0 > d1 + 1:
            intron_raw = fetch_inclusive(genome, chrom, d1 + 1, u0 - 1)
        else:
            intron_raw = b""
        exc = blk(u0, u1) + blk(d0, d1)
        inc = blk(u0, u1) + revcomp(intron_raw) + blk(d0, d1)
    return inc, exc, "ok"


def exon_end_matches(
    tid: str,
    exons_by_tx: Dict[str, List[Tuple[str, int, int, str]]],
    chrom: str,
    strand: str,
    want_end: int,
    tol: int = 2,
) -> bool:
    exons = exons_by_tx.get(tid)
    if not exons:
        return False
    for c, gs, ge, st in exons:
        if c == chrom and st == strand and abs(ge - want_end) <= tol:
            return True
    return False


def pick_tx_by_splice_end(
    gene: str,
    tx_to_gene: Dict[str, str],
    exons_by_tx: Dict[str, List[Tuple[str, int, int, str]]],
    chrom: str,
    strand: str,
    end_coord: int,
) -> Optional[str]:
    cands = []
    for tid in gene_transcripts(tx_to_gene, gene):
        if exon_end_matches(tid, exons_by_tx, chrom, strand, end_coord):
            cands.append(tid)
    if not cands:
        return None
    return cands[0]


def ass_inc_exc_from_tx_pair(
    tid_long: str, tid_short: str, fa: Dict[str, bytes]
) -> Tuple[Optional[bytes], Optional[bytes], str]:
    if tid_long not in fa or tid_short not in fa:
        return None, None, "ass_missing_fa"
    return fa[tid_long], fa[tid_short], "ok"


def mxe_inc_exc_seq(
    gene: str,
    row: dict,
    tx_to_gene: Dict[str, str],
    exons_by_tx: Dict[str, List[Tuple[str, int, int, str]]],
    fa: Dict[str, bytes],
) -> Tuple[Optional[bytes], Optional[bytes], str]:
    chrom = row["chr"].strip()
    a0 = fint(row.get("1stExonStart_0base", ""))
    a1 = fint(row.get("1stExonEnd", ""))
    b0 = fint(row.get("2ndExonStart_0base", ""))
    b1 = fint(row.get("2ndExonEnd", ""))
    if None in (a0, a1, b0, b1):
        return None, None, "mxe_missing_coords"
    ga0, ga1 = a0 + 1, a1
    gb0, gb1 = b0 + 1, b1

    rstrand = row["strand"].strip()

    def has_a_not_b(tid: str) -> bool:
        ex = exons_by_tx.get(tid)
        if not ex:
            return False
        ha = any(
            c == chrom and rstrand == stx and inclusive_genomic_overlap(gs, ge, ga0, ga1)
            for c, gs, ge, stx in ex
        )
        hb = any(
            c == chrom and rstrand == stx and inclusive_genomic_overlap(gs, ge, gb0, gb1)
            for c, gs, ge, stx in ex
        )
        return ha and not hb

    def has_b_not_a(tid: str) -> bool:
        ex = exons_by_tx.get(tid)
        if not ex:
            return False
        ha = any(
            c == chrom and rstrand == stx and inclusive_genomic_overlap(gs, ge, ga0, ga1)
            for c, gs, ge, stx in ex
        )
        hb = any(
            c == chrom and rstrand == stx and inclusive_genomic_overlap(gs, ge, gb0, gb1)
            for c, gs, ge, stx in ex
        )
        return hb and not ha

    ta = [t for t in gene_transcripts(tx_to_gene, gene) if has_a_not_b(t)]
    tb = [t for t in gene_transcripts(tx_to_gene, gene) if has_b_not_a(t)]
    if not ta or not tb:
        return None, None, "mxe_no_tx_pair"
    return fa.get(ta[0]), fa.get(tb[0]), "ok"


def a3ss_pair(
    gene: str,
    row: dict,
    tx_to_gene: Dict[str, str],
    exons_by_tx: Dict[str, List[Tuple[str, int, int, str]]],
    fa: Dict[str, bytes],
) -> Tuple[Optional[bytes], Optional[bytes], str]:
    chrom = row["chr"].strip()
    st = row["strand"].strip()
    le = fint(row.get("longExonEnd", ""))
    se = fint(row.get("shortEE", ""))
    if le is None or se is None:
        return None, None, "a3ss_missing_coords"
    tl = pick_tx_by_splice_end(gene, tx_to_gene, exons_by_tx, chrom, st, le)
    ts = pick_tx_by_splice_end(gene, tx_to_gene, exons_by_tx, chrom, st, se)
    if tl and ts:
        return ass_inc_exc_from_tx_pair(tl, ts, fa)
    return None, None, "a3ss_no_tx_pair"


def a5ss_pair(
    gene: str,
    row: dict,
    tx_to_gene: Dict[str, str],
    exons_by_tx: Dict[str, List[Tuple[str, int, int, str]]],
    fa: Dict[str, bytes],
) -> Tuple[Optional[bytes], Optional[bytes], str]:
    chrom = row["chr"].strip()
    st = row["strand"].strip()
    l0 = fint(row.get("longExonStart_0base", ""))
    s0 = fint(row.get("shortES", ""))
    if l0 is None or s0 is None:
        return None, None, "a5ss_missing_coords"
    g_long = l0 + 1
    g_short = s0 + 1

    def tx_with_exon_start(tid: str, gstart: int, tol: int = 2) -> bool:
        for c, gs, ge, st2 in exons_by_tx.get(tid, []):
            if c == chrom and st2 == st and abs(gs - gstart) <= tol:
                return True
        return False

    t_long = [t for t in gene_transcripts(tx_to_gene, gene) if tx_with_exon_start(t, g_long)]
    t_short = [t for t in gene_transcripts(tx_to_gene, gene) if tx_with_exon_start(t, g_short)]
    if t_long and t_short:
        return ass_inc_exc_from_tx_pair(t_long[0], t_short[0], fa)
    return None, None, "a5ss_no_tx_pair"


def protein_from_dna(dna: bytes) -> Tuple[Optional[bytes], str]:
    orf = longest_orf(dna)
    if orf is None:
        return None, "no_orf"
    prot = translate_orf(dna, orf)
    if len(prot) < MIN_PROT_AA:
        return None, "short_prot"
    return prot, "ok"


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
    keys = [g, ct, et]
    for col in (
        "exonStart_0base",
        "exonEnd",
        "upstreamES",
        "upstreamEE",
        "downstreamES",
        "downstreamEE",
        "longExonStart_0base",
        "longExonEnd",
        "shortES",
        "shortEE",
        "1stExonStart_0base",
        "1stExonEnd",
        "2ndExonStart_0base",
        "2ndExonEnd",
        "riExonStart_0base",
        "riExonEnd",
    ):
        keys.append((row.get(col) or "").strip())
    return tuple(keys)


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


def parse_domtbl(path: Path, ieval_max: float) -> Dict[str, Set[str]]:
    """query name -> set of Pfam accessions (no version)."""
    hits: Dict[str, Set[str]] = defaultdict(set)
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 14:
                continue
            query = parts[3]
            acc = parts[1].split(".")[0]
            try:
                ie = float(parts[12])
            except (ValueError, IndexError):
                continue
            if ie <= ieval_max:
                hits[query].add(acc)
    return dict(hits)


def ensure_tx_fa() -> None:
    TRANSCRIPT_FA.parent.mkdir(parents=True, exist_ok=True)
    if TRANSCRIPT_FA.is_file() and TRANSCRIPT_FA.stat().st_size > 1_000_000:
        return
    subprocess.run(
        [
            str(GFFREAD),
            "-w",
            str(TRANSCRIPT_FA),
            "-g",
            str(GENOME_FA),
            str(GTF_PATH),
        ],
        check=True,
    )


def summarize_from_domtbl(
    all_records: List[dict],
    domtbl_real: Path,
    ieval: float,
) -> None:
    hit_map = parse_domtbl(domtbl_real, ieval)
    by_ev: Dict[str, Tuple[Set[str], Set[str]]] = {}
    for name in hit_map:
        parts = name.split("|")
        if len(parts) < 4:
            continue
        base = "|".join(parts[:3])
        kind = parts[3]
        if base not in by_ev:
            by_ev[base] = (set(), set())
        inc_set, exc_set = by_ev[base]
        if kind == "inc":
            inc_set |= hit_map[name]
        elif kind == "exc":
            exc_set |= hit_map[name]

    for rec in all_records:
        contrast = rec["contrast"]
        eid = rec["event_id"]
        gene = rec["GeneID"]
        base = "%s|%s|%s" % (contrast, eid, gene)
        ps = rec.get("protein_status")
        if ps == "identical_prot":
            rec["domain_change"] = "false"
            rec["pfam_inc"] = ""
            rec["pfam_exc"] = ""
            continue
        if ps != "queued":
            continue
        inc_set, exc_set = by_ev.get(base, (set(), set()))
        rec["pfam_inc"] = ";".join(sorted(inc_set))
        rec["pfam_exc"] = ";".join(sorted(exc_set))
        ch = inc_set != exc_set
        rec["domain_change"] = "true" if ch else "false"

    ev_out = OUT_DIR / "das_pfam_domain_change_events.tsv"
    with ev_out.open("w", newline="") as fh:
        fields = [
            "contrast",
            "event_id",
            "GeneID",
            "Former_GeneID",
            "event_type",
            "build_status",
            "protein_status",
            "domain_change",
            "pfam_inc",
            "pfam_exc",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for rec in all_records:
            w.writerow(rec)

    by_ct: Dict[str, List[dict]] = defaultdict(list)
    for rec in all_records:
        by_ct[rec["contrast"]].append(rec)

    summary_path = OUT_DIR / "das_pfam_domain_change_summary.tsv"
    with summary_path.open("w", newline="") as fh:
        wo = csv.writer(fh, delimiter="\t")
        wo.writerow(
            [
                "contrast",
                "n_events_input",
                "n_events_queued_hmmscan",
                "n_events_domain_change",
                "n_genes_any_das",
                "n_genes_with_hmmscan_pair",
                "n_genes_domain_change",
                "pct_genes_domain_change_of_all_das",
                "pct_genes_domain_change_of_scanned_genes",
                "note",
            ]
        )
        for ct in sorted(by_ct.keys()):
            recs = by_ct[ct]
            genes_all = {r["GeneID"] for r in recs}
            queued = [r for r in recs if r.get("protein_status") == "queued"]
            idem = [r for r in recs if r.get("protein_status") == "identical_prot"]
            gq = {r["GeneID"] for r in queued}
            changed = [r for r in queued if r.get("domain_change") == "true"]
            gc = {r["GeneID"] for r in changed}
            n_g = len(genes_all)
            pct_all = (100.0 * len(gc) / n_g) if n_g else 0.0
            pct_q = (100.0 * len(gc) / len(gq)) if gq else 0.0
            wo.writerow(
                [
                    ct,
                    len(recs),
                    len(queued),
                    len(changed),
                    n_g,
                    len(gq),
                    len(gc),
                    "%.2f" % pct_all,
                    "%.2f" % pct_q,
                    "in_silico_Pfam_iEVAL<=%s; identical_prot=%d_events"
                    % (ieval, len(idem)),
                ]
            )
    print("Wrote:", ev_out)
    print("Wrote:", summary_path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu", type=int, default=8)
    ap.add_argument("--ieval", type=float, default=DOM_IEVAL)
    ap.add_argument("--limit", type=int, default=0, help="Max events per contrast (0=all)")
    ap.add_argument(
        "--summarize-only",
        action="store_true",
        help="Read stage1 pickle + existing hmmscan_domtbl.out; write tables only",
    )
    args = ap.parse_args()

    if args.summarize_only:
        if not STAGE1_PKL.is_file():
            print("Missing %s" % STAGE1_PKL, file=sys.stderr)
            sys.exit(1)
        domtbl_real = OUT_DIR / "hmmscan_domtbl.out"
        if not domtbl_real.is_file():
            print("Missing domtbl %s" % domtbl_real, file=sys.stderr)
            sys.exit(1)
        with STAGE1_PKL.open("rb") as fh:
            payload = pickle.load(fh)
        all_records = payload["records"]
        summarize_from_domtbl(all_records, domtbl_real, args.ieval)
        return

    hm = str(HMMSCAN)
    if not Path(hm).is_file() and not shutil.which(hm):
        print("Missing hmmscan (not on PATH and not a file): %s" % hm, file=sys.stderr)
        sys.exit(1)
    pfam_hmm = settings.require_pfam_hmm()

    ensure_tx_fa()
    print("Loading genome…", flush=True)
    genome = load_genome(GENOME_FA)
    print("Loading transcript FASTA…", flush=True)
    fa = load_fasta_dict(TRANSCRIPT_FA)
    print("Parsing GTF…", flush=True)
    exons_by_tx = parse_gtf_exons(GTF_PATH)
    tx_to_gene = tx_gene_map(GTF_PATH)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    event_tables = settings.event_tables()

    all_records: List[dict] = []
    fasta_pep: List[Tuple[str, bytes]] = []

    for tab in event_tables:
        contrast = settings.contrast_label_from_table(tab)
        rows = load_deduped_rows(tab)
        if args.limit:
            rows = rows[: args.limit]
        print("Contrast %s: %d events" % (contrast, len(rows)), flush=True)
        for row in rows:
            eid = (row.get("ID") or "").strip()
            gene = (row.get("GeneID") or "").strip().replace('"', "")
            et = (row.get("event_type") or "").strip()
            inc_dna: Optional[bytes] = None
            exc_dna: Optional[bytes] = None
            st_msg = "ok"

            if et == "SE":
                inc_dna, exc_dna, st_msg = se_inc_exc_seq(
                    gene, row, tx_to_gene, exons_by_tx, fa
                )
            elif et == "RI":
                inc_dna, exc_dna, st_msg = ri_inc_exc_seq(row, genome)
            elif et == "MXE":
                inc_dna, exc_dna, st_msg = mxe_inc_exc_seq(
                    gene, row, tx_to_gene, exons_by_tx, fa
                )
            elif et == "A3SS":
                inc_dna, exc_dna, st_msg = a3ss_pair(
                    gene, row, tx_to_gene, exons_by_tx, fa
                )
            elif et == "A5SS":
                inc_dna, exc_dna, st_msg = a5ss_pair(
                    gene, row, tx_to_gene, exons_by_tx, fa
                )
            else:
                st_msg = "unknown_etype"

            qbase = "%s|%s|%s" % (contrast, eid, gene)
            p_inc: Optional[bytes] = None
            p_exc: Optional[bytes] = None
            pr_msg = ""
            if inc_dna and exc_dna:
                p_inc, m1 = protein_from_dna(inc_dna)
                p_exc, m2 = protein_from_dna(exc_dna)
                if p_inc is None:
                    pr_msg = "inc_" + m1
                elif p_exc is None:
                    pr_msg = "exc_" + m2
                elif p_inc == p_exc:
                    pr_msg = "identical_prot"
                else:
                    fasta_pep.append((qbase + "|inc", p_inc))
                    fasta_pep.append((qbase + "|exc", p_exc))
                    pr_msg = "queued"
            else:
                pr_msg = st_msg

            all_records.append(
                {
                    "contrast": contrast,
                    "event_id": eid,
                    "GeneID": gene,
                    "Former_GeneID": (row.get("Former_GeneID") or "").strip(),
                    "event_type": et,
                    "build_status": st_msg,
                    "protein_status": pr_msg,
                    "domain_change": "",
                    "pfam_inc": "",
                    "pfam_exc": "",
                }
            )

    pep_path = OUT_DIR / "all_isoform_proteins.fa"
    domtbl_real = OUT_DIR / "hmmscan_domtbl.out"
    with pep_path.open("wb") as fh:
        for name, seq in fasta_pep:
            fh.write(b">" + name.encode() + b"\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + b"\n")

    print("Running hmmscan on %d proteins…" % len(fasta_pep), flush=True)
    if not fasta_pep:
        print("No proteins to scan.", file=sys.stderr)
        sys.exit(1)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with STAGE1_PKL.open("wb") as fh:
        pickle.dump({"records": all_records}, fh)
    if domtbl_real.is_file():
        domtbl_real.unlink()
    subprocess.run(
        [
            str(HMMSCAN),
            "--domtblout",
            str(domtbl_real),
            "-E",
            "10",
            "--cpu",
            str(args.cpu),
            str(pfam_hmm),
            str(pep_path),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )
    if not domtbl_real.is_file():
        print("hmmscan domtblout missing", file=sys.stderr)
        sys.exit(1)
    summarize_from_domtbl(all_records, domtbl_real, args.ieval)


if __name__ == "__main__":
    main()
