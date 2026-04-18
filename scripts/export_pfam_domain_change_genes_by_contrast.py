#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
One TSV per contrast: Former_GeneID with ≥1 event where
das_pfam_domain_change_events.tsv has domain_change == true.

Reads das_pfam_domain_change_events.tsv from the configured Pfam output directory;
writes genes_predicted_Pfam_domain_change_<contrast>.tsv alongside it.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

_PKG = Path(__file__).resolve().parent.parent
if str(_PKG) not in sys.path:
    sys.path.insert(0, str(_PKG))
import settings

IN_EVENTS = settings.OUT_DIR / "das_pfam_domain_change_events.tsv"
OUT_DIR = settings.OUT_DIR


def main():
    # contrast -> Former_GeneID -> list of event rows (dict)
    by_gene = defaultdict(lambda: defaultdict(list))

    with IN_EVENTS.open(newline="") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for row in rdr:
            if str(row.get("domain_change", "")).strip().lower() != "true":
                continue
            ct = row.get("contrast", "").strip()
            gid = row.get("Former_GeneID", "").strip()
            if not ct or not gid:
                continue
            by_gene[ct][gid].append(row)

    merged_cols = [
        "contrast",
        "Former_GeneID",
        "GeneID",
        "n_events_Pfam_domain_composition_change",
        "event_types",
        "event_ids",
    ]

    for ct in sorted(by_gene.keys()):
        genes = sorted(by_gene[ct].keys())
        out_rows = []
        for gid in genes:
            evs = by_gene[ct][gid]
            mstrg = sorted({e.get("GeneID", "").strip() for e in evs if e.get("GeneID")})
            gene_id = ";".join(mstrg) if len(mstrg) > 1 else (mstrg[0] if mstrg else "")
            types = sorted({e.get("event_type", "").strip() for e in evs if e.get("event_type")})
            eids = sorted({str(e.get("event_id", "")).strip() for e in evs if e.get("event_id")})
            out_rows.append(
                {
                    "contrast": ct,
                    "Former_GeneID": gid,
                    "GeneID": gene_id,
                    "n_events_Pfam_domain_composition_change": len(evs),
                    "event_types": ";".join(types),
                    "event_ids": ";".join(eids),
                }
            )

        out_path = OUT_DIR / "genes_predicted_Pfam_domain_change_{}.tsv".format(ct)
        with out_path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=merged_cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(out_rows)

        print("Wrote {} ({} genes)".format(out_path, len(out_rows)))


if __name__ == "__main__":
    main()
