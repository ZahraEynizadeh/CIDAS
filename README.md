# CDS overlap and in-silico Pfam change (rMATS-friendly)

Three Python scripts (stdlib only; **Python ≥ 3.6**) that:

1. **Predict coding footprint** — longest ORF per transcript (`gffread` on your GTF + genome), merged per gene; intersect alternative-splicing intervals from your event tables.
2. **Model isoform pairs + scan Pfam** — inclusion vs exclusion sequence per event → longest ORF → **`hmmscan`** against **Pfam-A**; report where domain accession sets differ (E-value cutoff configurable).
3. **Export gene lists** — one TSV per contrast for genes with ≥1 domain-composition–change event.

This is **computational screening**, not a substitute for proteomics or curated annotation.

---

## Prerequisites

| Component | Role |
|-----------|------|
| **gffread** | Transcript FASTA from GTF + genome ([StringTie / gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml)) |
| **hmmscan** | **HMMER 3** ([hmmer.org](http://hmmer.org/)) |
| **Pfam-A.hmm** | Download from [Pfam / InterPro](https://www.ebi.ac.uk/interpro/download/Pfam/) |

Clone this repo and work from its root (the directory that contains `settings.py` and `scripts/`).

---

## Data layout

Pick a **project root** directory on your machine (not part of this repo). Under it, defaults expect:

```text
(project root)/
  inputs/refs/annotation.gtf
  inputs/refs/genome.fa
  results/publication_tables/   ← event TSVs (match your settings)
  tmp/                          ← transcript FASTA cache (created)
  results/tables/               ← CDS + Pfam outputs (created)
```

Put **one tab-separated event file per contrast** in the event-table directory expected by `settings.py`. The **contrast** label usually comes from each file’s name; see `settings.py` for how that is derived.

---

## Event table format

- Tab-separated, header row, rMATS-style **event types**: `SE`, `RI`, `MXE`, `A3SS`, `A5SS`.
- Coordinates use rMATS conventions (0-based starts in `*_0base` columns, inclusive ends as in rMATS).

**Columns used**

| Purpose | Columns |
|---------|---------|
| Filtering | `is_significant` (optional; `true` / `1` / `yes`), `FDR` (≤ 0.05 if present), `abs_dPSI` (≥ 0.10 if present) |
| Identity | `ID`, `GeneID`, `Former_GeneID` (if you use a stable gene key), `event_type`, `chr`, `strand` |
| SE | `exonStart_0base`, `exonEnd` |
| RI | `riExonStart_0base`, `riExonEnd`, `upstreamES`, `upstreamEE`, `downstreamES`, `downstreamEE` |
| MXE | `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd` |
| A3SS / A5SS | `longExonStart_0base`, `longExonEnd`, `shortES`, `shortEE` |

If the table has no `contrast` column, the **contrast** is taken from the **filename** (see `settings.py`).

---

## Run

Set up paths and tools via **`settings.py`** (environment variables), then:

```bash
python3 scripts/summarize_das_predicted_cds_overlap.py
python3 scripts/das_pfam_domain_change.py --cpu 8
python3 scripts/export_pfam_domain_change_genes_by_contrast.py
```

**`das_pfam_domain_change.py` options**

- `--cpu N` — threads for `hmmscan`.
- `--ieval FLOAT` — override inclusion E-value (default is set in `settings.py`).
- `--limit N` — process only the first *N* events per contrast (testing).
- `--summarize-only` — rebuild event/summary TSVs from existing `stage1_records_and_headers.pkl` and `hmmscan_domtbl.out` (no new hmmscan).

---

## Main outputs

Exact paths depend on `settings.py`. Typical filenames:

| File | Content |
|------|---------|
| `das_overlap_predicted_cds_summary.tsv` | Per-contrast + pooled CDS overlap summary (CDS step). |
| `das_events_predicted_cds_overlap.tsv` | Per-event CDS overlap flags (CDS step). |
| `das_pfam_domain_change_events.tsv` | Per-event build/Pfam status and domain-change flag. |
| `das_pfam_domain_change_summary.tsv` | Contrast-level Pfam summary. |
| `genes_predicted_Pfam_domain_change_<contrast>.tsv` | Genes with ≥1 domain-change event (export step). |
| `all_isoform_proteins.fa`, `hmmscan_domtbl.out`, `stage1_records_and_headers.pkl` | Intermediate / reproducibility. |

---

## License

MIT — see `LICENSE`.
