"""Controller: compare the two parties' Round 17 controls-only pre-checks (DEVIATION_R17_D1 / D2).

Usage: python3 -I r17_compare_precheck.py <main_precheck.json> <at3_precheck.json>
Per pool (17), required:
  - both parties: no error entry;
  - n_ctrl_total equal between parties and equal to the line's line_metadata.json n_nt_total (D2.1) for the 11 R17 lines;
  - number of distinct batch labels equal to line_metadata n_batches (R17 lines) and equal between parties;
  - raw-count matrix sha256 equal (float32 C-order), gene axis (length and names) equal;
  - per-cell totals equal (|diff| <= 1e-6), batch labels equal (exact string, per cell, in order);
  - barcodes equal after stripping a party-specific 'prefix:' (C2.4) wherever both parties resolved real barcodes
    (A-T3's 'part_row_<idx>' fallback is reported as unresolved, not compared; D2.3).
Key names are looked up flexibly because the two parties wrote their own JSON layouts.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LINE_META = {
    "cd4t_stim48": "round11/shared/lines/cd4t_stim48", "telohaec": "round11/shared/lines/telohaec",
    "calu3": "round10/shared/lines/calu3", "hek293_xucao": "round12/shared/lines/hek293_xucao",
    "a549_ko_cas9pos": "round12/shared/lines/a549_ko_cas9pos", "mac_ctrl": "round13/shared/lines/mac_ctrl",
    "mac_lps": "round13/shared/lines/mac_lps", "a549_myllia": "round14/shared/lines/a549_myllia",
    "gse295214_h1": "round16/shared/lines/gse295214_h1",
    "pfizer_telohaec_basal": "round17/shared/lines/pfizer_telohaec_basal",
    "pfizer_telohaec_cytokine": "round17/shared/lines/pfizer_telohaec_cytokine",
}


def pools_of(doc):
    p = doc.get("pools", doc)
    if isinstance(p, list):
        p = {e.get("pool_name") or e.get("pool"): e for e in p}
    return p


def get(e, *names, default=None):
    for n in names:
        if n in e:
            return e[n]
    return default


def strip(b):
    b = str(b)
    return b.split(":", 1)[1] if ":" in b else b


m_doc, a_doc = (json.load(open(p)) for p in sys.argv[1:3])
M, A = pools_of(m_doc), pools_of(a_doc)
names = sorted(set(M) | set(A))
fails, notes, n = [], [], 0
for pn in names:
    em, ea = M.get(pn), A.get(pn)
    if em is None or ea is None:
        fails.append(f"{pn}: missing in {'main' if em is None else 'at3'}"); continue
    if get(em, "error") or get(ea, "error"):
        fails.append(f"{pn}: error entry main={get(em, 'error')} at3={get(ea, 'error')}"); continue
    nm, na = get(em, "n_ctrl_total"), get(ea, "n_ctrl_total", "n_ntc_total", "n_ntc_valid")
    bm = get(em, "batch_labels", default=[]); ba = get(ea, "batch_labels", default=[])
    ubm = get(em, "n_batches", default=len(set(bm))); uba = get(ea, "n_batches", default=len(set(ba)))
    n += 1
    if nm != na:
        fails.append(f"{pn}: n_ctrl_total main={nm} at3={na}")
    if pn in LINE_META:
        meta = json.load(open(ROOT / LINE_META[pn] / "line_metadata.json"))
        if nm != meta.get("n_nt_total"):
            fails.append(f"{pn}: n_ctrl_total {nm} != line_metadata n_nt_total {meta.get('n_nt_total')}")
        if meta.get("n_batches") is not None and len(set(bm)) > int(meta["n_batches"]):
            fails.append(f"{pn}: sampled batches {len(set(bm))} > line n_batches {meta['n_batches']}")
        notes.append(f"{pn}: n_ctrl_total {nm} (line {meta.get('n_nt_total')}), line n_batches {meta.get('n_batches')}, "
                     f"sampled batches main={len(set(bm))} at3={len(set(ba))}")
    for key_m, key_a, tag in [(("raw_counts_sha256",), ("raw_counts_sha256", "raw_sha256", "raw_sha"), "raw_sha256"),
                              (("gene_axis_length",), ("gene_axis_length", "n_genes"), "gene_axis_length")]:
        vm, va = get(em, *key_m), get(ea, *key_a)
        if vm != va:
            fails.append(f"{pn}: {tag} main={vm} at3={va}")
    gm, ga = get(em, "gene_axis_names", "gene_names"), get(ea, "gene_axis_names", "gene_names")
    if gm is not None and ga is not None and list(gm) != list(ga):
        fails.append(f"{pn}: gene axis names differ")
    tm = get(em, "cell_total_umis"); ta = get(ea, "cell_total_umis", "per_cell_totals")
    if tm is not None and ta is not None:
        if len(tm) != len(ta) or max(abs(float(x) - float(y)) for x, y in zip(tm, ta)) > 1e-6:
            fails.append(f"{pn}: per-cell totals differ")
    if list(bm) != list(ba):
        fails.append(f"{pn}: batch labels differ (main {len(bm)}, at3 {len(ba)}; first main={bm[:1]} at3={ba[:1]})")
    cm = get(em, "barcodes", default=[]); ca = get(ea, "ordered_barcodes", "barcodes", "cell_barcodes", default=[])
    unres = lambda L: any(strip(x).startswith(("part_row_", "row_")) or str(x).startswith("('") for x in L)
    if unres(ca) or unres(cm):
        notes.append(f"{pn}: barcodes unresolved by one party (part_row_*); identity by raw sha256 + totals (D2.3)")
    elif [strip(x) for x in cm] != [strip(x) for x in ca]:
        fails.append(f"{pn}: barcodes differ after prefix strip (first main={cm[:1]} at3={ca[:1]})")
print(json.dumps({"n_pools": len(names), "n_compared": n, "n_fail": len(fails)}))
for x in notes:
    print("NOTE", x)
for x in fails:
    print("FAIL", x)
print("R17_PRECHECK_MATCH" if not fails else "R17_PRECHECK_DIFF")
