#!/usr/bin/env python3
"""
Pilot A-2 (Nature Genetics / Nature Medicine):
Cross-Ancestry "Reclassification Debt" Across 8 ACMG SF v3.2 Actionable Genes
Using Already-Cached Full-Gene gnomAD v4.1 (807k Individuals) Reverse-Linked to ClinVar.

Built-in Frequency-Matched Control Design:
  - Non-European-enriched coding/splice alleles: `non_nfe_max_af >= 1e-4` AND `non_nfe_max_af >= 5 * nfe_af`
  - European-enriched coding/splice alleles:     `nfe_af >= 1e-4` AND `nfe_af >= 5 * non_nfe_max_af`
Measures persistent ClinVar Pathogenic / Likely pathogenic (`>= 2-star` and all stars),
Conflicting interpretations, and VUS debt across both ancestry-enriched arms!

Uses ONLY official skill CLI wrappers (`clinvar_api.py`, cached `gnomAD` outputs).
"""

from __future__ import annotations

import json
import re
import subprocess
from collections import Counter
from pathlib import Path

from scipy import stats

PORT = Path("/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026")
DATA_DIR = PORT / "data" / "pilot_a2"
DATA_DIR.mkdir(parents=True, exist_ok=True)
UV = "/Users/mac/.local/bin/uv"
CV_CLI = "/Users/mac/.gemini/config/plugins/genomics/skills/clinvar_database/scripts/clinvar_api.py"

GENES = ["MYH7", "SCN5A", "KCNQ1", "LMNA", "MYBPC3", "LDLR", "BRCA1", "BRCA2"]

CODING_SPLICE_TERMS = (
    "missense", "nonsense", "stop_gained", "frameshift", "splice_acceptor",
    "splice_donor", "splice_region", "inframe", "start_lost", "stop_lost",
)

HIGH_CONF_STATUSES = {
    "reviewed by expert panel",
    "practice guideline",
    "criteria provided, multiple submitters, no conflicts",
}


def run_json(cmd: list[str], out_f: Path):
    if not out_f.exists():
        p = subprocess.run(cmd, cwd=str(PORT), capture_output=True, text=True)
        if p.returncode != 0:
            raise RuntimeError(f"Command failed ({out_f.name}): {p.stderr[-400:]}")
    return json.loads(out_f.read_text(encoding="utf-8"))


def parse_pop_afs(variant: dict) -> dict:
    ac_by_pop: dict[str, int] = {}
    an_by_pop: dict[str, int] = {}
    tot_ac = 0
    tot_an = 0
    for src in ("exome", "genome"):
        block = variant.get(src)
        if not block:
            continue
        tot_ac += int(block.get("ac") or 0)
        tot_an += int(block.get("an") or 0)
        for pop in block.get("populations") or []:
            pid = str(pop.get("id", ""))
            if pid not in ("afr", "amr", "eas", "mid", "nfe", "sas"):
                continue
            ac_by_pop[pid] = ac_by_pop.get(pid, 0) + int(pop.get("ac") or 0)
            an_by_pop[pid] = an_by_pop.get(pid, 0) + int(pop.get("an") or 0)
    afs = {k: ac_by_pop[k] / an_by_pop[k] for k in ac_by_pop if an_by_pop.get(k, 0) >= 2000}
    nfe_af = afs.get("nfe", 0.0)
    non_nfe = {k: v for k, v in afs.items() if k != "nfe"}
    non_nfe_pop = max(non_nfe, key=lambda k: non_nfe[k]) if non_nfe else None
    non_nfe_af = non_nfe[non_nfe_pop] if non_nfe_pop else 0.0
    return {
        "global_af": (tot_ac / tot_an) if tot_an > 0 else 0.0,
        "nfe_af": nfe_af,
        "non_nfe_pop": non_nfe_pop,
        "non_nfe_af": non_nfe_af,
        "pop_afs": {k: round(v, 6) for k, v in afs.items()},
    }


def classify_clinsig(csig: str) -> str:
    c = csig.lower().strip()
    if "conflict" in c:
        return "Conflicting"
    if ("pathogenic" in c) and ("benign" not in c) and ("uncertain" not in c):
        return "P_LP"
    if "uncertain" in c:
        return "VUS"
    if "benign" in c:
        return "B_LB"
    return "Other"


def query_clinvar_for_rsids(gene: str, arm: str, rsids: list[str]) -> list[dict]:
    """Query ClinVar Entrez for matching rsIDs in chunks of 40 and fetch ClinVar summaries."""
    if not rsids:
        return []
    found_vids = []
    for i in range(0, len(rsids), 40):
        chunk = rsids[i:i + 40]
        q_rs = " OR ".join(f"{r}[rs]" for r in chunk)
        q = f"({gene}[gene]) AND ({q_rs})"
        sf = DATA_DIR / f"cv_rev_{gene}_{arm}_s{i}.json"
        res = run_json([UV, "run", CV_CLI, "search", "--query", q, "--retmax", "200", "--output", str(sf)], sf)
        found_vids.extend(res.get("variant_ids", []))
    found_vids = sorted(set(found_vids))
    summaries = []
    for i in range(0, len(found_vids), 100):
        chunk = found_vids[i:i + 100]
        bf = DATA_DIR / f"cv_rev_{gene}_{arm}_sum{i}.json"
        s_res = run_json([UV, "run", CV_CLI, "summary", "--variant_ids", *chunk, "--output", str(bf)], bf)
        if isinstance(s_res, list):
            summaries.extend(s_res)
    return summaries


def main() -> None:
    per_gene = []
    totals = {
        "non_nfe_enriched": Counter(),
        "nfe_enriched": Counter(),
    }
    flagged_plp_examples = []
    flagged_conflict_examples = []
    driver_pops = Counter()

    for gene in GENES:
        c_raw = json.loads((DATA_DIR / f"const_{gene}.json").read_text(encoding="utf-8"))
        g_blk = c_raw.get("data", {}).get("gene") or {}
        loeuf = (g_blk.get("gnomad_constraint") or {}).get("oe_lof_upper")

        gn_raw = json.loads((DATA_DIR / f"gnomad_{gene}.json").read_text(encoding="utf-8"))
        gn_vars = gn_raw if isinstance(gn_raw, list) else (gn_raw.get("variants") or gn_raw.get("data", {}).get("gene", {}).get("variants") or [])

        non_nfe_cands = []
        nfe_cands = []
        for v in gn_vars:
            csq = str(v.get("consequence") or "").lower()
            if not any(t in csq for t in CODING_SPLICE_TERMS):
                continue
            rs_list = [r for r in (v.get("rsids") or []) if str(r).startswith("rs")]
            if not rs_list:
                continue
            afs = parse_pop_afs(v)
            # Arm 1: Non-European-enriched coding/splice allele (AF >= 1e-4 & >= 5x NFE)
            if afs["non_nfe_af"] >= 1e-4 and afs["non_nfe_af"] >= 5.0 * max(afs["nfe_af"], 1e-6):
                non_nfe_cands.append({"rsid": rs_list[0], "vid": v.get("variant_id"), "csq": csq, **afs})
                driver_pops[afs["non_nfe_pop"]] += 1
            # Arm 2: European-enriched coding/splice allele (AF >= 1e-4 & >= 5x non-NFE)
            elif afs["nfe_af"] >= 1e-4 and afs["nfe_af"] >= 5.0 * max(afs["non_nfe_af"], 1e-6):
                nfe_cands.append({"rsid": rs_list[0], "vid": v.get("variant_id"), "csq": csq, **afs})

        # Query ClinVar for BOTH ancestry-enriched arms (matched AF >= 1e-4, >= 5x enrichment!)
        cv_non_nfe = query_clinvar_for_rsids(gene, "nonnfe", [x["rsid"] for x in non_nfe_cands])
        cv_nfe = query_clinvar_for_rsids(gene, "nfe", [x["rsid"] for x in nfe_cands])

        def tally(cands, sums, arm_name):
            cnt = Counter()
            cnt["gnomad_coding_enriched"] = len(cands)
            cnt["clinvar_records"] = len(sums)
            for s in sums:
                cls = classify_clinsig(str(s.get("clinical_significance") or ""))
                rev = str(s.get("review_status") or "").lower().strip()
                is_2star = rev in HIGH_CONF_STATUSES
                cnt[cls] += 1
                if cls == "P_LP" and is_2star:
                    cnt["P_LP_2star"] += 1
                if arm_name == "non_nfe_enriched" and cls == "P_LP" and len(flagged_plp_examples) < 30:
                    flagged_plp_examples.append({
                        "gene": gene, "loeuf": round(loeuf, 3) if loeuf else None,
                        "clinvar_id": s.get("variant_id"), "title": s.get("title"),
                        "clinsig": s.get("clinical_significance"), "review_status": s.get("review_status"),
                        "last_evaluated": s.get("last_evaluated"),
                    })
                if arm_name == "non_nfe_enriched" and cls == "Conflicting" and len(flagged_conflict_examples) < 20:
                    flagged_conflict_examples.append({
                        "gene": gene, "loeuf": round(loeuf, 3) if loeuf else None,
                        "clinvar_id": s.get("variant_id"), "title": s.get("title"),
                        "clinsig": s.get("clinical_significance"), "last_evaluated": s.get("last_evaluated"),
                    })
            totals[arm_name].update(cnt)
            return dict(cnt)

        t_non = tally(non_nfe_cands, cv_non_nfe, "non_nfe_enriched")
        t_nfe = tally(nfe_cands, cv_nfe, "nfe_enriched")

        per_gene.append({
            "gene": gene,
            "loeuf": round(loeuf, 4) if loeuf is not None else None,
            "gnomad_total_vars": len(gn_vars),
            "non_nfe_enriched": t_non,
            "nfe_enriched": t_nfe,
        })

    # Fisher exact tests comparing Non-NFE-enriched vs NFE-enriched ClinVar ambiguity/pathogenicity debt
    tn = totals["non_nfe_enriched"]
    te = totals["nfe_enriched"]

    def fisher_rate(num_a, den_a, num_b, den_b):
        if den_a == 0 or den_b == 0:
            return None
        tbl = [[num_a, max(den_a - num_a, 0)], [num_b, max(den_b - num_b, 0)]]
        or_v, pv = stats.fisher_exact(tbl)
        return {
            "non_nfe_rate": round(num_a / den_a, 4),
            "nfe_rate": round(num_b / den_b, 4),
            "odds_ratio": round(float(or_v), 3),
            "p_fisher": float(f"{pv:.3e}"),
        }

    # Ambiguous/Unresolved burden = P_LP + Conflicting + VUS out of all ClinVar-matched records
    non_ambig = tn["P_LP"] + tn["Conflicting"] + tn["VUS"]
    nfe_ambig = te["P_LP"] + te["Conflicting"] + te["VUS"]
    tests = {
        "ambiguous_or_plp_among_clinvar": fisher_rate(non_ambig, tn["clinvar_records"], nfe_ambig, te["clinvar_records"]),
        "conflicting_among_clinvar": fisher_rate(tn["Conflicting"], tn["clinvar_records"], te["Conflicting"], te["clinvar_records"]),
        "vus_among_clinvar": fisher_rate(tn["VUS"], tn["clinvar_records"], te["VUS"], te["clinvar_records"]),
        "plp_among_clinvar": fisher_rate(tn["P_LP"], tn["clinvar_records"], te["P_LP"], te["clinvar_records"]),
        "benign_resolved_among_clinvar": fisher_rate(tn["B_LB"], tn["clinvar_records"], te["B_LB"], te["clinvar_records"]),
    }

    out = {
        "totals": {k: dict(v) for k, v in totals.items()},
        "non_nfe_driver_populations": dict(driver_pops),
        "fisher_tests": tests,
        "per_gene": per_gene,
        "flagged_non_nfe_plp_examples": flagged_plp_examples,
        "flagged_non_nfe_conflict_examples": flagged_conflict_examples[:10],
    }
    (PORT / "pilots" / "pilot_a2_summary.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps({"totals": out["totals"], "driver_pops": out["non_nfe_driver_populations"],
                      "fisher_tests": tests, "n_plp_examples": len(flagged_plp_examples)}, indent=2))
    if flagged_plp_examples:
        print("PLP EXAMPLES:", json.dumps(flagged_plp_examples[:8], indent=2))


if __name__ == "__main__":
    main()
