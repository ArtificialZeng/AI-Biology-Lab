#!/usr/bin/env python3
"""
Pilot A-1 (Nature Biotechnology / Science Translational Medicine):
Modality-Decomposed Target Safety Frontier — Separating On-Target Germline Pleiotropy
& Evolutionary Constraint (OpenTargets + gnomAD LOEUF) from Off-Target Chemical
Promiscuity & Affinity Inversions (ChEMBL) Across FDA-Approved Monoclonal Antibodies
vs. Small Molecules (Black Box Warnings).

Uses ONLY official skill CLI wrappers (`chembl_api.py`, `get_gene_constraint.py`, `query_opentargets.py`).
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

PORT = Path("/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026")
DATA_DIR = PORT / "data" / "pilot_a1"
DATA_DIR.mkdir(parents=True, exist_ok=True)
UV = "/Users/mac/.local/bin/uv"
CHEMBL_CLI = "/Users/mac/.gemini/config/plugins/chemistry/skills/chembl_database/scripts/chembl_api.py"
GNOMAD_CONST = "/Users/mac/.gemini/config/plugins/genomics/skills/gnomad_database/scripts/get_gene_constraint.py"
OT_CLI = "/Users/mac/.gemini/config/plugins/genomics/skills/opentargets_database/scripts/query_opentargets.py"


def run_chembl(args: list[str], out_f: Path) -> dict:
    if not out_f.exists():
        p = subprocess.run([UV, "run", CHEMBL_CLI, *args, "--output", str(out_f)],
                           cwd=str(PORT), capture_output=True, text=True)
        if p.returncode != 0:
            return {"error": p.stderr[-300:]}
    return json.loads(out_f.read_text(encoding="utf-8"))


def get_loeuf(gene: str) -> dict:
    out_f = DATA_DIR / f"gnomad_const_{gene}.json"
    if not out_f.exists():
        p = subprocess.run([UV, "run", GNOMAD_CONST, "--gene", gene, "--output", str(out_f)],
                           cwd=str(PORT), capture_output=True, text=True)
        if p.returncode != 0:
            return {}
    try:
        d = json.loads(out_f.read_text(encoding="utf-8"))
        g = d.get("data", {}).get("gene") or d
        c = g.get("gnomad_constraint") or g.get("constraint") or {}
        return {
            "loeuf": c.get("oe_lof_upper"),
            "pli": c.get("pli"),
            "ensembl_id": g.get("gene_id") or g.get("ensembl_id"),
        }
    except Exception:
        return {}


def get_ot_pleiotropy(ensembl_id: str) -> int | None:
    if not ensembl_id:
        return None
    out_f = DATA_DIR / f"ot_dis_{ensembl_id}.json"
    if not out_f.exists():
        p = subprocess.run([UV, "run", OT_CLI, "--output", str(out_f), "--limit", "5",
                            "get-associated-diseases", ensembl_id],
                           cwd=str(PORT), capture_output=True, text=True)
        if p.returncode != 0:
            return None
    try:
        d = json.loads(out_f.read_text(encoding="utf-8"))
        # Handle nested or flat count structures
        if isinstance(d, dict):
            if "count" in d:
                return int(d["count"])
            t = d.get("data", {}).get("target") or d.get("target") or d
            ad = t.get("associatedDiseases") or {}
            if "count" in ad:
                return int(ad["count"])
    except Exception:
        pass
    return None


def main() -> None:
    # 1. Fetch approved Antibodies & approved Small Molecules from ChEMBL
    mab_d = run_chembl(["molecule", "--filter", "max_phase=4", "molecule_type=Antibody", "--limit", "250"],
                       DATA_DIR / "approved_mabs.json")
    sm_d = run_chembl(["molecule", "--filter", "max_phase=4", "molecule_type=Small molecule", "--limit", "500"],
                      DATA_DIR / "approved_sms.json")

    mabs = mab_d.get("molecules", [])
    sms = sm_d.get("molecules", [])

    def bbw_rate(mols):
        tot = len(mols)
        bbw = sum(1 for m in mols if m.get("black_box_warning") == 1)
        wd = sum(1 for m in mols if m.get("withdrawn_flag"))
        return {"n": tot, "bbw_count": bbw, "bbw_rate": round(bbw / tot, 4) if tot else None, "withdrawn": wd}

    # 2. Pick stratified samples (14 BBW + 14 non-BBW for mAbs; 14 BBW + 14 non-BBW for Small Molecules)
    #    with recent approvals (first_approval >= 2000 for cleaner ChEMBL target annotations)
    def pick_stratified(mols, n_each=14):
        pos = [m for m in mols if m.get("black_box_warning") == 1 and (m.get("first_approval") or 0) >= 2000]
        neg = [m for m in mols if m.get("black_box_warning") == 0 and (m.get("first_approval") or 0) >= 2000]
        pos = sorted(pos, key=lambda x: -(x.get("first_approval") or 0))[:n_each]
        neg = sorted(neg, key=lambda x: -(x.get("first_approval") or 0))[:n_each]
        return pos + neg

    sample_mabs = pick_stratified(mabs, 14)
    sample_sms = pick_stratified(sms, 14)

    # 3. Map each drug to its primary human target gene symbol via ChEMBL mechanism -> target
    def map_drug_profile(mol: dict, is_sm: bool) -> dict:
        cid = mol["molecule_chembl_id"]
        mech = run_chembl(["mechanism", "--filter", f"molecule_chembl_id={cid}", "--limit", "10"],
                          DATA_DIR / f"mech_{cid}.json")
        mechs = mech.get("mechanisms", [])
        target_ids = sorted({m["target_chembl_id"] for m in mechs if m.get("target_chembl_id")})
        primary_gene = None
        primary_tid = None
        for tid in target_ids:
            td = run_chembl(["target", "--id", tid], DATA_DIR / f"tgt_{tid}.json")
            if td.get("organism") != "Homo sapiens":
                continue
            for comp in td.get("target_components", []):
                for syn in comp.get("target_component_synonyms", []):
                    if syn.get("syn_type") == "GENE_SYMBOL":
                        primary_gene = syn.get("component_synonym")
                        primary_tid = tid
                        break
                if primary_gene:
                    break
            if primary_gene:
                break

        const = get_loeuf(primary_gene) if primary_gene else {}
        pleio = get_ot_pleiotropy(const.get("ensembl_id")) if const.get("ensembl_id") else None

        # For small molecules, quantify ChEMBL human sub-uM targets & affinity inversion!
        n_sub_um_targets = None
        affinity_inversion = None
        primary_best_nm = None
        offtarget_best_nm = None
        if is_sm:
            act_d = run_chembl(
                ["activity", "--filter", f"molecule_chembl_id={cid}", "target_organism=Homo sapiens",
                 "--normalize", "--limit", "200"],
                DATA_DIR / f"act_{cid}.json",
            )
            best_by_tgt = {}
            for a in act_d.get("activities", []):
                if a.get("standard_type") not in ("IC50", "Ki", "Kd", "EC50"):
                    continue
                if a.get("standard_relation") not in ("=", "<", "<="):
                    continue
                nm = a.get("normalized_value_nM")
                tid = a.get("target_chembl_id")
                if nm is None or tid is None or nm <= 0:
                    continue
                if tid not in best_by_tgt or nm < best_by_tgt[tid]:
                    best_by_tgt[tid] = float(nm)
            sub_um = {t: v for t, v in best_by_tgt.items() if v <= 1000.0}
            n_sub_um_targets = len(sub_um)
            if primary_tid and primary_tid in best_by_tgt:
                primary_best_nm = best_by_tgt[primary_tid]
                off_vals = [v for t, v in best_by_tgt.items() if t not in set(target_ids)]
                if off_vals:
                    offtarget_best_nm = min(off_vals)
                    # Strict affinity inversion: off-target <= primary_best_nm (equal or higher potency!)
                    affinity_inversion = bool(offtarget_best_nm <= primary_best_nm)

        return {
            "chembl_id": cid,
            "pref_name": mol.get("pref_name"),
            "modality": "Small molecule" if is_sm else "Antibody",
            "first_approval": mol.get("first_approval"),
            "black_box_warning": int(mol.get("black_box_warning") or 0),
            "primary_target_chembl": primary_tid,
            "primary_gene": primary_gene,
            "loeuf": const.get("loeuf"),
            "pli": const.get("pli"),
            "ot_disease_breadth": pleio,
            "n_sub_um_human_targets": n_sub_um_targets,
            "primary_best_nm": primary_best_nm,
            "offtarget_best_nm": offtarget_best_nm,
            "affinity_inversion": affinity_inversion,
        }

    mab_profiles = [map_drug_profile(m, is_sm=False) for m in sample_mabs]
    sm_profiles = [map_drug_profile(m, is_sm=True) for m in sample_sms]

    def grp_stats(profiles, bbw_val):
        sub = [p for p in profiles if p["black_box_warning"] == bbw_val]
        loeufs = [p["loeuf"] for p in sub if p["loeuf"] is not None]
        pleios = [p["ot_disease_breadth"] for p in sub if p["ot_disease_breadth"] is not None]
        subums = [p["n_sub_um_human_targets"] for p in sub if p["n_sub_um_human_targets"] is not None]
        invs = [p["affinity_inversion"] for p in sub if p["affinity_inversion"] is not None]
        return {
            "n_total": len(sub),
            "n_mapped_gene": sum(1 for p in sub if p["primary_gene"]),
            "n_loeuf": len(loeufs),
            "median_loeuf": round(sorted(loeufs)[len(loeufs)//2], 4) if loeufs else None,
            "mean_loeuf": round(sum(loeufs)/len(loeufs), 4) if loeufs else None,
            "n_pleio": len(pleios),
            "median_ot_diseases": sorted(pleios)[len(pleios)//2] if pleios else None,
            "mean_ot_diseases": round(sum(pleios)/len(pleios), 1) if pleios else None,
            "n_subum": len(subums),
            "mean_sub_um_targets": round(sum(subums)/len(subums), 2) if subums else None,
            "n_inv_eval": len(invs),
            "affinity_inversion_rate": round(sum(1 for x in invs if x)/len(invs), 3) if invs else None,
        }

    out = {
        "population_chembl_approved": {
            "antibodies": bbw_rate(mabs),
            "small_molecules_page1": bbw_rate(sms),
        },
        "stratified_summary": {
            "Antibody_BBW1": grp_stats(mab_profiles, 1),
            "Antibody_BBW0": grp_stats(mab_profiles, 0),
            "SmallMolecule_BBW1": grp_stats(sm_profiles, 1),
            "SmallMolecule_BBW0": grp_stats(sm_profiles, 0),
        },
        "profiles": mab_profiles + sm_profiles,
    }
    (PORT / "pilots" / "pilot_a1_summary.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(out["population_chembl_approved"], indent=2))
    print(json.dumps(out["stratified_summary"], indent=2))


if __name__ == "__main__":
    main()
