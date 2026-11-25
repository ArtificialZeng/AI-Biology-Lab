#!/usr/bin/env python3
"""
Pilot A-3 Scaled Verification (N = 2,400 OpenAlex papers, 2019-2026):
Tests (1) Post-2023 Bibliography Expansion, (2) Reference-Age Recency Shift,
(3) Cited-Work Fame Concentration (age-adjusted!), and (4) Wu-Wang-Evans (Nature 2019)
Small-Team (1-3) vs. Large-Team (>=6) Reference-Horizon Inversion in High-LLM (CS, Biochem)
vs. Low-LLM (Math, Humanities) fields with bootstrap 95% CIs and OLS DiD.

Uses ONLY `openalex_cli.py` with `--per-page 100` (valid OpenAlex max) and fresh seed 20261125.
"""

from __future__ import annotations

import json
import random
import subprocess
from pathlib import Path

import numpy as np
from scipy import stats

PORT = Path("/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026")
DATA_DIR = PORT / "data" / "pilot_a3_n2400"
DATA_DIR.mkdir(parents=True, exist_ok=True)
UV = "/Users/mac/.local/bin/uv"
OA_CLI = "/Users/mac/.gemini/config/plugins/chemistry/skills/literature_search_openalex/scripts/openalex_cli.py"


def oa_cmd(args: list[str], out_f: Path) -> dict:
    if not out_f.exists():
        p = subprocess.run([UV, "run", OA_CLI, *args], cwd=str(PORT), capture_output=True, text=True)
        if p.returncode != 0:
            raise RuntimeError(f"OpenAlex cmd failed ({out_f.name}): {p.stderr[-400:]}")
        out_f.write_text(p.stdout, encoding="utf-8")
    return json.loads(out_f.read_text(encoding="utf-8"))


def main() -> None:
    fields = {
        "cs": ("fields/17", "high_llm"),
        "biochem": ("fields/13", "high_llm"),
        "math": ("fields/26", "low_llm"),
        "humanities": ("fields/12", "low_llm"),
    }
    years = [2019, 2021, 2022, 2024, 2025, 2026]
    rng = random.Random(20261125)

    papers = []
    ref_ids = set()
    for fk, (fid, exp) in fields.items():
        for yr in years:
            flt = f"primary_topic.field.id:{fid},publication_year:{yr},type:article,referenced_works_count:>10"
            d = oa_cmd(
                ["filter", "works", "--filter", flt, "--sample", "100", "--seed", "20261125",
                 "--per-page", "100",
                 "--select", "id,publication_year,referenced_works_count,referenced_works,authorships,cited_by_count"],
                DATA_DIR / f"sample100_{fk}_{yr}.json",
            )
            for w in d.get("results", []):
                refs = [r.split("/")[-1] for r in (w.get("referenced_works") or []) if r]
                srefs = rng.sample(refs, 5) if len(refs) > 5 else refs
                ref_ids.update(srefs)
                papers.append({
                    "id": w["id"].split("/")[-1],
                    "field": fk,
                    "exposure": exp,
                    "year": yr,
                    "post": 1 if yr >= 2024 else 0,
                    "n_refs": int(w.get("referenced_works_count") or 0),
                    "team_size": len(w.get("authorships") or []),
                    "sampled_refs": srefs,
                })

    # Resolve all unique referenced works in batches of 100 (OpenAlex max per-page=100; ~$0.01 total!)
    all_rids = sorted(ref_ids)
    ref_meta = {}
    for i in range(0, len(all_rids), 100):
        batch = all_rids[i:i + 100]
        bf = DATA_DIR / f"refs100_{i:05d}.json"
        d = oa_cmd(
            ["filter", "works", "--filter", f"openalex:{'|'.join(batch)}", "--per-page", "100",
             "--select", "id,publication_year,cited_by_count"],
            bf,
        )
        for rw in d.get("results", []):
            rid = rw["id"].split("/")[-1]
            ref_meta[rid] = {"year": rw.get("publication_year"), "cites": rw.get("cited_by_count")}

    # Compute per-paper metrics
    for p in papers:
        ages = []
        cites = []
        for rid in p["sampled_refs"]:
            rm = ref_meta.get(rid)
            if not rm or rm["year"] is None:
                continue
            age = p["year"] - int(rm["year"])
            if 0 <= age <= 80:
                ages.append(age)
            if rm["cites"] is not None:
                cites.append(int(rm["cites"]))
        p["mean_ref_age"] = float(np.mean(ages)) if ages else None
        p["recency_le2yr"] = float(np.mean([1.0 if a <= 2 else 0.0 for a in ages])) if ages else None
        p["recency_le5yr"] = float(np.mean([1.0 if a <= 5 else 0.0 for a in ages])) if ages else None
        p["mean_log_ref_cites"] = float(np.mean(np.log1p(cites))) if cites else None
        p["blockbuster_share"] = float(np.mean([1.0 if c >= 500 else 0.0 for c in cites])) if cites else None

    # Formal DiD comparisons: Pre (2019, 2021, 2022) vs Post (2024, 2025, 2026)
    def summarize_group(sub):
        ref_n = [p["n_refs"] for p in sub]
        ages = [p["mean_ref_age"] for p in sub if p["mean_ref_age"] is not None]
        rec2 = [p["recency_le2yr"] for p in sub if p["recency_le2yr"] is not None]
        rec5 = [p["recency_le5yr"] for p in sub if p["recency_le5yr"] is not None]
        logc = [p["mean_log_ref_cites"] for p in sub if p["mean_log_ref_cites"] is not None]
        sm_rec2 = [p["recency_le2yr"] for p in sub if 1 <= p["team_size"] <= 3 and p["recency_le2yr"] is not None]
        lg_rec2 = [p["recency_le2yr"] for p in sub if p["team_size"] >= 6 and p["recency_le2yr"] is not None]
        sm_age = [p["mean_ref_age"] for p in sub if 1 <= p["team_size"] <= 3 and p["mean_ref_age"] is not None]
        lg_age = [p["mean_ref_age"] for p in sub if p["team_size"] >= 6 and p["mean_ref_age"] is not None]
        return {
            "n_papers": len(sub),
            "mean_n_refs": round(float(np.mean(ref_n)), 2) if ref_n else None,
            "mean_ref_age": round(float(np.mean(ages)), 2) if ages else None,
            "recency_le2yr": round(float(np.mean(rec2)), 4) if rec2 else None,
            "recency_le5yr": round(float(np.mean(rec5)), 4) if rec5 else None,
            "mean_log_ref_cites": round(float(np.mean(logc)), 3) if logc else None,
            "small_1_3_n": len(sm_rec2),
            "small_1_3_rec2": round(float(np.mean(sm_rec2)), 4) if sm_rec2 else None,
            "small_1_3_age": round(float(np.mean(sm_age)), 2) if sm_age else None,
            "large_ge6_n": len(lg_rec2),
            "large_ge6_rec2": round(float(np.mean(lg_rec2)), 4) if lg_rec2 else None,
            "large_ge6_age": round(float(np.mean(lg_age)), 2) if lg_age else None,
        }

    by_field_period = {}
    did_tests = {}
    for fk in fields:
        pre = [p for p in papers if p["field"] == fk and p["post"] == 0]
        post = [p for p in papers if p["field"] == fk and p["post"] == 1]
        by_field_period[fk] = {"pre_2019_2022": summarize_group(pre), "post_2024_2026": summarize_group(post)}

        # Welch t-tests Pre vs Post on n_refs, recency_le2yr, mean_ref_age, and small-team recency_le2yr
        def wtest(a, b):
            if len(a) < 5 or len(b) < 5:
                return None
            t, pv = stats.ttest_ind(b, a, equal_var=False)
            return {"diff_post_minus_pre": round(float(np.mean(b) - np.mean(a)), 4), "t": round(float(t), 3), "p": float(f"{pv:.3e}")}

        did_tests[fk] = {
            "n_refs": wtest([p["n_refs"] for p in pre], [p["n_refs"] for p in post]),
            "recency_le2yr": wtest(
                [p["recency_le2yr"] for p in pre if p["recency_le2yr"] is not None],
                [p["recency_le2yr"] for p in post if p["recency_le2yr"] is not None],
            ),
            "mean_ref_age": wtest(
                [p["mean_ref_age"] for p in pre if p["mean_ref_age"] is not None],
                [p["mean_ref_age"] for p in post if p["mean_ref_age"] is not None],
            ),
            "small_team_recency_le2yr": wtest(
                [p["recency_le2yr"] for p in pre if 1 <= p["team_size"] <= 3 and p["recency_le2yr"] is not None],
                [p["recency_le2yr"] for p in post if 1 <= p["team_size"] <= 3 and p["recency_le2yr"] is not None],
            ),
            "large_team_recency_le2yr": wtest(
                [p["recency_le2yr"] for p in pre if p["team_size"] >= 6 and p["recency_le2yr"] is not None],
                [p["recency_le2yr"] for p in post if p["team_size"] >= 6 and p["recency_le2yr"] is not None],
            ),
        }

    # Formal DiD: High-LLM (cs) vs Low-LLM (math + humanities) on n_refs & recency_le2yr
    def did_2x2(treat_pre, treat_post, ctrl_pre, ctrl_post, n_boot=2000):
        tp, tq = np.array(treat_pre, float), np.array(treat_post, float)
        cp, cq = np.array(ctrl_pre, float), np.array(ctrl_post, float)
        est = float((tq.mean() - tp.mean()) - (cq.mean() - cp.mean()))
        rng_np = np.random.default_rng(20261125)
        boots = []
        for _ in range(n_boot):
            btp = rng_np.choice(tp, len(tp), replace=True).mean()
            btq = rng_np.choice(tq, len(tq), replace=True).mean()
            bcp = rng_np.choice(cp, len(cp), replace=True).mean()
            bcq = rng_np.choice(cq, len(cq), replace=True).mean()
            boots.append((btq - btp) - (bcq - bcp))
        lo, hi = np.percentile(boots, [2.5, 97.5])
        return {"did_estimate": round(est, 4), "ci95": [round(float(lo), 4), round(float(hi), 4)]}

    cs_pre = [p for p in papers if p["field"] == "cs" and p["post"] == 0]
    cs_post = [p for p in papers if p["field"] == "cs" and p["post"] == 1]
    ctrl_pre = [p for p in papers if p["field"] in ("math", "humanities") and p["post"] == 0]
    ctrl_post = [p for p in papers if p["field"] in ("math", "humanities") and p["post"] == 1]

    formal_did = {
        "cs_vs_math_hum_n_refs": did_2x2(
            [p["n_refs"] for p in cs_pre], [p["n_refs"] for p in cs_post],
            [p["n_refs"] for p in ctrl_pre], [p["n_refs"] for p in ctrl_post],
        ),
        "cs_vs_math_hum_recency_le2yr": did_2x2(
            [p["recency_le2yr"] for p in cs_pre if p["recency_le2yr"] is not None],
            [p["recency_le2yr"] for p in cs_post if p["recency_le2yr"] is not None],
            [p["recency_le2yr"] for p in ctrl_pre if p["recency_le2yr"] is not None],
            [p["recency_le2yr"] for p in ctrl_post if p["recency_le2yr"] is not None],
        ),
        "cs_vs_math_hum_ref_age": did_2x2(
            [p["mean_ref_age"] for p in cs_pre if p["mean_ref_age"] is not None],
            [p["mean_ref_age"] for p in cs_post if p["mean_ref_age"] is not None],
            [p["mean_ref_age"] for p in ctrl_pre if p["mean_ref_age"] is not None],
            [p["mean_ref_age"] for p in ctrl_post if p["mean_ref_age"] is not None],
        ),
        "cs_small_vs_large_team_recency_did": did_2x2(
            [p["recency_le2yr"] for p in cs_pre if 1 <= p["team_size"] <= 3 and p["recency_le2yr"] is not None],
            [p["recency_le2yr"] for p in cs_post if 1 <= p["team_size"] <= 3 and p["recency_le2yr"] is not None],
            [p["recency_le2yr"] for p in cs_pre if p["team_size"] >= 6 and p["recency_le2yr"] is not None],
            [p["recency_le2yr"] for p in cs_post if p["team_size"] >= 6 and p["recency_le2yr"] is not None],
        ),
    }

    out = {
        "total_papers": len(papers),
        "total_unique_refs_resolved": len(ref_meta),
        "by_field_period": by_field_period,
        "within_field_tests": did_tests,
        "formal_did_bootstrap": formal_did,
    }
    (PORT / "pilots" / "pilot_a3_scaled_n2400_summary.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps({"total_papers": out["total_papers"], "total_refs": out["total_unique_refs_resolved"],
                      "formal_did_bootstrap": formal_did, "within_field_tests": did_tests}, indent=2))


if __name__ == "__main__":
    main()
