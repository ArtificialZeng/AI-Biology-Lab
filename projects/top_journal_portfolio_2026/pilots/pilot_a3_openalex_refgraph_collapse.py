#!/usr/bin/env python3
"""
Pilot A-3 (Nature / Science / Nature Human Behaviour):
Reference-Graph Horizon Collapse & Team-Size Disruption Compression in the
Generative-AI Era (OpenAlex 2019-2026 Multi-Field DiD Panel).

Measures REAL published reference-list structures (reference age, recency <=2yr share,
cited-work fame concentration, and Small-Team 1-3 vs Large-Team 8+ divergence) across
High-LLM-Exposure fields (Computer Science, Molecular Biology) vs. Low-LLM-Exposure
control fields (Mathematics, History).

Uses ONLY the official literature-search-openalex CLI wrapper (`openalex_cli.py`).
"""

from __future__ import annotations

import json
import random
import statistics
import subprocess
from pathlib import Path

PORT = Path("/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026")
DATA_DIR = PORT / "data" / "pilot_a3"
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
    # 1. Resolve canonical OpenAlex Field IDs first (per OpenAlex Core Rule #2: resolve before filter!)
    field_queries = {
        "cs": ("Computer Science", "high_llm"),
        "biochem": ("Biochemistry, Genetics and Molecular Biology", "high_llm"),
        "math": ("Mathematics", "low_llm"),
        "humanities": ("Arts and Humanities", "low_llm"),
    }
    resolved_fields = {}
    for fk, (qname, exposure) in field_queries.items():
        d = oa_cmd(["resolve", "fields", qname, "--per-page", "3"], DATA_DIR / f"resolve_{fk}.json")
        cands = d if isinstance(d, list) else d.get("results", [])
        if not cands:
            raise RuntimeError(f"Could not resolve field {qname}")
        fid = cands[0]["id"].replace("https://openalex.org/", "")
        resolved_fields[fk] = {"id": fid, "display_name": cands[0].get("display_name", qname), "exposure": exposure}

    # 2. Sample 120 published journal articles per (field, year) for years 2019, 2021, 2022, 2024, 2025, 2026
    years = [2019, 2021, 2022, 2024, 2025, 2026]
    papers_by_fy = {}
    ref_ids_to_resolve = set()
    rng = random.Random(20261124)

    for fk, finfo in resolved_fields.items():
        fid = finfo["id"]
        for yr in years:
            flt = f"primary_topic.field.id:{fid},publication_year:{yr},type:article,referenced_works_count:>10"
            d = oa_cmd(
                ["filter", "works", "--filter", flt, "--sample", "120", "--seed", "20261124",
                 "--per-page", "120",
                 "--select", "id,publication_year,referenced_works_count,referenced_works,authorships,cited_by_count"],
                DATA_DIR / f"sample_{fk}_{yr}.json",
            )
            works = d.get("results", [])
            papers_by_fy[(fk, yr)] = works
            # Pick up to 6 random references per paper for compact batch resolution (budget-friendly!)
            for w in works:
                refs = w.get("referenced_works") or []
                short_refs = [r.split("/")[-1] for r in refs if r]
                if len(short_refs) > 6:
                    short_refs = rng.sample(short_refs, 6)
                w["_sampled_refs"] = short_refs
                ref_ids_to_resolve.update(short_refs)

    # 3. Batch-resolve sampled referenced works (50 per filter call, only $0.0001 each!)
    all_ref_ids = sorted(ref_ids_to_resolve)
    # Cap at 6,000 unique refs (120 calls * $0.0001 = $0.012 total!)
    if len(all_ref_ids) > 6000:
        all_ref_ids = rng.sample(all_ref_ids, 6000)
    ref_meta = {}
    for i in range(0, len(all_ref_ids), 50):
        batch = all_ref_ids[i:i + 50]
        pipe_ids = "|".join(batch)
        bf = DATA_DIR / f"refbatch_{i:05d}.json"
        d = oa_cmd(
            ["filter", "works", "--filter", f"openalex:{pipe_ids}", "--per-page", "50",
             "--select", "id,publication_year,cited_by_count"],
            bf,
        )
        for rw in d.get("results", []):
            rid = rw["id"].split("/")[-1]
            ref_meta[rid] = {
                "year": rw.get("publication_year"),
                "cites": rw.get("cited_by_count"),
            }

    # 4. Aggregate structural reference-graph & team-size metrics by (field, year)
    panel_stats = {}
    for (fk, yr), works in sorted(papers_by_fy.items()):
        ages = []
        recent_2yr_flags = []
        ref_cites = []
        ref_counts = []
        team_sizes = []
        small_team_recent = []
        large_team_recent = []
        small_team_cites = []
        large_team_cites = []

        for w in works:
            n_auth = len(w.get("authorships") or [])
            team_sizes.append(n_auth)
            ref_counts.append(int(w.get("referenced_works_count") or 0))
            p_ages = []
            p_cites = []
            for rid in w.get("_sampled_refs", []):
                rm = ref_meta.get(rid)
                if not rm or rm["year"] is None:
                    continue
                age = yr - int(rm["year"])
                if 0 <= age <= 80:
                    p_ages.append(age)
                    ages.append(age)
                    recent_2yr_flags.append(1 if age <= 2 else 0)
                if rm["cites"] is not None:
                    c = int(rm["cites"])
                    p_cites.append(c)
                    ref_cites.append(c)
            if p_ages:
                rec_share = sum(1 for a in p_ages if a <= 2) / len(p_ages)
                if 1 <= n_auth <= 3:
                    small_team_recent.append(rec_share)
                elif n_auth >= 6:
                    large_team_recent.append(rec_share)
            if p_cites:
                med_c = statistics.median(p_cites)
                if 1 <= n_auth <= 3:
                    small_team_cites.append(med_c)
                elif n_auth >= 6:
                    large_team_cites.append(med_c)

        panel_stats.setdefault(fk, {})[str(yr)] = {
            "n_papers": len(works),
            "n_refs_resolved": len(ages),
            "mean_ref_count": round(statistics.mean(ref_counts), 1) if ref_counts else None,
            "mean_team_size": round(statistics.mean(team_sizes), 2) if team_sizes else None,
            "median_ref_age": round(statistics.median(ages), 2) if ages else None,
            "recency_le_2yr_share": round(sum(recent_2yr_flags) / len(recent_2yr_flags), 4) if recent_2yr_flags else None,
            "median_cited_work_cites": round(statistics.median(ref_cites), 1) if ref_cites else None,
            "blockbuster_gt500_share": round(sum(1 for c in ref_cites if c >= 500) / len(ref_cites), 4) if ref_cites else None,
            "small_team_1_3_n": len(small_team_recent),
            "small_team_recency_le2yr": round(statistics.mean(small_team_recent), 4) if small_team_recent else None,
            "large_team_ge6_n": len(large_team_recent),
            "large_team_recency_le2yr": round(statistics.mean(large_team_recent), 4) if large_team_recent else None,
            "small_team_med_ref_cites": round(statistics.median(small_team_cites), 1) if small_team_cites else None,
            "large_team_med_ref_cites": round(statistics.median(large_team_cites), 1) if large_team_cites else None,
        }

    out = {
        "resolved_fields": resolved_fields,
        "total_unique_refs_resolved": len(ref_meta),
        "panel_by_field_year": panel_stats,
    }
    (PORT / "pilots" / "pilot_a3_summary.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(out, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
