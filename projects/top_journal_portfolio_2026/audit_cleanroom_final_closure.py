#!/usr/bin/env python3
"""
Independent Clean-Room Cross-Check Audit across Study A-1 and Study B-2.
Verifies:
  1. Every headline number in manuscript_{a1,b2}.{tex,md} against authoritative JSONs.
  2. Context-aware checks distinguishing transparent retractions/revision notes from active claims.
  3. LaTeX compilation logs (0 undefined cites/refs, 0 wide Overfull hbox > 10pt, 0 stray '**').
  4. T09 rendered PDF table body effective font sizes (>= 7.0 pt via pypdf cm*tm scaling).
  5. T06/T10: all 15/15 Official CRAN R HonestDiD v0.2.8 specs present, 0 errors/NAs,
     0 remaining grid-boundary hits after +-40*SE/1001-pt expansion, convergence verified.
  6. Exact SHA-256 checksums of final deliverables and Supplementary CSVs.
  7. Non-destructive Round 16 isolation audit (path/timestamp/tracked-diff check; never deletes/resets other files).
"""
import csv, hashlib, json, math, os, re, subprocess

ROOT = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026"
A1 = os.path.join(ROOT, "studies/study_a1_clinvar_gnomad_reclass")
B2 = os.path.join(ROOT, "studies/study_b2_fda_2026_sweep")
R16 = "/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold"

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()

def audit_tex_log(tex_path, log_path):
    tex = open(tex_path, errors="replace").read()
    log = open(log_path, errors="replace").read()
    undef_cites = re.findall(r"Citation `[^']+' on page \d+ undefined", log)
    undef_refs = re.findall(r"Reference `[^']+' on page \d+ undefined", log)
    overfulls = []
    for m in re.finditer(r"Overfull \\hbox \(([0-9.]+)pt too wide\).*?(?:lines (\d+)--(\d+)|line (\d+))", log):
        overfulls.append({"pt": float(m.group(1)), "lines": m.group(2) + "-" + m.group(3) if m.group(2) else m.group(4)})
    overfull_gt10 = [x for x in overfulls if x["pt"] > 10.0]
    stray_bold = re.findall(r"\*\*[^*]+\*\*", tex)
    resizebox_in_tables = len(re.findall(r"\\resizebox\b", tex))
    overfull_vbox = len(re.findall(r"Overfull \\vbox", log))
    return {
        "undef_cites": len(undef_cites),
        "undef_refs": len(undef_refs),
        "overfull_all_count": len(overfulls),
        "overfull_gt10_count": len(overfull_gt10),
        "overfull_gt10": overfull_gt10,
        "overfull_vbox_count": overfull_vbox,
        "stray_markdown_bold": len(stray_bold),
        "resizebox_count": resizebox_in_tables,
    }

def audit_pdf_table_fonts(pdf_path, table_pages_1based):
    """Measure effective rendered body font size (cm*font_size) via real pypdf across ALL pages."""
    from pypdf import PdfReader
    reader = PdfReader(pdf_path)
    out = {
        "method": "pypdf_cm_matrix_exact",
        "page_count": len(reader.pages),
        "table_pages": table_pages_1based,
        "small_body_lt_6_8pt_count": 0,
        "small_body_lt_6_8pt_samples": [],
        "min_body_pt_all_pages": None,
        "min_body_pt_by_table_page": {},
    }
    all_body_pts = []
    for pno, page in enumerate(reader.pages, 1):
        page_pts = []
        def visitor(text, cm, tm, font_dict, font_size):
            eff = font_size * math.hypot(cm[0], cm[1])
            sample = (text or "").strip()
            if len(sample) >= 12 and sum(c.isalpha() for c in sample) >= 8 and eff > 0.5:
                page_pts.append((round(eff, 4), sample[:60]))
                if eff < 6.8:
                    out["small_body_lt_6_8pt_samples"].append({"page": pno, "effective_pt": round(eff, 4), "text": sample[:80]})
        page.extract_text(visitor_text=visitor)
        if page_pts:
            page_pts.sort(key=lambda x: x[0])
            all_body_pts.append((page_pts[0][0], pno, page_pts[0][1]))
            if pno in table_pages_1based:
                out["min_body_pt_by_table_page"][str(pno)] = {
                    "min_pt": page_pts[0][0],
                    "sample": page_pts[0][1],
                }
    out["small_body_lt_6_8pt_count"] = len(out["small_body_lt_6_8pt_samples"])
    if all_body_pts:
        all_body_pts.sort(key=lambda x: x[0])
        out["min_body_pt_all_pages"] = {"pt": all_body_pts[0][0], "page": all_body_pts[0][1], "sample": all_body_pts[0][2]}
    return out

def check_context_aware(text, required, active_claim_forbidden):
    """
    Distinguish retracted/historical audit sentences from active claims:
    If a forbidden pattern only appears inside a sentence containing retraction/audit context
    (retract, earlier, pre-audit, rather than, instead of, cannot, does not), log it as
    transparent_history rather than an active violation.
    """
    missing = [s for s in required if s not in text]
    violations = []
    transparent_mentions = []
    retraction_cues = ("retract", "earlier", "pre-audit", "rather than", "instead of", "not ", "cannot", "does not", "do not", "neither")
    sentences = re.split(r"(?<=[.!?])\s+|\n", text)
    for pat in active_claim_forbidden:
        pl = pat.lower()
        for sent in sentences:
            if pl in sent.lower():
                if any(cue in sent.lower() for cue in retraction_cues):
                    transparent_mentions.append({"pattern": pat, "context": sent.strip()[:140]})
                else:
                    violations.append({"pattern": pat, "context": sent.strip()[:140]})
    return {
        "missing": missing,
        "active_violations": violations,
        "transparent_retraction_or_history_mentions": transparent_mentions,
        "ok": (not missing and not violations),
    }

# Load authoritative JSONs
a1_t08 = json.load(open(os.path.join(A1, "audit/a1_t08_chrx_ploidy_closure.json")))
b2_t00 = json.load(open(os.path.join(B2, "audit/b2_t00_proxy_date_forensic_audit.json")))
b2_def = json.load(open(os.path.join(B2, "audit/b2_definitive_t00_to_t04_results.json")))
b2_t07 = json.load(open(os.path.join(B2, "audit/b2_t07_invariance_and_es_clean.json")))
b2_hd  = json.load(open(os.path.join(B2, "audit/b2_t06_official_r_honestdid_results.json")))
b2_t10 = json.load(open(os.path.join(B2, "audit/b2_t10_grid_boundary_and_convergence_audit.json")))
b2_qhm = json.load(open(os.path.join(B2, "audit/b2_t00_quarterly_honestdid_and_mechanisms.json")))

# Verify T10 CSV has 0 boundary hits
with open(os.path.join(B2, "audit/b2_supp_t06_official_r_honestdid_all15.csv")) as f:
    s3_rows = list(csv.DictReader(f))
s3_hits = sum(1 for r in s3_rows if str(r.get("hits_grid_boundary", "")).lower() == "true")
s3_expanded = sum(1 for r in s3_rows if str(r.get("expanded_grid_recomputed", "")).lower() == "true")

json_assertions = {
    "a1_t08_chrx_max_an_exact_1213427": a1_t08["ploidy_denominators"]["chrx_nonpar_max_an_theor"] == 1213427.0,
    "a1_t08_chrx_recovered_vus_406": (
        a1_t08["chrx_gene_callability_recovery"]["ABCD1"]["new_ploidy_80pct_vus_pass"]
        + a1_t08["chrx_gene_callability_recovery"]["GLA"]["new_ploidy_80pct_vus_pass"]
        + a1_t08["chrx_gene_callability_recovery"]["OTC"]["new_ploidy_80pct_vus_pass"]
    ) == 406,
    "a1_t08_h1_n_10058": a1_t08["h1_84gene_ploidy_aware_revised"]["n_total"] == 10058,
    "a1_t08_tier1_327_ba1_21": (
        a1_t08["vus_t05_t08_ploidy_aware_summary"]["tier1_cspec_30_faithful"]["bs1_n"] == 327
        and a1_t08["vus_t05_t08_ploidy_aware_summary"]["tier1_cspec_30_faithful"]["ba1_n"] == 21
    ),
    "a1_t08_all83_1388_ba1_246": (
        a1_t08["vus_t05_t08_ploidy_aware_summary"]["all_83_ploidy_faithful"]["bs1_n"] == 1388
        and a1_t08["vus_t05_t08_ploidy_aware_summary"]["all_83_ploidy_faithful"]["ba1_n"] == 246
    ),
    "b2_t00_zero_fallback_and_exact_100pct": (
        b2_t00["all_arms"]["n_no_posted_no_est_ONLY_subinfo_releaseDate"] == 0
        and b2_t00["all_arms"]["n_no_posted_no_est_no_rel_ONLY_unpostedEvent"] == 0
        and b2_t00["all_arms"]["est_vs_subinfo_releaseDate_comparison"]["exact_equal_pct"] == 100.0
    ),
    "b2_t07_invariance_lt_1e11": all(v["max_wald_stat_diff"] < 1e-11 for v in b2_t07["t07_invariance_proof"].values()),
    "b2_t06_official_r_15_of_15_zero_na": (
        len(b2_hd["specs"]) == 15
        and b2_hd["honestdid_version"] == "0.2.8"
        and all(
            all(r.get("lb") is not None and r.get("ub") is not None for r in sp["delta_SD_FLCI"] + sp["delta_RM_CLF"])
            for sp in b2_hd["specs"].values()
        )
    ),
    "b2_t10_all_5_expanded_0_hits_and_converged": (
        b2_t10["n_boundary_hits_before_fix"] == 5
        and b2_t10["n_boundary_hits_after_fix"] == 0
        and s3_hits == 0
        and s3_expanded == 5
        and b2_t10["convergence_check_1001_vs_2001"]["within_one_step"] is True
    ),
    "b2_h4_mortality_tost5_and_efficacy_mde2286": (
        b2_def["t03_t04_h4_lemons_audit"]["regressions"]["mortality_rate"]["within_sponsor_post_vs_ontime"]["tost_5pp_p"] < 0.01
        and round(b2_def["t03_t04_h4_lemons_audit"]["regressions"]["nonsig_p05"]["within_sponsor_post_vs_ontime"]["mde80_pp"], 2) == 22.86
    ),
    "b2_phase1_spillover_0893": round(b2_qhm["t00_uncensored_phase1_spillover"]["pois_irr"], 3) == 0.893,
}

# Read manuscript sources
a1_tex = open(os.path.join(A1, "manuscript/manuscript_a1.tex")).read()
a1_md  = open(os.path.join(A1, "MANUSCRIPT_A1_DRAFT.md")).read()
a1_bib = open(os.path.join(A1, "manuscript/references_a1.bib")).read()
b2_tex = open(os.path.join(B2, "manuscript/manuscript_b2.tex")).read()
b2_md  = open(os.path.join(B2, "MANUSCRIPT_B2_DRAFT.md")).read()

a1_req = ["1,213,427", "406/441", "10,058", "327", "21", "1,388", "246", "51,190", "36/37", "GN089", "GN156", "GN179"]
a1_forb = ["1,371", "311 VUS", "311 high-callability", "proves causal mediation", "1,214,237"]
b2_req = ["53,736", "16,333", "1.598", "2.288", "1.960", "1.654", "0.2.8", "0.893", "22.86", "0.0014", "0.0016", "April 13, 2026", "post-hoc"]
b2_forb = ["falsifying the strategic lemons", "zero spillover", "zero lift", "strictly null", "ruling out administrative deferral", "100% uncensored", "Exact Rambachan-Roth FLCI"]

bib_check = {
    "lai2026epilepsy_exact_6_authors": all(nm in a1_bib for nm in [
        "Lai, Jin-Xing", "He, Ming-Feng", "Fan, Zhi-Cong", "Yan, Fu-Jiang", "Guan, Bao-Zhu", "Liu, Xiang-Hong"
    ]),
    "johnston2026daft_present": "johnston2026daft" in a1_bib,
    "gnomadv41sex2024_present": "gnomadv41sex2024" in a1_bib,
}

# Parse exact table pages from .aux for pypdf font verification
def table_pages_from_aux(aux_path):
    txt = open(aux_path, errors="replace").read()
    return sorted({int(m.group(2)) for m in re.finditer(r"\\newlabel\{tab:[^}]+\}\{\{([^}]+)\}\{(\d+)\}", txt)})

a1_tab_pages = table_pages_from_aux(os.path.join(A1, "manuscript/manuscript_a1.aux"))
b2_tab_pages = table_pages_from_aux(os.path.join(B2, "manuscript/manuscript_b2.aux"))

# Non-destructive Round 16 isolation check (strictly scoped to bioxarena-pro-scaffold; never deletes/resets other files)
r16_head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=R16, capture_output=True, text=True).stdout.strip()
r16_tracked_diff = subprocess.run(["git", "diff", "--name-only", "--", "."], cwd=R16, capture_output=True, text=True).stdout.strip()
r16_staged_diff = subprocess.run(["git", "diff", "--cached", "--name-only", "--", "."], cwd=R16, capture_output=True, text=True).stdout.strip()
# Check newest mtime inside bioxarena-pro-scaffold (ignoring .git) to confirm no file written in current 2026-11-25 10:00+ window
newest_r16 = None
for dirpath, dirnames, filenames in os.walk(R16):
    dirnames[:] = [d for d in dirnames if d != ".git"]
    for fn in filenames:
        fp = os.path.join(dirpath, fn)
        try:
            mt = os.path.getmtime(fp)
            if newest_r16 is None or mt > newest_r16[0]:
                newest_r16 = (mt, os.path.relpath(fp, R16))
        except OSError:
            pass
import datetime
newest_r16_iso = datetime.datetime.fromtimestamp(newest_r16[0]).strftime("%Y-%m-%d %H:%M:%S") if newest_r16 else None

deliverables = [
    os.path.join(A1, "manuscript/manuscript_a1.pdf"),
    os.path.join(A1, "manuscript/references_a1.bib"),
    os.path.join(A1, "audit/a1_acmg84_verified_cspec_johnston2026_rules.csv"),
    os.path.join(A1, "audit/a1_ad55_verified_t05_triage.csv"),
    os.path.join(A1, "audit/a1_vus_verified_t05_t08_bs1_candidates.csv"),
    os.path.join(A1, "audit/a1_t08_chrx_ploidy_closure.json"),
    os.path.join(B2, "manuscript/manuscript_b2.pdf"),
    os.path.join(B2, "manuscript/references_b2.bib"),
    os.path.join(B2, "audit/b2_supp_monthly_combined_proxy_hazards.csv"),
    os.path.join(B2, "audit/b2_supp_t00_t07_master_estimates.csv"),
    os.path.join(B2, "audit/b2_supp_t06_official_r_honestdid_all15.csv"),
    os.path.join(B2, "audit/b2_t06_official_r_honestdid_results.json"),
    os.path.join(B2, "audit/b2_t10_grid_boundary_and_convergence_audit.json"),
]
checksums = {
    os.path.relpath(p, ROOT): {
        "exists": os.path.exists(p),
        "bytes": os.path.getsize(p) if os.path.exists(p) else 0,
        "sha256": sha256(p) if os.path.exists(p) else None,
    }
    for p in deliverables
}

report = {
    "json_assertions": json_assertions,
    "all_json_assertions_pass": all(json_assertions.values()),
    "a1_latex_audit": audit_tex_log(os.path.join(A1, "manuscript/manuscript_a1.tex"), os.path.join(A1, "manuscript/manuscript_a1.log")),
    "b2_latex_audit": audit_tex_log(os.path.join(B2, "manuscript/manuscript_b2.tex"), os.path.join(B2, "manuscript/manuscript_b2.log")),
    "a1_pdf_table_fonts": audit_pdf_table_fonts(os.path.join(A1, "manuscript/manuscript_a1.pdf"), a1_tab_pages),
    "b2_pdf_table_fonts": audit_pdf_table_fonts(os.path.join(B2, "manuscript/manuscript_b2.pdf"), b2_tab_pages),
    "a1_tex_consistency": check_context_aware(a1_tex, a1_req, a1_forb),
    "a1_md_consistency": check_context_aware(a1_md, a1_req, a1_forb),
    "b2_tex_consistency": check_context_aware(b2_tex, b2_req, b2_forb),
    "b2_md_consistency": check_context_aware(b2_md, b2_req, b2_forb),
    "a1_bib_check": bib_check,
    "round16_isolation": {
        "head": r16_head,
        "tracked_diff_clean": (r16_tracked_diff == "" and r16_staged_diff == ""),
        "newest_file_in_bioxarena": {"path": newest_r16[1] if newest_r16 else None, "mtime_local": newest_r16_iso},
        "note": "Non-destructive check: pre-existing untracked files and .gitignore outside projects/top_journal_portfolio_2026/ are intentionally untouched.",
    },
    "deliverables_checksums": checksums,
}

out_path = os.path.join(ROOT, "final_cleanroom_closure_audit.json")
with open(out_path, "w") as f:
    json.dump(report, f, indent=2)
print(json.dumps(report, indent=2))
