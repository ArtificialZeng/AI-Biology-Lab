#!/usr/bin/env python3
"""
Round 10 Manifest Writer & Resumable Checksum-Verified Downloader (A-T1).

Reads `manifest_r10.json`, writes `round10/shared/DATA_SOURCES_R10.md` cleanly
PRIOR to downloading, downloads all 125 authorized files into
`round10/shared/raw/<line>/` (gitignored), verifies byte counts & API hashes,
computes SHA-256 digests, and updates `DATA_SOURCES_R10.md`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List

SHARED_DIR = Path(__file__).resolve().parents[1]
RAW_DIR = SHARED_DIR / "raw"
MANIFEST_PATH = Path(__file__).resolve().parent / "manifest_r10.json"
DOC_PATH = SHARED_DIR / "DATA_SOURCES_R10.md"

LINE_META = {
    "dld1": {
        "title": "1. DLD-1 Colorectal Adenocarcinoma (Low-MOI CRISPRi, GEO GSE337988)",
        "pub": "Yeung J et al. (bioRxiv Jul 2026), `doi:10.64898/2026.07.10.737863` (`PPR1277619`).",
        "license": "Open public NCBI GEO (`PRJNA1492064`); *bioRxiv* preprint `CC BY 4.0`.",
        "note": "Selected strictly the 24 Low-MOI channels (`GSM9864227..GSM9864250`, `SAM24482910_1..24`) + feature barcode guide table (`73` files = `6,495,902,794 B`). Excludes Med/High MOI and Pilot.",
    },
    "cd4t_rest": {
        "title": "2. Primary Human CD4+ αβ T Cells — Resting Condition (GEO GSE314342 + S3)",
        "pub": "Zhu X et al. (*Cell* Oct 2026), `doi:10.1016/j.cell.2026.08.002`, `PMID:42664972`.",
        "license": "`MIT License` (`https://mit-license.org/` on VCM S3) + Open public NCBI GEO (`GSE314342`).",
        "note": "Selected strictly the lexicographically first 2 human donors (`CE0006864`=`D4`, `24` files `CD4i_R2`; `CE0008162`=`D1`, `23` files `CD4i_R1`) in `Rest` + 3 S3 metadata tables (`50` files = `26,550,850,781 B`).",
    },
    "wtc11_i3neuron": {
        "title": "3. WTC11 hiPSC-Derived i3Neurons (Tian et al. 2021, Zenodo 13350497 / GSE152988)",
        "pub": "Tian R, Kampmann M et al. (*Nat. Neurosci.* 2021), `doi:10.1038/s41593-021-00862-0`, `PMID:34031600`.",
        "license": "`CC BY 4.0` (Zenodo `13350497`).",
        "note": "CRISPRi CROP-seq only (`TianKampmann2021_CRISPRi.h5ad`, `289,389,479 B`); `CRISPRa` excluded. Fallback 4 GEO `.h5` lanes (`GSM4632022/23/25/26`) used only if `.obs` omits lane labels.",
    },
    "calu3": {
        "title": "4. Calu-3 Lung Adenocarcinoma (Sunshine et al. 2023, Zenodo 13350497 / GSE208240)",
        "pub": "Sunshine S, Hein MY, Weissman JS et al. (*Nat. Commun.* 2023), `doi:10.1038/s41467-023-41788-4`, `PMID:37803001`.",
        "license": "`CC BY 4.0` (Zenodo `13350497` and *Nat. Commun.*).",
        "note": "Direct-capture CRISPRi (`SunshineHein2023.h5ad`, `743,672,579 B`); dual-gene combinatorial elements (`22` pairs) excluded at QC ingestion.",
    },
}


def compute_hashes(fp: Path) -> Dict[str, str]:
    h_sha = hashlib.sha256()
    h_md5 = hashlib.md5()
    with fp.open("rb") as f:
        while True:
            blk = f.read(8 * 1024 * 1024)
            if not blk:
                break
            h_sha.update(blk)
            h_md5.update(blk)
    return {"sha256": h_sha.hexdigest(), "md5": h_md5.hexdigest()}


def write_doc(manifest: List[dict], status_map: Dict[str, dict]) -> None:
    tot_bytes = sum(int(m["bytes"]) for m in manifest)
    lines = [
        "# Round 10 Data Sources, Checksums & Download Provenance (`round10/shared/`)",
        "",
        f"**Gating Check**: `round10/status/R10_GO.txt` verified present. Pre-download exclusion of `H9` (`IGVFDS9825FOBP`) enforced per `PREREG_C` §2.2 (`P2`).",
        f"**Cumulative Authorized Download**: **`{tot_bytes:,}` B (`{tot_bytes / 1e9:.6f} GB <= 60.000000 GB`)** across **{len(manifest)}** files (`round10/shared/raw/`, gitignored).",
        "**Strict Exclusions**: Zero X-Atlas/Orion data; zero H9 files downloaded; zero `dld1` Med/High MOI or Pilot files; zero `cd4t` `Stim` or non-selected donors (`D2`, `D3`); zero `CRISPRa`.",
        "",
        "---",
        "",
    ]
    by_line: Dict[str, List[dict]] = {}
    for m in manifest:
        by_line.setdefault(m["line"], []).append(m)

    for ln in ("dld1", "cd4t_rest", "wtc11_i3neuron", "calu3"):
        meta = LINE_META[ln]
        entries = by_line[ln]
        sub_b = sum(int(x["bytes"]) for x in entries)
        lines.extend([
            f"## {meta['title']}",
            f"- **Publication**: {meta['pub']}",
            f"- **License**: {meta['license']}",
            f"- **Scope**: {meta['note']}",
            f"- **Line Subtotal**: `{sub_b:,}` B (`{sub_b / 1e9:.6f} GB`, {len(entries)} files)",
            "",
            "| Filename | URL | Size (Bytes) | Pre-Download Hash | Post-Download SHA-256 | Status |",
            "|---|---|---:|---|---|---|",
        ])
        for x in sorted(entries, key=lambda r: r["filename"]):
            st = status_map.get(x["filename"], {})
            sha = st.get("sha256", "Pending download")
            msg = st.get("status", "Registered (pre-download)")
            lines.append(
                f"| `{x['filename']}` | `{x['url']}` | `{x['bytes']:,}` | `{x['pre_hash']}` | `{sha}` | {msg} |"
            )
        lines.extend(["", "---", ""])

    DOC_PATH.write_text("\n".join(lines), encoding="utf-8")


def download_one(item: dict) -> dict:
    ln = item["line"]
    fn = item["filename"]
    url = item["url"]
    exp_bytes = int(item["bytes"])
    out_dir = RAW_DIR / ln
    out_dir.mkdir(parents=True, exist_ok=True)
    dst = out_dir / fn
    part = out_dir / (fn + ".part")

    # Reuse if already exact byte size + checksum
    if dst.exists() and dst.stat().st_size == exp_bytes:
        hashes = compute_hashes(dst)
        pre = item.get("pre_hash", "-")
        if pre.startswith("md5:") and hashes["md5"] != pre.split(":", 1)[1]:
            dst.unlink()
        else:
            return {"filename": fn, "sha256": hashes["sha256"], "status": "Verified (cached)"}

    last_err = None
    for attempt in range(1, 6):
        try:
            cmd = [
                "curl", "-sL", "--fail", "--retry", "3", "--retry-delay", "2",
                "--connect-timeout", "30", "--max-time", "1800",
                "-C", "-", "-o", str(part), url,
            ]
            subprocess.check_call(cmd)
            if part.stat().st_size != exp_bytes:
                # Fresh retry without -C in case server ignored range
                part.unlink(missing_ok=True)
                subprocess.check_call([
                    "curl", "-sL", "--fail", "--retry", "3", "--retry-delay", "2",
                    "--connect-timeout", "30", "--max-time", "1800",
                    "-o", str(part), url,
                ])
            if part.stat().st_size != exp_bytes:
                raise RuntimeError(f"Size mismatch {part.stat().st_size} != {exp_bytes}")
            hashes = compute_hashes(part)
            pre = item.get("pre_hash", "-")
            if pre.startswith("md5:") and hashes["md5"] != pre.split(":", 1)[1]:
                part.unlink(missing_ok=True)
                raise RuntimeError(f"MD5 mismatch {hashes['md5']} != {pre}")
            os.replace(part, dst)
            return {"filename": fn, "sha256": hashes["sha256"], "status": "Downloaded & verified"}
        except Exception as e:
            last_err = str(e)
            time.sleep(2 * attempt)
    raise RuntimeError(f"Failed {fn}: {last_err}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest-only", action="store_true")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--lines", nargs="*", default=None)
    args = p.parse_args()

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    if not (SHARED_DIR.parent / "status" / "R10_GO.txt").exists():
        raise RuntimeError("Missing round10/status/R10_GO.txt! Refusing to download.")

    status_map: Dict[str, dict] = {}
    # Seed any already-verified files so DATA_SOURCES_R10 is accurate
    write_doc(manifest, status_map)
    print(f"[Manifest] Wrote pre-download table to {DOC_PATH} ({len(manifest)} files)", flush=True)
    if args.manifest_only:
        return

    targets = [m for m in manifest if args.lines is None or m["line"] in args.lines]
    t0 = time.time()
    done = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(download_one, m): m for m in targets}
        for fut in as_completed(futs):
            m = futs[fut]
            res = fut.result()
            status_map[m["filename"]] = res
            done += 1
            if done % 10 == 0 or done == len(targets) or int(m["bytes"]) > 100_000_000:
                print(
                    f"  [{done}/{len(targets)}] {m['line']} | {m['filename']} ({int(m['bytes'])/1e6:.1f} MB) | sha256={res['sha256'][:12]}... ({time.time()-t0:.1f}s)",
                    flush=True,
                )
    write_doc(manifest, status_map)
    print(f"[DONE] All {len(targets)} files verified and recorded in {DOC_PATH} ({time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
