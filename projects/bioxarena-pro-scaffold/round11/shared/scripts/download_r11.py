#!/usr/bin/env python3
"""
Round 11 Manifest Writer & Resumable Checksum-Verified Downloader (A-T1).

Strictly gated on `round11/status/R11_GO.txt`. Reads `manifest_r11.json`,
writes `round11/shared/DATA_SOURCES_R11.md` pre-download, downloads all
approved files into `round11/shared/raw/<line>/` (gitignored, total <= 80 GB),
verifies exact byte sizes and MD5s (plus structural HDF5 readability for GEO .h5),
and records SHA-256 checksums or exact failure diagnostics without substituting files.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List

SHARED_DIR = Path(__file__).resolve().parents[1]
RAW_DIR = SHARED_DIR / "raw"
MANIFEST_PATH = Path(__file__).resolve().parent / "manifest_r11.json"
DOC_PATH = SHARED_DIR / "DATA_SOURCES_R11.md"
GO_PATH = SHARED_DIR.parent / "status" / "R11_GO.txt"

LINE_META = {
    "telohaec": {
        "title": "1. TeloHAEC_2285 Arterial Endothelial (IGVF IGVFDS7943IIWZ / GEO GSE210681)",
        "pub": "Schnitzler GR, Engreitz JM et al. (*Nature* 2024), `PMC10921916` (`GSE210681` / `IGVFDS7943IIWZ`).",
        "license": "Open public IGVF (`IGVFFI8645YMGX.h5ad`, `IGVFFI9754AGFB.tsv.gz`) + NCBI GEO (`PRJNA867016`).",
        "note": "Confirmatory Candidate 1 (resting DOX-inducible KRAB-dCas9; new endothelial lineage & Engreitz lab). `3` files = `1,486,778,253` B (`1.487 GB`).",
    },
    "hcasmc": {
        "title": "2. HCASMC-hTERT Coronary Artery Smooth Muscle (IGVF IGVFDS6248NDUW, Parse v1.2.1)",
        "pub": "Quertermous Lab / IGVF Consortium (`IGVFDS6248NDUW`, Parse Biosciences `v1.2.1`, MOI `0.1-0.2`, ZIM3-KRAB).",
        "license": "Open public IGVF (`IGVFFI5999NXAF.h5`, `IGVFFI4575UMXX.tsv.gz`, `IGVFFI9546UECF.tsv.gz`).",
        "note": "Confirmatory Candidate 2 (smooth muscle lineage; `IGVFFI9007DBZO.h5ad` `87.29 GB` & `IGVFFI7271ZXXZ.h5` pilot excluded; subject to `Q3` pre-eval streaming & guide feasibility gate). `3` files = `42,885,287,884` B (`42.885 GB`).",
    },
    "cd4t_stim48": {
        "title": "3. Primary Human CD4+ alpha-beta T Cells — Stim48hr (`GSE314342` `CD4i_R2` Donors `CE0008678` & `CE0010866`)",
        "pub": "Zhu X et al. (*Cell* Oct 2026), `doi:10.1016/j.cell.2026.08.002`, `PMID:42664972`.",
        "license": "`MIT License` (VCM S3 metadata CSVs) + Open public NCBI GEO (`GSE314342`).",
        "note": "Confirmatory Candidate 3 (`48` `Stim48hr` `.h5` lanes `GSM9393875..GSM9393898`, `CD4i_R2L25..R2L48` + `3` verified shared S3 CSVs reused locally = `27,397,331,332` B / `27.397 GB`; donor/state disjoint from R10 `cd4t_rest`, weaker lab/lineage independence explicitly disclosed per `Q2`/`Q7`).",
    },
    "feng_ipsc": {
        "title": "4. Feng et al. 2026 HipSci iPSC TargetedScreen (Figshare `27989294 v2`, excl. `kolf_2` & `kolf_3`)",
        "pub": "Feng E et al. (*Nature Communications* / `PMC12903452`), `doi:10.6084/m9.figshare.27989294.v2`.",
        "license": "`MIT License` (Figshare `27989294 v2`).",
        "note": "Mandatory Secondary Only (`H-D5` / `SA-KOLF`; never enters `K_conf`). `17` non-KOLF HipSci lines (`9` donors) after excluding `kolf_2`/`kolf_3`. `3` files (incl. optional `Guide-UMI-Counts.csv.gz`) = `4,568,239,476` B (`4.568 GB`).",
    },
}

_DOC_LOCK = threading.Lock()


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
        "# Round 11 Data Sources, Checksums & Download Provenance (`round11/shared/`)",
        "",
        f"**Gating Check**: `round11/status/R11_GO.txt` verified present (`PREREG_D` commit `72fc0c9`).",
        f"**Cumulative Approved Manifest**: **`{tot_bytes:,}` B (`{tot_bytes / 1e9:.6f} GB <= 80.000000 GB`)** across **{len(manifest)}** entries (`round11/shared/raw/`, gitignored).",
        "**Strict Exclusions**: Zero `X-Atlas/Orion` files; zero `Cas9 KO` datasets; zero `IGVFFI9007DBZO.h5ad` (`87.29 GB`) or `IGVFFI7271ZXXZ.h5` pilot; zero R10 `cd4t_rest`/`CE0006864`/`CE0008162`/`Stim8hr`; zero `Fitness`/`NonFitness` iPSC screens; strict exclusion of any non-manifest fallback replacement.",
        "",
        "---",
        "",
    ]
    by_line: Dict[str, List[dict]] = {}
    for m in manifest:
        by_line.setdefault(m["line"], []).append(m)

    for ln in ("telohaec", "hcasmc", "cd4t_stim48", "feng_ipsc"):
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
            key = f"{x['line']}/{x['filename']}"
            st = status_map.get(key, {})
            sha = st.get("sha256", "Pending download")
            msg = st.get("status", "Registered (pre-download)")
            lines.append(
                f"| `{x['filename']}` | `{x['url']}` | `{x['bytes']:,}` | `{x['pre_hash']}` | `{sha}` | {msg} |"
            )
        lines.extend(["", "---", ""])

    with _DOC_LOCK:
        DOC_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _verify_geo_h5_header(fp: Path) -> None:
    import h5py

    with h5py.File(fp, "r") as f:
        if "matrix" not in f:
            raise RuntimeError(f"Missing /matrix group in {fp.name}")
        m = f["matrix"]
        for req in ("shape", "barcodes", "indptr", "indices", "data", "features"):
            if req not in m:
                raise RuntimeError(f"Missing /matrix/{req} in {fp.name}")


def download_one(item: dict) -> dict:
    ln = item["line"]
    fn = item["filename"]
    url = item["url"]
    exp_bytes = int(item["bytes"])
    out_dir = RAW_DIR / ln
    out_dir.mkdir(parents=True, exist_ok=True)
    dst = out_dir / fn
    part = out_dir / (fn + ".part")

    # Local verified cache copy for shared R10 CD4i CSVs
    if url.startswith("cached://"):
        src = Path(url[len("cached://") :])
        if not src.exists():
            return {"sha256": "-", "status": f"FAILED (missing cached {src})"}
        if not (dst.exists() and dst.stat().st_size == exp_bytes):
            shutil.copy2(src, part)
            if part.stat().st_size != exp_bytes:
                part.unlink(missing_ok=True)
                return {"sha256": "-", "status": f"FAILED (cached size mismatch)"}
            os.replace(part, dst)
        hashes = compute_hashes(dst)
        pre = item.get("pre_hash", "-")
        if pre.startswith("md5:") and hashes["md5"] != pre.split(":", 1)[1]:
            dst.unlink(missing_ok=True)
            return {"sha256": "-", "status": f"FAILED (MD5 mismatch {hashes['md5']} != {pre})"}
        return {"sha256": hashes["sha256"], "status": "Verified (copied from verified R10 CSV)"}

    # Reuse if already exact byte size + checksum
    if dst.exists() and dst.stat().st_size == exp_bytes:
        hashes = compute_hashes(dst)
        pre = item.get("pre_hash", "-")
        if pre.startswith("md5:") and hashes["md5"] != pre.split(":", 1)[1]:
            dst.unlink()
        else:
            if ln == "cd4t_stim48" and fn.endswith(".h5"):
                _verify_geo_h5_header(dst)
            return {"sha256": hashes["sha256"], "status": "Verified (cached)"}

    last_err = None
    max_t = "21600" if exp_bytes > 5_000_000_000 else "3600"
    for attempt in range(1, 6):
        try:
            cmd = [
                "curl", "-sL", "--fail", "--retry", "3", "--retry-delay", "5",
                "--connect-timeout", "30", "--speed-limit", "20000", "--speed-time", "120",
                "--max-time", max_t, "-C", "-", "-o", str(part), url,
            ]
            subprocess.check_call(cmd)
            if part.stat().st_size != exp_bytes:
                part.unlink(missing_ok=True)
                subprocess.check_call([
                    "curl", "-sL", "--fail", "--retry", "3", "--retry-delay", "5",
                    "--connect-timeout", "30", "--speed-limit", "20000", "--speed-time", "120",
                    "--max-time", max_t, "-o", str(part), url,
                ])
            if part.stat().st_size != exp_bytes:
                raise RuntimeError(f"Size mismatch {part.stat().st_size} != {exp_bytes}")
            hashes = compute_hashes(part)
            pre = item.get("pre_hash", "-")
            if pre.startswith("md5:") and hashes["md5"] != pre.split(":", 1)[1]:
                part.unlink(missing_ok=True)
                raise RuntimeError(f"MD5 mismatch {hashes['md5']} != {pre}")
            if ln == "cd4t_stim48" and fn.endswith(".h5"):
                _verify_geo_h5_header(part)
            os.replace(part, dst)
            st_lbl = "Downloaded & verified (MD5+SHA256)" if pre.startswith("md5:") else "Downloaded & verified (Bytes+H5+SHA256)"
            return {"sha256": hashes["sha256"], "status": st_lbl}
        except Exception as e:
            last_err = str(e)
            time.sleep(3 * attempt)
    return {"sha256": "-", "status": f"FAILED ({last_err})"}


def main() -> None:
    p = argparse.ArgumentParser(description="Download & verify locked R11 datasets.")
    p.add_argument("--manifest-only", action="store_true")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--lines", nargs="*", default=None)
    args = p.parse_args()

    if not GO_PATH.exists():
        raise RuntimeError(f"Missing {GO_PATH}! Refusing to download before R11_GO.txt.")

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    status_map: Dict[str, dict] = {}
    write_doc(manifest, status_map)
    print(f"[Manifest] Wrote pre-download table to {DOC_PATH} ({len(manifest)} entries)", flush=True)
    if args.manifest_only:
        return

    targets = [m for m in manifest if args.lines is None or m["line"] in args.lines]
    # Tiny metadata first, then start 42.89 GB hcasmc immediately alongside telohaec/feng/cd4t
    targets.sort(key=lambda r: (0 if int(r["bytes"]) < 100_000_000 else 1, -int(r["bytes"]) if r["line"] != "cd4t_stim48" else 0, r["filename"]))

    t0 = time.time()
    done = 0
    failed = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(download_one, m): m for m in targets}
        for fut in as_completed(futs):
            m = futs[fut]
            key = f"{m['line']}/{m['filename']}"
            res = fut.result()
            status_map[key] = res
            done += 1
            if res["status"].startswith("FAILED"):
                failed.append((key, res["status"]))
            write_doc(manifest, status_map)
            if res["status"].startswith("FAILED") or done % 5 == 0 or done == len(targets) or int(m["bytes"]) > 100_000_000 or m["line"] != "cd4t_stim48":
                print(
                    f"  [{done}/{len(targets)}] {key} ({int(m['bytes'])/1e6:.2f} MB) | {res['status']} | sha256={res['sha256'][:12]} ({time.time()-t0:.1f}s)",
                    flush=True,
                )
    if failed:
        raise RuntimeError(f"{len(failed)} files failed verification: {failed}")
    print(f"[DONE] All {len(targets)} files verified and recorded in {DOC_PATH} ({time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
