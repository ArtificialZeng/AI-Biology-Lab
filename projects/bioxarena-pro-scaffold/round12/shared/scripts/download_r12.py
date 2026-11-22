#!/usr/bin/env python3
"""
Round 12 Locked Data Downloader & Checksum Verifier (`PREREG_E`).

Downloads ONLY the 5 manifest-locked files from `round12/shared/R12_DOWNLOAD_PLAN.md`
into `round12/shared/raw/{a549_ko,hek293_xucao}/` after confirming
`round12/status/PREREG_E_COMMIT.txt` and `round12/status/R12_GO.txt`.

Verifies exact byte sizes, `gzip -t` for `.gz` archives, published Zenodo MD5 +
HDF5 readability for `XuCao2023.h5ad`, and records SHA-256 digests into
`round12/shared/R12_CHECKSUMS.sha256` and `round12/shared/download_manifest_r12.json`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional

import h5py

SCRIPTS_DIR = Path(__file__).resolve().parent
SHARED_DIR = SCRIPTS_DIR.parent
ROUND12_ROOT = SHARED_DIR.parent
RAW_DIR = SHARED_DIR / "raw"
STATUS_DIR = ROUND12_ROOT / "status"
PREREG_COMMIT_FILE = STATUS_DIR / "PREREG_E_COMMIT.txt"
R12_GO_FILE = STATUS_DIR / "R12_GO.txt"
CHECKSUM_FILE = SHARED_DIR / "R12_CHECKSUMS.sha256"
MANIFEST_JSON_FILE = SHARED_DIR / "download_manifest_r12.json"

MAX_CUMULATIVE_BYTES = 15_000_000_000
CHUNK_SIZE = 8 * 1024 * 1024


@dataclass(frozen=True)
class ManifestEntry:
    line: str
    rel_path: str
    url: str
    expected_bytes: int
    expected_md5: Optional[str]
    kind: str  # "gz" or "h5ad"
    source_desc: str


# Order smallest -> largest so metadata/h5ad fail fast before the 10.96 GB MTX transfer.
MANIFEST: List[ManifestEntry] = [
    ManifestEntry(
        line="a549_ko",
        rel_path="a549_ko/GSE345058_sgRNA.csv.gz",
        url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345058/suppl/GSE345058_sgRNA.csv.gz",
        expected_bytes=68_194,
        expected_md5=None,
        kind="gz",
        source_desc="NCBI GEO GSE345058 (A549 Cas9 KO) sgRNA table",
    ),
    ManifestEntry(
        line="a549_ko",
        rel_path="a549_ko/GSE345058_gene_annotations.csv.gz",
        url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345058/suppl/GSE345058_gene_annotations.csv.gz",
        expected_bytes=378_455,
        expected_md5=None,
        kind="gz",
        source_desc="NCBI GEO GSE345058 (A549 Cas9 KO) gene annotations",
    ),
    ManifestEntry(
        line="a549_ko",
        rel_path="a549_ko/GSE345058_cell_annotations.csv.gz",
        url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345058/suppl/GSE345058_cell_annotations.csv.gz",
        expected_bytes=40_152_708,
        expected_md5=None,
        kind="gz",
        source_desc="NCBI GEO GSE345058 (A549 Cas9 KO) cell annotations",
    ),
    ManifestEntry(
        line="hek293_xucao",
        rel_path="hek293_xucao/XuCao2023.h5ad",
        url="https://zenodo.org/api/records/13350497/files/XuCao2023.h5ad/content",
        expected_bytes=336_637_209,
        expected_md5="3efb25a2927c5c228f04f783f18004d8",
        kind="h5ad",
        source_desc="Zenodo 13350497 XuCao2023.h5ad (HEK293 dCas9-KRAB-MeCP2)",
    ),
    ManifestEntry(
        line="a549_ko",
        rel_path="a549_ko/GSE345058_raw_counts.mtx.gz",
        url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345058/suppl/GSE345058_raw_counts.mtx.gz",
        expected_bytes=10_956_707_618,
        expected_md5=None,
        kind="gz",
        source_desc="NCBI GEO GSE345058 (A549 Cas9 KO) raw counts MTX",
    ),
]


def assert_gate_unlocked() -> None:
    if not PREREG_COMMIT_FILE.exists() or not PREREG_COMMIT_FILE.read_text(encoding="utf-8").strip():
        raise RuntimeError(f"Blocked: missing or empty {PREREG_COMMIT_FILE}")
    if not R12_GO_FILE.exists() or not R12_GO_FILE.read_text(encoding="utf-8").strip():
        raise RuntimeError(f"Blocked: missing or empty {R12_GO_FILE}")
    total_bytes = sum(e.expected_bytes for e in MANIFEST)
    if total_bytes > MAX_CUMULATIVE_BYTES:
        raise RuntimeError(f"Manifest exceeds 15 GB cap: {total_bytes} > {MAX_CUMULATIVE_BYTES}")


def probe_remote_content_length(url: str) -> Optional[int]:
    proc = subprocess.run(
        ["curl", "-4", "-fsSIL", "--connect-timeout", "20", url],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        return None
    last_cl: Optional[int] = None
    for line in proc.stdout.splitlines():
        low = line.strip().lower()
        if low.startswith("content-length:"):
            val = low.split(":", 1)[1].strip()
            if val.isdigit():
                last_cl = int(val)
    return last_cl


def hash_file(path: Path) -> Dict[str, str]:
    h_sha256 = hashlib.sha256()
    h_md5 = hashlib.md5()
    with path.open("rb") as fh:
        while True:
            buf = fh.read(CHUNK_SIZE)
            if not buf:
                break
            h_sha256.update(buf)
            h_md5.update(buf)
    return {"sha256": h_sha256.hexdigest(), "md5": h_md5.hexdigest()}


def verify_file_integrity(entry: ManifestEntry, dest: Path) -> Dict[str, object]:
    if not dest.exists():
        raise FileNotFoundError(f"Missing downloaded file: {dest}")
    actual_bytes = dest.stat().st_size
    if actual_bytes != entry.expected_bytes:
        raise ValueError(
            f"Size mismatch for {entry.rel_path}: got {actual_bytes}, expected {entry.expected_bytes}"
        )
    if entry.kind == "gz":
        t0 = time.time()
        proc = subprocess.run(["gzip", "-t", str(dest)], capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            raise ValueError(f"gzip -t failed for {entry.rel_path}: {proc.stderr.strip()}")
        extra = f"gzip -t OK ({time.time() - t0:.1f}s)"
    elif entry.kind == "h5ad":
        with h5py.File(dest, "r") as f:
            top_keys = sorted(list(f.keys()))
            for req in ("X", "obs", "var"):
                if req not in f:
                    raise ValueError(f"h5ad missing /{req}: keys={top_keys}")
        extra = f"h5py keys={top_keys}"
    else:
        raise ValueError(f"Unsupported kind {entry.kind}")

    digests = hash_file(dest)
    if entry.expected_md5 and digests["md5"] != entry.expected_md5:
        raise ValueError(
            f"MD5 mismatch for {entry.rel_path}: got {digests['md5']}, expected {entry.expected_md5}"
        )
    return {
        "rel_path": entry.rel_path,
        "line": entry.line,
        "url": entry.url,
        "expected_bytes": entry.expected_bytes,
        "actual_bytes": actual_bytes,
        "expected_md5": entry.expected_md5,
        "md5": digests["md5"],
        "sha256": digests["sha256"],
        "integrity": extra,
    }


def download_entry(entry: ManifestEntry, verify_only: bool = False) -> Dict[str, object]:
    dest = RAW_DIR / entry.rel_path
    dest.parent.mkdir(parents=True, exist_ok=True)

    remote_cl = probe_remote_content_length(entry.url)
    if remote_cl is not None and remote_cl != entry.expected_bytes:
        raise RuntimeError(
            f"Remote Content-Length mismatch for {entry.rel_path}: remote={remote_cl} != locked={entry.expected_bytes}"
        )

    if verify_only:
        return verify_file_integrity(entry, dest)

    cur_size = dest.stat().st_size if dest.exists() else 0
    if cur_size == entry.expected_bytes:
        print(f"[SKIP] {entry.rel_path} already exact size ({cur_size:,} B); verifying...", flush=True)
        try:
            return verify_file_integrity(entry, dest)
        except Exception as exc:
            print(f"[RETRY] Existing {entry.rel_path} failed check ({exc}); removing and re-downloading.", flush=True)
            dest.unlink()
            cur_size = 0
    elif cur_size > entry.expected_bytes:
        print(f"[RESET] Oversized partial {entry.rel_path} ({cur_size:,} > {entry.expected_bytes:,}); removing.", flush=True)
        dest.unlink()
        cur_size = 0

    for attempt in range(1, 5):
        cmd = [
            "curl",
            "-4",
            "-fSL",
            "--connect-timeout",
            "30",
            "--retry",
            "5",
            "--retry-delay",
            "5",
        ]
        if dest.exists() and 0 < dest.stat().st_size < entry.expected_bytes:
            print(
                f"[RESUME attempt {attempt}] {entry.rel_path} from {dest.stat().st_size:,} / {entry.expected_bytes:,} B",
                flush=True,
            )
            cmd.extend(["-C", "-"])
        else:
            print(f"[DOWNLOAD attempt {attempt}] {entry.rel_path} ({entry.expected_bytes:,} B)", flush=True)
        cmd.extend(["-o", str(dest), entry.url])
        t0 = time.time()
        proc = subprocess.run(cmd, check=False)
        elapsed = time.time() - t0
        if proc.returncode == 0 and dest.exists() and dest.stat().st_size == entry.expected_bytes:
            print(
                f"[FETCHED] {entry.rel_path} ({entry.expected_bytes:,} B) in {elapsed:.1f}s; verifying integrity...",
                flush=True,
            )
            return verify_file_integrity(entry, dest)
        got = dest.stat().st_size if dest.exists() else 0
        print(
            f"[WARN attempt {attempt}] rc={proc.returncode}, size={got:,}/{entry.expected_bytes:,} in {elapsed:.1f}s",
            flush=True,
        )
        if got > entry.expected_bytes and dest.exists():
            dest.unlink()
        time.sleep(3 * attempt)

    raise RuntimeError(f"Failed to download exact file after 4 attempts: {entry.rel_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Round 12 locked downloader")
    parser.add_argument("--verify-only", action="store_true", help="Only verify existing files & rewrite checksums")
    args = parser.parse_args()

    assert_gate_unlocked()
    total_expected = sum(e.expected_bytes for e in MANIFEST)
    print(
        f"[GATE OK] Locked manifest: {len(MANIFEST)} files, {total_expected:,} B ({total_expected / 1e9:.6f} GB <= 15.0 GB)",
        flush=True,
    )

    records: List[Dict[str, object]] = []
    for entry in MANIFEST:
        rec = download_entry(entry, verify_only=args.verify_only)
        records.append(rec)
        print(f"[VERIFIED] {rec['rel_path']} | bytes={rec['actual_bytes']:,} | sha256={rec['sha256']} | md5={rec['md5']}", flush=True)

    # Sort checksum manifest by rel_path for deterministic output.
    records_sorted = sorted(records, key=lambda r: str(r["rel_path"]))
    checksum_lines = [f"{r['sha256']}  shared/raw/{r['rel_path']}\n" for r in records_sorted]
    CHECKSUM_FILE.write_text("".join(checksum_lines), encoding="utf-8")

    manifest_out = {
        "total_files": len(records_sorted),
        "total_bytes": sum(int(r["actual_bytes"]) for r in records_sorted),
        "files": records_sorted,
    }
    MANIFEST_JSON_FILE.write_text(json.dumps(manifest_out, indent=2) + "\n", encoding="utf-8")
    print(f"[DONE] Wrote {CHECKSUM_FILE} and {MANIFEST_JSON_FILE}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
