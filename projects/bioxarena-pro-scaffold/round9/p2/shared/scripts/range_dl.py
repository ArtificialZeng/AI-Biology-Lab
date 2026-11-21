#!/usr/bin/env python3
"""Resumable multi-range HTTP downloader with prefix preservation and MD5/SHA-256 verification."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import socket
import threading
import time
import urllib.request
from pathlib import Path
from typing import List, Set, Tuple

socket.setdefaulttimeout(25)


def compute_hashes(path: Path) -> Tuple[str, str]:
    md5 = hashlib.md5()
    sha = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(4 << 20), b""):
            md5.update(chunk)
            sha.update(chunk)
    return md5.hexdigest(), sha.hexdigest()


def download_range_chunk(
    url: str,
    start: int,
    end: int,
    max_retries: int = 12,
    timeout: int = 90,
) -> bytes:
    import subprocess

    expected = end - start + 1
    for attempt in range(max_retries):
        try:
            res = subprocess.run(
                [
                    "curl", "-sL", "--fail",
                    "--connect-timeout", "15",
                    "--speed-time", "15",
                    "--speed-limit", "30000",
                    "--max-time", str(timeout),
                    "-r", f"{start}-{end}",
                    url,
                ],
                capture_output=True,
                timeout=timeout + 10,
            )
            if res.returncode == 0 and len(res.stdout) == expected:
                return res.stdout
            raise IOError(f"curl rc={res.returncode}, got={len(res.stdout)} != {expected}")
        except Exception as e:
            if attempt == max_retries - 1:
                raise RuntimeError(f"Failed range {start}-{end} after {max_retries} retries: {e}") from e
            time.sleep(min(1.5 * (attempt + 1), 10.0))
    raise RuntimeError("Unreachable")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--size", type=int, required=True)
    ap.add_argument("--md5", default=None)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--chunk-mb", type=int, default=4)
    args = ap.parse_args()

    out: Path = args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    state_path = Path(str(out) + ".dlstate.json")
    total = int(args.size)
    chunk_bytes = int(args.chunk_mb) * 1024 * 1024

    if out.exists() and out.stat().st_size == total and not state_path.exists():
        print(f"[ALREADY COMPLETE] {out} ({total} B)", flush=True)
        got_md5, got_sha = compute_hashes(out)
        if args.md5 and got_md5 != args.md5:
            raise ValueError(f"MD5 mismatch: {got_md5} != {args.md5}")
        print(f"[VERIFIED] md5={got_md5} sha256={got_sha}", flush=True)
        return

    existing_len = out.stat().st_size if out.exists() else 0
    if existing_len > total:
        raise ValueError(f"Existing file {existing_len} > expected {total}")

    done_ranges: Set[Tuple[int, int]] = set()
    covered_prefix = 0
    if state_path.exists():
        st = json.loads(state_path.read_text())
        if st.get("total") == total:
            done_ranges = {(int(a), int(b)) for a, b in st.get("done", [])}
            covered_prefix = int(st.get("prefix", 0))
    elif existing_len > 0:
        safe_prefix = max(0, existing_len - (1 << 20))
        if safe_prefix > 0:
            done_ranges.add((0, safe_prefix - 1))
            covered_prefix = safe_prefix

    todo: List[Tuple[int, int]] = []
    pos = covered_prefix
    while pos < total:
        end = min(pos + chunk_bytes - 1, total - 1)
        if (pos, end) not in done_ranges:
            todo.append((pos, end))
        pos = end + 1

    # Always persist state BEFORE truncating file to full target size
    state_path.write_text(json.dumps({"total": total, "prefix": covered_prefix, "done": sorted(done_ranges)}))

    if not out.exists():
        with out.open("wb") as f:
            f.truncate(total)
    else:
        with out.open("r+b") as f:
            f.truncate(total)

    if todo:
        print(
            f"[START] {out.name}: total={total/1e9:.3f} GB, preserved_prefix={covered_prefix/1e9:.3f} GB, "
            f"todo_chunks={len(todo)} ({args.workers} workers, {args.chunk_mb} MB/chunk)",
            flush=True,
        )
        lock = threading.Lock()
        bytes_done = sum(b - a + 1 for a, b in done_ranges)
        t0 = time.time()
        last_print = t0
        f_handle = out.open("r+b")

        def worker(rng: Tuple[int, int]) -> int:
            s, e = rng
            data = download_range_chunk(args.url, s, e)
            with lock:
                f_handle.seek(s)
                f_handle.write(data)
                done_ranges.add((s, e))
                nonlocal bytes_done, last_print
                bytes_done += len(data)
                now = time.time()
                if (now - last_print >= 5.0) or bytes_done == total or len(done_ranges) % 5 == 0:
                    f_handle.flush()
                    state_path.write_text(
                        json.dumps({"total": total, "prefix": covered_prefix, "done": sorted(done_ranges)})
                    )
                    el = max(now - t0, 1e-6)
                    spd = (bytes_done - covered_prefix) / el / 1e6
                    print(
                        f"  [{out.name}] {bytes_done/1e9:.3f}/{total/1e9:.3f} GB "
                        f"({100.0*bytes_done/total:.1f}%) | speed={spd:.2f} MB/s",
                        flush=True,
                    )
                    last_print = now
            return len(data)

        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
                futs = [ex.submit(worker, rng) for rng in todo]
                for fut in concurrent.futures.as_completed(futs):
                    fut.result()
        finally:
            f_handle.flush()
            f_handle.close()
            state_path.write_text(
                json.dumps({"total": total, "prefix": covered_prefix, "done": sorted(done_ranges)})
            )

    assert out.stat().st_size == total, f"Size mismatch {out.stat().st_size} != {total}"
    got_md5, got_sha = compute_hashes(out)
    if args.md5 and got_md5 != args.md5:
        raise ValueError(f"MD5 mismatch for {out.name}: got {got_md5}, expected {args.md5}")
    if state_path.exists():
        state_path.unlink()
    print(f"[DONE] {out.name}: size={total} B, md5={got_md5}, sha256={got_sha}", flush=True)


if __name__ == "__main__":
    main()
