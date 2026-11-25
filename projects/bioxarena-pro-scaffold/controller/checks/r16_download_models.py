"""Controller download of the approved Arc STATE files (explicit list only; pinned revisions; sha256 verified).
Usage: python3 -I r16_download_models.py <plan.json> <out_root> <checksums_out>
"""
import hashlib, json, subprocess, sys
from pathlib import Path
plan = json.load(open(sys.argv[1])); out = Path(sys.argv[2]); ck = Path(sys.argv[3])
lines, bad = [], []
for repo, rev, name, size, lfs_sha in plan:
    assert '/eval_final.ckpt/' not in name and not name.endswith(('metrics.csv', 'last.ckpt')), name
    dst = out / repo.split('/')[1] / name
    dst.parent.mkdir(parents=True, exist_ok=True)
    if not (dst.exists() and dst.stat().st_size == size):
        url = f"https://huggingface.co/{repo}/resolve/{rev}/{name}"
        for attempt in range(30):  # resume partial files across transient SSL/network errors
            r = subprocess.run(["curl", "-sSLf", "-C", "-", "--retry", "8", "--retry-delay", "5", "--retry-all-errors", "-o", str(dst), url])
            if r.returncode == 0 and dst.stat().st_size == size:
                break
            if dst.exists() and dst.stat().st_size > size:
                dst.unlink()
        else:
            raise RuntimeError("download failed after 30 attempts: " + name)
    h = hashlib.sha256()
    with open(dst, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    d = h.hexdigest(); ok_size = dst.stat().st_size == size; ok_sha = (lfs_sha is None) or (d == lfs_sha)
    if not (ok_size and ok_sha):
        bad.append(name)
    lines.append(f"{d}  {dst.relative_to(out.parent.parent.parent)}")
    print(("OK " if ok_size and ok_sha else "BAD ") + repo + ":" + name, size, flush=True)
ck.write_text("\n".join(lines) + "\n")
print("DONE files=%d bad=%d" % (len(plan), len(bad)), bad)
