"""Controller: build FREEZE_MANIFEST_R16_DATA.txt for Round 16 Groups 1-2 (line data + raw control files).

Scope = every input file the frozen Round 16 code (round16/src and round16/audit) reads for Groups 1-2:
  - Nadig raw h5ad (local copy, verified against the GEO sha256 sidecars copied from the sync share)
  - HGNC bridge table (round10 DEFAULT_HGNC_FILE)
  - round9/p2/shared/lines/* (R9 source lines; FCT sources and Group 1 targets)
  - round10/shared/lines/{dld1,cd4t_rest}, round11/shared/lines/feng_ipsc, round12/shared/lines/a549_ko
  - raw control dirs: round10/shared/raw/{dld1,cd4t_rest}, round11/shared/raw/feng_ipsc, round12/shared/raw/a549_ko
Group 3 (GSE295214) is out of scope and receives its own data freeze (PREREG_I section 4.3).

Usage: python3 -I make_data_manifest_r16.py <out_manifest> [--verify <existing_manifest>]
"""
import hashlib
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

P = Path("/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold")
NADIG = Path("/Users/mac/.cache/r9_p2_raw/nadig2025")
EXTRA = [
    NADIG / "GSE264667_hepg2_raw_singlecell_01.h5ad",
    NADIG / "GSE264667_jurkat_raw_singlecell_01.h5ad",
    Path("/Users/mac/.cache/r9_p2_raw/hgnc_complete_set.txt"),
]
DIRS = [
    "round9/p2/shared/lines",
    "round10/shared/lines/dld1", "round10/shared/lines/cd4t_rest",
    "round11/shared/lines/feng_ipsc", "round12/shared/lines/a549_ko",
    "round10/shared/raw/dld1", "round10/shared/raw/cd4t_rest",
    "round11/shared/raw/feng_ipsc", "round12/shared/raw/a549_ko",
]


def files():
    out = []
    for d in DIRS:
        for f in (P / d).rglob("*"):
            if f.is_file() and not f.is_symlink() and f.name != ".DS_Store" and not f.name.startswith("._"):
                out.append(f)
    out += EXTRA
    return sorted(set(out), key=lambda x: str(x))


def sha(f):
    h = hashlib.sha256()
    with open(f, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def label(f):
    try:
        return str(f.relative_to(P))
    except ValueError:
        return str(f)


def main():
    out = Path(sys.argv[1])
    fs = files()
    with ThreadPoolExecutor(8) as ex:
        hs = list(ex.map(sha, fs))
    rows = [f"{h}  {label(f)}  {f.stat().st_size}" for h, f in zip(hs, fs)]
    if len(sys.argv) > 3 and sys.argv[2] == "--verify":
        ref = [l for l in Path(sys.argv[3]).read_text().splitlines() if l and not l.startswith("#")]
        ok = ref == rows
        print("VERIFY", "OK" if ok else "MISMATCH", len(rows), "files")
        if not ok:
            print("\n".join(sorted(set(ref) ^ set(rows))[:40]))
        sys.exit(0 if ok else 1)
    out.write_text("\n".join(rows) + "\n")
    print("WROTE", out, len(rows), "files", sum(f.stat().st_size for f in fs), "bytes")


if __name__ == "__main__":
    main()
