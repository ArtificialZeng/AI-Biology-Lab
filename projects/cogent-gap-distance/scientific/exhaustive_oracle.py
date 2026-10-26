"""Exact direct-column oracle; finite coverage is stated in the JSON result."""
import hashlib
import itertools
import json
import math
import pathlib
import sys
from fractions import Fraction

import cogent3
import cogent3.evolve.fast_distance as distance

alphabet = "AC-"
rows = ["".join(x) for x in itertools.product(alphabet, repeat=2)]
checks = []
for method in ("pdist", "hamming"):
    errors = []
    checked = 0
    for ordered in itertools.permutations(rows, 3):
        seqs = dict(zip(("x", "y", "z"), ordered, strict=True))
        aln = cogent3.make_aligned_seqs(seqs, moltype="dna")
        calc = distance.get_distance_calculator(method, moltype="dna")
        calc(aln)
        got = calc.get_pairwise_distances()
        for a, b in itertools.combinations(seqs, 2):
            valid = [(s, t) for s, t in zip(seqs[a], seqs[b], strict=True) if s in "ACGT" and t in "ACGT"]
            num_valid = len(valid)
            differences = sum(s != t for s, t in valid)
            expected = None if not num_valid else Fraction(differences, num_valid) if method == "pdist" else Fraction(differences)
            actual = float(got[a, b])
            checked += 1
            ok = math.isnan(actual) if expected is None else actual == float(expected)
            if not ok:
                errors.append({"sequences": seqs, "pair": [a, b], "exact": str(expected), "actual": actual})
    checks.append({"method": method, "ordered_triples": math.perm(len(rows), 3), "checked_pair_values": checked, "mismatch_count": len(errors), "first_mismatches": errors[:6]})
source = pathlib.Path(distance.__file__)
print(json.dumps({"version": cogent3.__version__, "source_module": str(source), "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(), "python": sys.version, "scope": "All ordered triples of distinct length-two rows over A,C,-; pairwise deletion; no statement about untested domains", "checks": checks}, indent=2))
