"""Independent public-API audit: analytic uniform values and direct Decimal formula."""
import argparse
from decimal import Decimal, localcontext
import hashlib
import importlib
import json
import math
from pathlib import Path
import platform
import sys
import warnings

import numpy as np
import scipy
import skbio
import skbio.diversity.alpha._base as base
from skbio.diversity.alpha import hill, renyi


def oracle(counts, q):
    with localcontext() as ctx:
        ctx.prec = 100
        c = [Decimal(x) for x in counts if x > 0]
        total = sum(c)
        p = [v / total for v in c]
        q = Decimal(str(q))
        h = -sum(v * v.ln() for v in p) if q == 1 else sum(v**q for v in p).ln() / (1 - q)
        return float(h.exp()), float(h), float(h / Decimal(2).ln())


def close(got, expected):
    if len(got) != len(expected) or not all(math.isclose(a, b, rel_tol=2e-13, abs_tol=2e-13) for a, b in zip(got, expected)):
        raise AssertionError((got, expected))


parser = argparse.ArgumentParser()
parser.add_argument('--mode', choices=['baseline', 'candidate'], required=True)
parser.add_argument('--source', type=Path, required=True)
parser.add_argument('--expected-hash', required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
module = Path(base.__file__).resolve()
expected = (args.source / 'skbio/diversity/alpha/_base.py').resolve()
if module != expected:
    raise AssertionError(('wrong import provenance', str(module), str(expected)))
source_hash = hashlib.sha256(module.read_bytes()).hexdigest()
if source_hash != args.expected_hash:
    raise AssertionError(('source hash mismatch', source_hash))
extension = importlib.import_module('skbio.diversity._phylogenetic')
ext_path = Path(extension.__file__).resolve()
if not ext_path.is_relative_to(args.source.resolve()) or not ext_path.name.endswith(('.so', '.pyd')):
    raise AssertionError(('wrong native extension provenance', str(ext_path)))
rows = []
for n, q in [(2, 0), (2, 1), (2, 2), (2, 10), (2, math.inf), (2, 2000), (1024, 128)]:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        got = [float(hill([1] * n, order=q)), float(renyi([1] * n, order=q, base=2))]
    known_failure = (n, q) in [(2, 2000), (1024, 128)]
    if args.mode == 'baseline' and known_failure:
        if not all(math.isinf(x) and x > 0 for x in got):
            raise AssertionError(('baseline discrepancy no longer reproduced', n, q, got))
    else:
        close(got, [n, math.log2(n)])
        if caught:
            raise AssertionError([str(w.message) for w in caught])
    rows.append({'richness': n, 'order': str(q), 'observed': [str(x) for x in got], 'analytic_expected': [n, math.log2(n)], 'warnings': [str(w.message) for w in caught]})
ordinary = extreme = 0
if args.mode == 'candidate':
    communities = [[1], [0, 1, 0], [1, 1], [1, 3, 4], [1, 10, 10], [1, 2, 3, 4, 5], [0, 2, 2, 0, 2], [1, 100, 999], list(range(1, 51)), [1] * 1024]
    for c in communities:
        for q in [0, .5, .9, 1, 1.1, 1.9, 2, 2.01, 10, 128, 2000]:
            ref = oracle(c, q)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                with np.errstate(divide='raise', invalid='raise', over='raise'):
                    got = [float(hill(c, order=q)), float(renyi(c, order=q)), float(renyi(c, order=q, base=2))]
            close(got, ref)
            if caught:
                raise AssertionError([str(w.message) for w in caught])
            ordinary += 1
        pmax = max(c) / sum(c)
        limit = -math.log(pmax)
        for q in [1e308, np.finfo(float).max]:
            # For S positive/zero slots, |H_q - H_inf| <=
            # (|log(pmax)| + log(S))/(q-1). At these q this is below
            # the comparison tolerance. This bound is not a direct finite-q oracle.
            bound = (abs(math.log(pmax)) + math.log(len(c))) / (q - 1)
            if bound > 1e-13:
                raise AssertionError(bound)
            with np.errstate(divide='raise', invalid='raise', over='raise'):
                close([float(hill(c, order=q)), float(renyi(c, order=q))], [1 / pmax, limit])
            extreme += 1
    for f in [hill, renyi]:
        if not (math.isnan(f([])) and math.isnan(f([0, 0]))):
            raise AssertionError('empty-community contract changed')
report = {'mode': args.mode, 'all_pass': True, 'python': platform.python_version(), 'executable': sys.executable, 'skbio': skbio.__version__, 'numpy': np.__version__, 'scipy': scipy.__version__, 'module': str(module), 'source_sha256': source_hash, 'native_module': str(ext_path), 'native_sha256': hashlib.sha256(ext_path.read_bytes()).hexdigest(), 'uniform_rows': rows, 'decimal_cases': ordinary, 'extreme_limit_bound_cases': extreme, 'oracle': 'analytic uniform values; direct 100-digit Decimal probability powers; independent limiting bound'}
args.output.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({'mode': args.mode, 'all_pass': True, 'decimal_cases': ordinary, 'extreme_limit_bound_cases': extreme}))
