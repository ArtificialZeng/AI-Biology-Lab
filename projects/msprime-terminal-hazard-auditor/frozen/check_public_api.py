"""Independent scalar/isolated-population oracles against a selected source tree."""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path
import sys
import warnings

parser = argparse.ArgumentParser()
parser.add_argument("--source", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
source = args.source.resolve()
sys.path.insert(0, str(source))
import numpy as np
import msprime

assert Path(msprime.__file__).resolve().is_relative_to(source)
assert Path(msprime._msprime.__file__).resolve().is_relative_to(source)
assert Path(msprime.demography.__file__).resolve() == source / "msprime/demography.py"
expected = json.loads((Path(__file__).parent / "expected_sources.json").read_text())[source.name]
assert hashlib.sha256((source / "msprime/demography.py").read_bytes()).hexdigest() == expected["demography_sha256"]
assert hashlib.sha256(Path(msprime._msprime.__file__).read_bytes()).hexdigest() == expected["native_sha256"]

def scalar(t, sizes, clamp=1):
    # Integrate the mathematical piecewise-constant scalar hazard exactly.
    accumulated = 0.0
    current = sizes[0][1]
    for i, (start, size) in enumerate(sizes):
        if t < start:
            break
        current = size
        end = sizes[i + 1][0] if i + 1 < len(sizes) else math.inf
        accumulated += max(0, min(t, end) - start) / (2 * max(clamp, size))
    return 1 / (2 * max(clamp, current)), math.exp(-accumulated)

def one_pop(sizes):
    d = msprime.Demography()
    d.add_population(name="A", initial_size=sizes[0][1])
    for time, size in sizes[1:]:
        d.add_population_parameters_change(time=time, initial_size=size, population="A")
    return d

records = []
def check(name, d, lineages, steps, oracle, validation, clamp=1):
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        rates, survival = d.debug().coalescence_rate_trajectory(
            steps, lineages, min_pop_size=clamp, double_step_validation=validation
        )
    expected = [oracle(t) for t in steps]
    er = np.array([r for r, p in expected])
    ep = np.array([p for r, p in expected])
    # Zero floating-point survival has an explicit NaN-rate convention.
    er = np.where(ep == 0, np.nan, er)
    ok_r = np.allclose(rates, er, rtol=2e-12, atol=2e-15, equal_nan=True)
    ok_p = np.allclose(survival, ep, rtol=2e-12, atol=5e-11)
    records.append(dict(name=name, validation=validation, steps=list(steps),
        rate=rates.tolist(), survival=survival.tolist(), expected_rate=er.tolist(),
        expected_survival=ep.tolist(), rate_ok=bool(ok_r), survival_ok=bool(ok_p),
        warnings=[str(w.message) for w in ws]))

for validation in [False, True]:
    for sizes in [[(0, 1000), (100, 10000)], [(0, 10000), (100, 1000)],
                  [(0, 1000), (30, 2000), (100, 7000)]]:
        for steps in [[0, 100], [100], [0, 99, 100], [0, 100, 101]]:
            check("scalar-size-history", one_pop(sizes), {"A": 2}, steps,
                  lambda t, sizes=sizes: scalar(t, sizes), validation)
    for steps in [[0], [100], [0, 100], [0, 100, 101]]:
        check("constant-control", one_pop([(0, 1234)]), {"A": 2}, steps,
              lambda t: scalar(t, [(0, 1234)]), validation)
    sizes = [(0, 10), (100, 0.25)]
    for steps in [[0, 100], [0, 100, 101]]:
        check("minimum-size-clamp", one_pop(sizes), {"A": 2}, steps,
              lambda t: scalar(t, sizes, 2), validation, clamp=2)
    # A future epoch boundary means requested100 is already an interior point.
    sizes = [(0, 1000), (100, 10000), (200, 100000)]
    check("future-epoch-control", one_pop(sizes), {"A": 2}, [0, 100],
          lambda t: scalar(t, sizes), validation)
    for prob in [0.25, 0.5, 1.0]:
        d = msprime.Demography()
        d.add_population(name="A", initial_size=1000)
        d.add_population(name="B", initial_size=10000)
        d.add_mass_migration(time=100, source="A", dest="B", proportion=prob)
        def movement_oracle(t, prob=prob):
            if t < 100:
                return scalar(t, [(0, 1000)])
            dt = t - 100
            a = (1 - prob) ** 2 * math.exp(-dt / 2000)
            b = prob ** 2 * math.exp(-dt / 20000)
            cross = 2 * prob * (1 - prob)
            p = math.exp(-100 / 2000) * (a + b + cross)
            return (a / 2000 + b / 20000) / (a + b + cross), p
        for steps in [[0, 100], [100], [0, 100, 101]]:
            check("positive-time-movement", d, {"A": 2}, steps,
                  movement_oracle, validation)
    d = msprime.Demography()
    for name, size in [("A", 1000), ("B", 10000), ("C", 100000)]:
        d.add_population(name=name, initial_size=size)
    d.add_mass_migration(time=100, source="A", dest="B", proportion=0.5)
    d.add_mass_migration(time=100, source="B", dest="C", proportion=0.5)
    terminal_hazard = 0.5**2 / 2000 + 0.25**2 / 20000 + 0.25**2 / 200000
    check("ordered-two-movements", d, {"A": 2}, [0, 100],
          lambda t: (terminal_hazard, math.exp(-0.05)) if t == 100 else (0.0005, 1),
          validation)
    d = msprime.Demography()
    d.add_population(name="A", initial_size=1000)
    d.add_population(name="B", initial_size=10000)
    d.add_population_split(time=100, derived=["A"], ancestral="B")
    check("positive-time-population-split", d, {"A": 2}, [0, 100],
          lambda t: scalar(t, [(0, 1000), (100, 10000)]), validation)
    d = msprime.Demography()
    d.add_population(name="A", initial_size=1000)
    d.add_population(name="B", initial_size=10000)
    d.add_mass_migration(time=100, source="A", dest="B", proportion=1)
    check("cross-population-pair-merger", d, {"A": 1, "B": 1}, [0, 100],
          lambda t: (1 / 20000, 1) if t == 100 else (0, 1), validation)
    d.add_population_parameters_change(time=100, population="B", initial_size=20000)
    d.sort_events()
    check("size-change-and-movement", d, {"A": 2}, [0, 100],
          lambda t: (1 / 40000, math.exp(-0.05)) if t == 100 else (0.0005, 1),
          validation)
    d = one_pop([(0, 1000), (100, 10000)])
    d.add_population(name="B", initial_size=2000)
    def mixture_oracle(t):
        ra, pa = scalar(t, [(0, 1000), (100, 10000)])
        rb, pb = scalar(t, [(0, 2000)])
        p = (pa + pb) / 6 + 2 / 3
        return (ra * pa + rb * pb) / (6 * p), p
    check("isolated-population-mixture", d, {"A": 2, "B": 2}, [0, 100],
          mixture_oracle, validation)
    check("floating-zero-survival", one_pop([(0, 1)]), {"A": 2}, [0, 2000],
          lambda t: scalar(t, [(0, 1)]), validation)

# Boundary observations are reported, not retroactively declared supported APIs.
boundary = []
for steps, lineages in [([], {"A": 2}), ([0], {"A": 2}), ([100], {"A": 2}),
                         ([0, 100], {"A": 1}), ([0, 100], {})]:
    try:
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("always")
            r, p = one_pop([(0, 1000)]).debug().coalescence_rate_trajectory(steps, lineages)
        boundary.append(dict(steps=steps, lineages=lineages, rates=r.tolist(),
                             survival=p.tolist(), warnings=[str(w.message) for w in ws]))
    except Exception as error:
        boundary.append(dict(steps=steps, lineages=lineages,
                             error=type(error).__name__, message=str(error)))

payload = dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    module=msprime.__file__, native=msprime._msprime.__file__,
    source_sha256=hashlib.sha256((source / "msprime/demography.py").read_bytes()).hexdigest(),
    records=records, boundary_observations=boundary,
    cases=len(records), failures=sum(not (r["rate_ok"] and r["survival_ok"]) for r in records))
args.output.write_text(json.dumps(payload, indent=2))
print(json.dumps({"source":str(source),"cases":payload["cases"],"failures":payload["failures"]}))
