"""Read-only exact finite checks; the infinite analytic proof is in proof.md."""
import sys
sys.dont_write_bytecode = True

from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib
import json
import subprocess

import types

cm = None

HERE = Path(__file__).resolve().parent
MEMBERS = {
    "README.md", "proof.md", "sources.md", "count_moran.py", "verify.py",
    "examples.json", "certificate.json", "LICENSE",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def archive_check():
    manifest_path = HERE / "archive_manifest.json"
    require(manifest_path.is_file() and not manifest_path.is_symlink(), "manifest file kind")
    manifest = json.loads(manifest_path.read_text())
    require(manifest.get("schema") == 1 and set(manifest.get("files", {})) == MEMBERS,
            "archive manifest membership")
    entries = {p.name: p for p in HERE.iterdir()}
    cache = entries.pop("__pycache__", None)
    if cache is not None:
        require(cache.is_dir() and not cache.is_symlink(), "bytecode cache kind")
        require(all(f.is_file() and not f.is_symlink() and f.suffix == ".pyc"
                    for f in cache.iterdir()), "ordinary bytecode cache contents")
    require(set(entries) == MEMBERS | {"archive_manifest.json"}, "archive closure")
    for name, expected in manifest["files"].items():
        path = HERE / name
        require(path.is_file() and not path.is_symlink(), "archive file kind")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
                "archive hash: " + name)


def load_current_source():
    # Ignore bytecode, including timestamp-valid stale caches from ordinary API use.
    global cm
    path = HERE / "count_moran.py"
    module = types.ModuleType("count_moran_verified_source")
    module.__file__ = str(path)
    exec(compile(path.read_bytes(), str(path), "exec"), module.__dict__)
    cm = module


def fine_kernel(parts, r, lam):
    """Independent ordered events on original vertices, not count-rate formulas."""
    block_of = []
    for i, n in enumerate(parts):
        block_of.extend([i] * n)
    N = len(block_of)
    neighbors = [[v for v in range(N) if block_of[v] != block_of[u]]
                 for u in range(N)]
    counts = {}
    kernel = {}
    for state in range(1 << N):
        bits = [(state >> (N - 1 - v)) & 1 for v in range(N)]
        counts[state] = tuple(sum(bits[v] for v in range(N) if block_of[v] == i)
                              for i in range(len(parts)))
        weights = [r if bit else Q(1) for bit in bits]
        total = sum(weights, Q(0))
        neighbor_fitness = [sum((weights[v] for v in neighbors[u]), Q(0))
                            for u in range(N)]
        row = defaultdict(Q)
        for parent in range(N):
            for target in neighbors[parent]:
                prob = (lam * weights[parent] / (total * len(neighbors[parent]))
                        + (1 - lam) * weights[parent] / (N * neighbor_fitness[target]))
                dest = state
                if bits[parent] != bits[target]:
                    dest ^= 1 << (N - 1 - target)
                row[dest] += prob
        require(sum(row.values(), Q(0)) == 1 and all(v >= 0 for v in row.values()),
                "fine native probability row")
        kernel[state] = dict(row)
    return kernel, counts


def boundary_reachability(kernel, boundary):
    reverse = {s: set() for s in kernel}
    for s, row in kernel.items():
        for t, weight in row.items():
            if weight > 0:
                reverse[t].add(s)
    reached = {boundary}
    stack = [boundary]
    while stack:
        t = stack.pop()
        for s in reverse[t] - reached:
            reached.add(s)
            stack.append(s)
    return reached


def check_fixture(parts, r, lam):
    fine, counts = fine_kernel(parts, r, lam)
    coarse = cm.native_count_kernel(parts, r, lam)
    for state, row in fine.items():
        aggregated = defaultdict(Q)
        for target, weight in row.items():
            aggregated[counts[target]] += weight
        expected = coarse[counts[state]]
        require(all(aggregated.get(t, Q(0)) == expected.get(t, Q(0)) for t in coarse),
                "strong lumpability/native count rates")
    zero, full = tuple(0 for _ in parts), tuple(parts)
    all_fine = (1 << sum(parts)) - 1
    require(boundary_reachability(fine, 0) == set(fine) - {all_fine}
            and boundary_reachability(fine, all_fine) == set(fine) - {0},
            "fine absorption reachability")
    require(boundary_reachability(coarse, zero) == set(coarse) - {full}
            and boundary_reachability(coarse, full) == set(coarse) - {zero},
            "count absorption reachability")
    fine_values = cm.harmonic_values(fine, 0, all_fine)
    values = cm.harmonic_values(coarse, zero, full)
    require(all(fine_values[s] == values[counts[s]] for s in fine),
            "fine/count exact fixation agreement")
    if r == 1:
        singles = [1 << v for v in range(sum(parts))]
        require(sum((fine_values[s] for s in singles), Q(0)) / sum(parts)
                == Q(1, sum(parts)), "neutral uniform-singleton mean")
    if cm.universal_count_attractivity_criterion(parts, lam):
        for x in coarse:
            bx, ax = cm.count_rates(parts, x, r, lam)
            for y in coarse:
                if not all(a <= b for a, b in zip(x, y)):
                    continue
                by, ay = cm.count_rates(parts, y, r, lam)
                for i in range(len(parts)):
                    if x[i] == y[i]:
                        require(bx[i] <= by[i] and ax[i] >= ay[i],
                                "finite sufficient-class generator order")
    return values, len(fine), len(coarse)


def exact_obstruction(data):
    expected = {"part_sizes": [1, 1, 2], "fitness": "20", "lambda_Bd": "1/2",
                "lower_counts": [0, 0, 1], "upper_counts": [1, 0, 1], "target_part": 1}
    require(data == expected, "audited obstruction fixture")
    parts, r, lam = tuple(data["part_sizes"]), Q(data["fitness"]), Q(data["lambda_Bd"])
    x, y, i = tuple(data["lower_counts"]), tuple(data["upper_counts"]), data["target_part"]
    fine, count_map = fine_kernel(parts, r, lam)
    coarse = cm.native_count_kernel(parts, r, lam)
    rates = []
    for k in (x, y):
        representative = next(s for s in fine if count_map[s] == k)
        direct = sum((weight for t, weight in fine[representative].items()
                      if count_map[t][i] >= 1), Q(0))
        birth, _ = cm.count_rates(parts, k, r, lam)
        require(direct == birth[i], "native obstruction increasing-function rate")
        rates.append(direct)
    require(rates == [Q(335, 1012), Q(1655, 5166)], "exact birth-rate witness")
    require(rates[1] - rates[0] == Q(-27875, 2613996), "strict birth-rate reverse order")
    fine_values = cm.harmonic_values(fine, 0, (1 << sum(parts)) - 1)
    values = cm.harmonic_values(coarse, (0, 0, 0), parts)
    require(all(fine_values[s] == values[count_map[s]] for s in fine),
            "obstruction fine/count fixation")
    require(values[x] == Q(1015463620779209651900, 1243178519205529924639),
            "lower true fixation reference")
    require(values[y] == Q(3688254807710756250595, 3729535557616589773917),
            "upper true fixation reference")
    require(values[y] > values[x], "witness must not be a fixation counterexample")
    require(not cm.universal_count_attractivity_criterion(parts, lam),
            "universal criterion on excluded graph")
    return {"lower_birth_rate": str(rates[0]), "upper_birth_rate": str(rates[1]),
            "upper_minus_lower_birth_rate": str(rates[1] - rates[0]),
            "fixation_lower": str(values[x]), "fixation_upper": str(values[y]),
            "upper_minus_lower_fixation": str(values[y] - values[x]),
            "is_fixation_counterexample": False}


def cli_check(expected):
    cmd = [sys.executable]
    if sys.flags.optimize:
        cmd.append("-O")
    cmd += [str(HERE / "count_moran.py"), "--parts", "2,3", "--fitness", "2",
            "--bd", "1/2", "--initial", "1,0"]
    result = subprocess.run(cmd, text=True, capture_output=True, check=True)
    data = json.loads(result.stdout)
    require(data["fixation_probability_exact"] == str(expected),
            "standalone CLI exact probability")
    require(data["count_states"] == 12 and data["binary_vertex_states"] == 32,
            "standalone state reduction")
    require(data["theorem_mechanically_certified"] is False,
            "CLI theorem certification scope")
    return data


def compute_checks():
    load_current_source()
    examples = json.loads((HERE / "examples.json").read_text())
    require(examples.get("schema") == 1, "input schema")
    require(examples["part_sizes"] == [[2, 3], [2, 2, 2], [1, 1, 2]]
            and examples["fitness"] == ["1/2", "1", "2"]
            and examples["lambda_Bd"] == ["0", "1/2", "1"],
            "audited finite fixture coverage")
    fixture_count = fine_rows = coarse_rows = order_comparisons = 0
    cached = {}
    rs = list(map(Q, examples["fitness"]))
    for sizes in examples["part_sizes"]:
        parts = tuple(sizes)
        for ls in examples["lambda_Bd"]:
            lam = Q(ls)
            curves = {}
            for r in rs:
                values, nf, nc = check_fixture(parts, r, lam)
                curves[r] = values
                cached[(parts, r, lam)] = values
                fixture_count += 1
                fine_rows += nf
                coarse_rows += nc
            # Finite corroboration, including the unequal graph:
            # it is not an infinite extension of the analytic sufficient class.
            for r1 in rs:
                for r2 in rs:
                    if r1 > r2:
                        continue
                    for x, vx in curves[r1].items():
                        for y, vy in curves[r2].items():
                            if all(a <= b for a, b in zip(x, y)):
                                require(vx <= vy, "finite joint fixation order")
                                order_comparisons += 1
    obstruction = exact_obstruction(examples["obstruction"])
    cli = cli_check(cached[((2, 3), Q(2), Q(1, 2))][(1, 0)])
    return {"fixture_count": fixture_count, "fine_state_rows": fine_rows,
            "count_state_rows": coarse_rows, "finite_joint_order_comparisons": order_comparisons,
            "obstruction": obstruction, "CLI_example": cli,
            "theorem_mechanically_certified": False,
            "general_graph_fixation_question_resolved": False}


def main():
    archive_check()
    expected = json.loads((HERE / "certificate.json").read_text())
    require(expected.get("schema") == 1, "certificate schema")
    require(expected.get("theorem_mechanically_certified") is False
            and expected.get("general_graph_fixation_question_resolved") is False,
            "certificate theorem scope")
    actual = compute_checks()
    require(actual == {k: v for k, v in expected.items() if k != "schema"},
            "finite exact certificate")
    print(json.dumps({"status": "PASS", "fixtures": actual["fixture_count"],
                      "fine_state_rows": actual["fine_state_rows"],
                      "count_state_rows": actual["count_state_rows"],
                      "finite_joint_order_comparisons": actual["finite_joint_order_comparisons"],
                      "obstruction_is_fixation_counterexample": False,
                      "theorem_mechanically_certified": False,
                      "general_graph_fixation_question_resolved": False}, sort_keys=True))


if __name__ == "__main__":
    main()
