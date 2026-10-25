"""Offline, source-pinned real-MEMOTE integration test for a direction-only fix.

Run: .venv-integration/bin/python verify_direction_patch.py
Only the negative-only branch is added. Installed packages are never modified.
The candidate is imported from a temporary copy of the real MEMOTE package;
the complete installed ATP detection implementation and real COBRApy/GLPK run.
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import importlib.metadata
import importlib.util
import inspect
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile

BASE = Path(__file__).resolve().parent
PIN = "bf76477b4087feefc8ad429ebd492e30ef10eebc"
HELPERS_SHA256 = "bbb3112dd94f8e6b92853a38483e9c623396d03dd8091e240a24420822963311"
OLD = "        if rxn.reversibility:\n            rxn.bounds = -1, 1\n        else:\n            rxn.bounds = 0, 1\n"
NEW = "        if rxn.reversibility:\n            rxn.bounds = -1, 1\n        elif rxn.upper_bound <= 0 and rxn.lower_bound < 0:\n            rxn.bounds = -1, 0\n        else:\n            rxn.bounds = 0, 1\n"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def child_run():
    import cobra
    import swiglpk
    from memote.support import consistency, helpers

    specs = [
        ("MNXM3", "ATP", "C10H12N5O13P3", -4, -1),
        ("MNXM2", "water", "H2O", 0, -1),
        ("MNXM7", "ADP", "C10H12N5O10P2", -3, 1),
        ("MNXM1", "proton", "H", 1, 1),
        ("MNXM9", "phosphate", "HO4P", -2, 1),
    ]

    def build(sigma, bounds):
        model = cobra.Model("direction_invariance_toy")
        model.solver = "glpk"
        h = {}
        for mnx, name, formula, charge, coefficient in specs:
            met = cobra.Metabolite(mnx + "_c", name=name, formula=formula,
                                   charge=charge, compartment="c")
            h[met] = coefficient
        reaction = cobra.Reaction("R", lower_bound=bounds[0], upper_bound=bounds[1])
        reaction.add_metabolites({met: sigma * value for met, value in h.items()})
        model.add_reactions([reaction])
        model.compartments = {"c": "cytosol"}
        assert reaction.check_mass_balance() == {}
        assert not model.boundary
        return model, h

    def summarize(model, support=None):
        sol = model.optimize(raise_error=True)
        residuals = [sum(r.metabolites.get(m, 0) * float(sol.fluxes[r.id])
                         for r in model.reactions) for m in model.metabolites]
        assert all(r.check_mass_balance() == {} for r in model.reactions)
        assert not model.boundary
        assert max(map(abs, residuals)) < 1e-9
        return {"objective": float(sol.objective_value), "status": sol.status,
                "R_bounds": list(model.reactions.R.bounds),
                "diagnostic_bounds": list(model.reactions[-1].bounds),
                "fluxes": {k: float(v) for k, v in sol.fluxes.items()},
                "max_stoichiometric_residual": max(map(abs, residuals)),
                "all_reactions_element_charge_balanced": True,
                "boundary_count": 0, "support": support}

    def manual(sigma, bounds, normalized):
        model, h = build(sigma, bounds)
        if normalized:
            helpers.close_boundaries_sensibly(model)
        diagnostic = cobra.Reaction("D", lower_bound=0, upper_bound=1)
        diagnostic.add_metabolites(h)
        model.add_reactions([diagnostic])
        model.objective = diagnostic
        return summarize(model)

    cases = [
        ("hydrolysis_forward", 1, (0, 1000), 0),
        ("hydrolysis_reverse", -1, (-1000, 0), 0),
        ("charging_forward", -1, (0, 1000), 1),
        ("charging_reverse", 1, (-1000, 0), 1),
        ("reversible_forward", 1, (-1000, 1000), 1),
        ("reversible_reverse", -1, (-1000, 1000), 1),
    ]
    rows = []
    for name, sigma, bounds, expected_native in cases:
        native = manual(sigma, bounds, False)
        assert abs(native["objective"] - expected_native) < 1e-9
        normalized = manual(sigma, bounds, True)
        model, _ = build(sigma, bounds)
        support = consistency.detect_energy_generating_cycles(model, "MNXM3")
        entry = summarize(model, support)
        assert abs(entry["objective"] - normalized["objective"]) < 1e-9
        assert support == (["R"] if normalized["objective"] > 0 else [])
        rows.append({"case": name, "sigma": sigma, "original_bounds": list(bounds),
                     "native_manual": native, "normalized_manual": normalized,
                     "complete_ATP_entry": entry})

    # The helper deliberately removes fixed/nonzero minimum magnitudes. For
    # these intervals compare equivalent orientations only AFTER this policy;
    # the native closed-network LP can be infeasible and is not the baseline.
    relaxed_rows = []
    for interval_name, forward_bounds, reverse_bounds in [
        ("strict_one_way", (2, 1000), (-1000, -2)),
        ("fixed_nonzero", (5, 5), (-5, -5)),
    ]:
        for direction, sigma_forward in [("hydrolysis", 1), ("charging", -1)]:
            for orientation, sigma, bounds in [
                ("forward", sigma_forward, forward_bounds),
                ("reverse", -sigma_forward, reverse_bounds),
            ]:
                model, _ = build(sigma, bounds)
                support = consistency.detect_energy_generating_cycles(model, "MNXM3")
                relaxed_rows.append({"case": interval_name + "_" + direction + "_" + orientation,
                                     "original_bounds": list(bounds), "sigma": sigma,
                                     "complete_ATP_entry": summarize(model, support)})

    branch_rows = []
    bounds_cases = [(-1000, 0), (-1000, -2), (-5, -5), (-1000, 1000),
                    (0, 1000), (5, 1000), (5, 5), (0, 0)]
    for bounds in bounds_cases:
        model, _ = build(1, bounds)
        helpers.close_boundaries_sensibly(model)
        branch_rows.append({"input": list(bounds), "output": list(model.reactions.R.bounds)})
    for bounds in [(-1000, 0), (0, 1000), (-1000, 1000), (0, 0)]:
        model, h = build(1, (0, 1000))
        boundary = cobra.Reaction("EX_water", lower_bound=bounds[0], upper_bound=bounds[1])
        water = next(m for m in h if m.id == "MNXM2_c")
        boundary.add_metabolites({water: -1})
        model.add_reactions([boundary])
        assert boundary in model.boundary
        helpers.close_boundaries_sensibly(model)
        assert boundary.bounds == (0, 0)
        branch_rows.append({"boundary_input": list(bounds), "output": list(boundary.bounds)})

    return {"python": platform.python_version(), "versions": {
                name: importlib.metadata.version(name)
                for name in ["memote", "cobra", "optlang", "swiglpk"]},
            "GLPK_version": swiglpk.glp_version(),
            "helpers_sha256": sha(Path(inspect.getsourcefile(helpers)).read_bytes()),
            "consistency_sha256": sha(Path(inspect.getsourcefile(consistency)).read_bytes()),
            "ATP_function_sha256": sha(inspect.getsource(consistency.detect_energy_generating_cycles).encode()),
            "cases": rows, "relaxed_equivalent_cases": relaxed_rows,
            "normalization_cases": branch_rows}


def run(output):
    pinned = BASE / "fixtures/src/memote/support/helpers.py"
    original = pinned.read_bytes()
    assert sha(original) == HELPERS_SHA256
    text = original.decode()
    assert text.count(OLD) == 1
    candidate = text.replace(OLD, NEW)
    patch = "".join(difflib.unified_diff(text.splitlines(keepends=True),
                    candidate.splitlines(keepends=True),
                    fromfile="a/src/memote/support/helpers.py",
                    tofile="b/src/memote/support/helpers.py"))
    if (BASE / 'patches/preserve-negative-direction.patch').read_text() != patch:
        raise RuntimeError('Published candidate patch mismatch')
    installed_root = Path(importlib.util.find_spec("memote").origin).parent
    assert importlib.metadata.version("memote") == "0.17.0"
    assert (installed_root / "support/helpers.py").read_bytes() == original
    assert (installed_root / "support/consistency.py").read_bytes() == (BASE / "fixtures/src/memote/support/consistency.py").read_bytes()

    # Both variants run in fresh processes. Only the candidate's real package
    # copy changes, by the two added source lines; no monkeypatch is used.
    variants = {}
    for label in ["original", "candidate"]:
        with tempfile.TemporaryDirectory(prefix="direction-check-", dir=output.parent) as temp:
            root = Path(temp)
            shutil.copytree(installed_root, root / "memote", ignore=shutil.ignore_patterns("__pycache__"))
            helper_path = root / "memote/support/helpers.py"
            if label == "candidate":
                helper_path.write_text(candidate)
            environment = os.environ.copy()
            environment["PYTHONPATH"] = str(root)
            child = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--child"],
                                   cwd=BASE, env=environment, capture_output=True, text=True,
                                   timeout=90, check=True)
            variants[label] = json.loads(child.stdout)

    assert variants["original"]["helpers_sha256"] == HELPERS_SHA256
    assert variants["candidate"]["helpers_sha256"] == sha(candidate.encode())

    old = {r["case"]: r for r in variants["original"]["cases"]}
    new = {r["case"]: r for r in variants["candidate"]["cases"]}
    assert old["hydrolysis_reverse"]["native_manual"]["objective"] == 0
    assert old["hydrolysis_reverse"]["complete_ATP_entry"]["objective"] == 1
    assert old["charging_reverse"]["native_manual"]["objective"] == 1
    assert old["charging_reverse"]["complete_ATP_entry"]["objective"] == 0
    for row in new.values():
        assert row["complete_ATP_entry"]["objective"] == row["native_manual"]["objective"]
    pairs = [("hydrolysis_forward", "hydrolysis_reverse"),
             ("charging_forward", "charging_reverse"),
             ("reversible_forward", "reversible_reverse")]
    for forward, reverse in pairs:
        assert new[forward]["complete_ATP_entry"]["objective"] == new[reverse]["complete_ATP_entry"]["objective"]
    relaxed = {r["case"]: r for r in variants["candidate"]["relaxed_equivalent_cases"]}
    for interval_name in ["strict_one_way", "fixed_nonzero"]:
        for direction, expected in [("hydrolysis", 0), ("charging", 1)]:
            forward = interval_name + "_" + direction + "_forward"
            reverse = interval_name + "_" + direction + "_reverse"
            assert relaxed[forward]["complete_ATP_entry"]["objective"] == expected
            assert relaxed[reverse]["complete_ATP_entry"]["objective"] == expected
    for before, after in zip(variants["original"]["normalization_cases"],
                             variants["candidate"]["normalization_cases"]):
        if "input" in before and before["input"][0] < 0 and before["input"][1] <= 0:
            assert before["output"] == [0, 1] and after["output"] == [-1, 0]
        else:
            assert before == after
    assert variants["original"]["consistency_sha256"] == variants["candidate"]["consistency_sha256"]
    assert variants["original"]["ATP_function_sha256"] == variants["candidate"]["ATP_function_sha256"]
    return {"status": "passed", "source_commit": PIN,
            "source_URL": f"https://github.com/opencobra/memote/blob/{PIN}/src/memote/support/helpers.py#L676-L703",
            "original_helpers_sha256": HELPERS_SHA256,
            "candidate_helpers_sha256": sha(candidate.encode()), "patch_sha256": sha(patch.encode()),
            "variants": variants, "equivalent_orientation_pairs": pairs,
            "assertions": ["Original reverse-only false positive 0 -> 1 reproduced",
                           "Original reverse-only false negative 1 -> 0 reproduced",
                           "Candidate restores all six original-model ATP-LP optima under diagnostic cap one",
                           "All three equivalent-orientation ATP outcome pairs agree after fix",
                           "Four strict-one-way/fixed-nonzero equivalent-orientation pairs agree after documented magnitude relaxation",
                           "Positive, reversible, fixed-zero and boundary normalization unchanged",
                           "Fixed-negative magnitudes relaxed while retaining negative direction"],
            "scope": ["Real COBRApy models, real complete MEMOTE ATP function, GLPK",
                      "Installed packages untouched; temporary full-package copy imports candidate source",
                      "Fixed-zero reopening remains the documented intentional normalization policy",
                      "No full MEMOTE test suite, CLI report, genome-scale claim or maintainer confirmation",
                      "LP floating-point checks supplement earlier exact integer toy certificates",
                      "Original native diagnostic bounds [0,1]; MEMOTE diagnostic [0,1000], normalized R caps objective at one"]}


if __name__ == "__main__":
    if "--child" in sys.argv:
        print(json.dumps(child_run()))
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--output', type=Path, default=BASE.parents[1]/'artifacts/memote-direction-auditor/verification.json')
        args = parser.parse_args()
        if not __debug__:
            parser.error('Assertions must be enabled')
        args.output.parent.mkdir(parents=True, exist_ok=True)
        result = run(args.output)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
        print(json.dumps({"status": result["status"], "assertions": result["assertions"],
                          "patch_sha256": result["patch_sha256"]}, indent=2))
