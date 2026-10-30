#!/usr/bin/env python3
"""Read-only exact finite checks. Infinite-state probability claims need proof.md."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path
import re
import sys

import sympy as sp

PROJECT = Path(__file__).resolve().parent
MANIFEST = "archive_manifest.json"
ARCHIVE = "evidence/certificate.json"


class VerificationError(Exception):
    """A failed verification or an unsafe regeneration destination."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        require(key not in result, f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def load_json(path: Path) -> object:
    try:
        return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique_object)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise VerificationError(f"Cannot read valid JSON from {path.name}: {exc}") from exc


def verify_manifest() -> None:
    manifest = load_json(PROJECT / MANIFEST)
    require(isinstance(manifest, dict), "Manifest must be a JSON object")
    require(manifest.get("schema") == 1, "Unsupported manifest schema")
    hashes = manifest.get("sha256")
    require(isinstance(hashes, dict) and bool(hashes), "Missing manifest file map")
    actual = {
        path.relative_to(PROJECT).as_posix()
        for path in PROJECT.rglob("*") if path.is_file() and path != PROJECT / MANIFEST
    }
    require(set(hashes) == actual,
            f"Manifest inventory mismatch: missing={sorted(actual-set(hashes))}, "
            f"stale={sorted(set(hashes)-actual)}")
    for name, expected in sorted(hashes.items()):
        require(isinstance(name, str) and Path(name).as_posix() == name
                and not Path(name).is_absolute() and ".." not in Path(name).parts,
                f"Invalid manifest path: {name}")
        require(isinstance(expected, str) and re.fullmatch(r"[0-9a-f]{64}", expected) is not None,
                f"Invalid SHA-256 for {name}")
        require(digest(PROJECT / name) == expected, f"Integrity mismatch: {name}")


def compute_certificate() -> dict[str, object]:
    require(sp.__version__ == "1.13.3", f"Expected SymPy 1.13.3; found {sp.__version__}")
    C,B,e,s,le3,le2,ls = sp.symbols("C B e s le3 le2 ls", positive=True)
    A = s/(C+2)**2
    D = C/(B+2)**2
    H = 2*B/((e+3)*le3*le2)
    J = e/(4*ls**2)
    expected = (s/sp.Integer(16))*(C/(C+2))**2*(B/(B+2))**4*(e/(e+3))**4/(le3**4*le2**4*ls**8)
    require(sp.cancel(A*D**2*H**4*J**4-expected) == 0, "Product identity failed")
    weights = sp.Matrix([1,2,4,4])
    powers = sp.Matrix([[-2,1,0,0],[0,-2,1,0],[0,0,-1,1],[1,0,0,0]])
    require(powers*weights == sp.Matrix([0,0,0,1]), "Exponent balance failed")
    require(sum(weights) == 11, "Total exponent failed")
    require(sp.Rational(1,16)*sp.Rational(1,2)**10 == sp.Rational(1,2**14), "Product lower-bound constant failed")
    require(sp.Rational(1,2**14)*2**36 == 4**11, "Power comparison constant failed")
    r = sp.symbols("r", positive=True)
    L = 4096/r**2
    require(sp.simplify(r*L-16*(12+2/r)-38-(4064/r-230)) == 0, "Log-threshold rearrangement failed")
    require(4064-230 > 0 and 4096-16 > 0, "Threshold positivity failed")
    require(4096//2-2*12-6 == 2018, "Second-case constant failed")
    require(sp.Rational(1,4)-sp.Rational(3,32) > 0, "Initial-function bound failed")
    require(1-sp.Rational(3,32)-sp.Rational(3,32) == sp.Rational(13,16), "Probability arithmetic failed")
    beta = sp.Rational(11,12)
    r0 = 11*beta-10
    L0 = 4096/r0**2
    a0 = 1/(L0+1)
    M0 = sp.ceiling(64/a0)
    require((r0,L0,a0,M0) == (sp.Rational(1,12),589824,sp.Rational(1,589825),37748800), "Explicit parameter witness failed")

    cases = 0
    for k in range(5):
        for cold in combinations_with_replacement(range(4), k):
            c,b = k+1,sum(n>0 for n in cold)
            rates = Counter()
            rates["B_birth"] = sum(n==0 for n in cold)
            for i in range(c):
                for j in range(i+1,c):
                    rates["C_death"] += 2
                    if i == 0 and cold[j-1] > 0:
                        rates["B_death"] += 2
                        rates["positive_import_event"] += 2
                    elif i > 0 and cold[i-1] > 0 and cold[j-1] > 0:
                        rates["B_death"] += 2
            require(rates["B_birth"] == c-1-b, f"Carrier birth rate failed for {cold}")
            require(rates["C_death"] == c*(c-1), f"Coagulation rate failed for {cold}")
            require(rates["B_death"] == b*(b+1), f"Carrier death rate failed for {cold}")
            require(rates["positive_import_event"] == 2*b, f"Positive import rate failed for {cold}")
            cases += 1
    require(cases == 70, "Finite configuration coverage changed")
    return {
        "status": "pass",
        "arithmetic": "exact rational/integer/symbolic",
        "identity": "A D^2 H^4 J^4; weights (1,2,4,4), total exponent 11",
        "finite_event_configurations": cases,
        "coverage": {"outside_compartments":[0,4], "outside_enzymes_per_compartment":[0,3]},
        "example": {"beta":"11/12", "r":"1/12", "L":589824, "a":"1/589825", "M":37748800, "N":"2^37748800"},
        "limitations": "Finite algebra and event checks only. Logarithmic inequalities, infinite-state localization, compensator, marked CTMC coupling and reachability depend on proof.md; this script does not prove explosion.",
        "proof_sha256": digest(PROJECT / "proof.md"),
        "sympy": sp.__version__,
    }


def encoded(value: object) -> str:
    return json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"


def regenerate(output: str, result: dict[str, object]) -> None:
    requested = Path(output).expanduser()
    require(requested.is_absolute(), "Regeneration needs an explicit absolute output path")
    target = requested.resolve()
    require(not target.is_relative_to(PROJECT), "Regeneration output must be outside the project")
    require(target.parent.is_dir(), "Regeneration output parent directory must already exist")
    try:
        with target.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(encoded(result))
    except OSError as exc:
        raise VerificationError(f"Cannot exclusively create regeneration output: {exc}") from exc


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command")
    sub.add_parser("verify", help="Read-only archive and exact finite checks (default)")
    regen = sub.add_parser("regenerate", help="Exclusively write a new certificate outside the project")
    regen.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        if args.command == "regenerate":
            result = compute_certificate()
            regenerate(args.output, result)
        else:
            verify_manifest()
            result = compute_certificate()
            require(load_json(PROJECT / ARCHIVE) == result, "Recomputed certificate differs from archived JSON")
        sys.stdout.write(encoded(result))
        return 0
    except (VerificationError, OSError) as exc:
        print(f"verification failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
