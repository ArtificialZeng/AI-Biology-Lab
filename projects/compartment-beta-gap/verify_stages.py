#!/usr/bin/env python3
"""Read-only stage checks. Infinite-state probability claims need stage-proof.md."""
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
ARCHIVE = "evidence/stage-certificate.json"


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


def positive_shift(expr, variable, lower):
    u = sp.symbols('u', nonnegative=True)
    poly = sp.Poly(sp.expand(expr.subs(variable,u+lower)),u)
    coefficients = poly.all_coeffs()
    require(all(c>=0 for c in coefficients),f'Negative shifted coefficient: {expr}')
    require(poly.eval(0)>0,f'Nonpositive boundary: {expr}')
    return [str(c) for c in coefficients]


def compute_certificate():
    require(sp.__version__=='1.13.3','Pin SymPy 1.13.3')
    x,r,beta=sp.symbols('x r beta',positive=True)
    certificates={
      'carrier_birth_x_ge_4':positive_shift(x**4/8-x**2/8-1,x,4),
      'carrier_death_x_ge_4':positive_shift(x**4/64-x**2/8,x,4),
      'fill_before_half_window_x_ge_16':positive_shift(x**3-8*x**2-4,x,16),
      'four_power_induction_r_ge_4':positive_shift(4*r**4-(r+1)**4,r,4),
      'C_tail_x_ge_16':positive_shift(x**7-64,x,16),
      'B_tail_x_ge_16':positive_shift(x**3-2048,x,16),
      'summed_stage_error_x_ge_16':positive_shift(64*x**7-42*x**6-10*x**4-2,x,16),
    }
    require(sp.simplify(8*(x**4/4)/x**8-2/x**4)==0,'C filling time')
    require(sp.simplify(8*(x**2/16)/(x**4/8)-4/x**2)==0,'B filling time')
    require(sp.simplify((x**2/8)*(1/(2*x))-x/16)==0,'Import Poisson mean')
    require(sp.simplify(4/(x/16)-64/x)==0,'Import Chebyshev')
    require(sp.Rational(1,32)>sp.Rational(1,1024),'Imports exceed enzyme target')
    require(sp.Rational(1,512)<sp.Rational(1,2),'Yule exponent comparison')
    require(sp.Rational(2,254**2)<1,'Early Yule Chebyshev')
    require(sp.simplify((x/1024)*x**8*(2**20/x)-1024*x**8)==0,'Late Poisson mean')
    require(sp.Rational(1024,769**2)<1,'Late Poisson Chebyshev')
    require(4**4==4**4,'Power induction anchor')
    require(sp.Integer(4)*4**4>5**4,'Power induction ratio')
    require(sp.simplify(8*(1-beta)-1-(7-8*beta))==0,'Occupation exponent')
    require(1-sp.Rational(256,1024)==sp.Rational(3,4),'Stage success probability')
    require(2**21//1024==2048,'Explosion time sum')
    cases=0
    for k in range(5):
        for cold in combinations_with_replacement(range(4),k):
            C,B=k+1,sum(v>0 for v in cold)
            rates=Counter(B_birth=sum(v==0 for v in cold))
            for i in range(C):
                for j in range(i+1,C):
                    rates['C_death']+=2
                    if i==0 and cold[j-1]>0:
                        rates['B_death']+=2
                        rates['import']+=2
                    elif i>0 and cold[i-1]>0 and cold[j-1]>0:
                        rates['B_death']+=2
            require(rates['B_birth']==C-1-B,'Finite carrier birth identity')
            require(rates['B_death']==B*(B+1),'Finite carrier death identity')
            require(rates['C_death']==C*(C-1),'Finite coagulation identity')
            require(rates['import']==2*B,'Finite import event rate')
            cases+=1
    require(cases==70,'Finite coverage')
    base=Path(__file__).resolve().parent
    return {'status':'pass','arithmetic':'exact integer/rational/symbolic','shifted_positive_coefficients':certificates,
      'finite_event_configurations':cases,'finite_coverage':{'outside_compartments':[0,4],'outside_enzymes_per_compartment':[0,3]},
      'stage':{'x0':1024,'initial_tag_enzymes':1,'initial_tag_substrate':'2^80','substrate_exponent':8,'enzyme_level':'x/1024','next_enzyme_level':'x/512','length':'2^20/x','failure_bound':'128/x','all_stages_success_probability_at_least':'3/4','all_stages_time_at_most':2048},
      'occupation':'delta<1/8; sum of x_n^(8delta-1) converges','beta_range':'7/8<beta<1',
      'proof_sha256':hashlib.sha256((base/'stage-proof.md').read_bytes()).hexdigest(),
      'limitations':'These finite checks do not prove the guarded couplings, stage iteration, finite-compensator implication, conditional mark product, reachability or infinite-state explosion theorem; those require stage-proof.md and mathematical review.'}



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
