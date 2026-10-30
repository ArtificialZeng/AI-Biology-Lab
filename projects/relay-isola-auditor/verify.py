#!/usr/bin/env python3
"""Read-only exact verification; regeneration requires exclusive external output."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parent
CERTIFICATE = "evidence/exact-certificate.json"
MANIFEST = "archive_manifest.json"


class VerificationError(Exception):
    """A failed exact identity, archive check or regeneration precondition."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def equal(left, right, message: str) -> None:
    require(sp.expand(left - right) == 0, message)


def encode(value) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n").encode()


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def public_paths() -> list[str]:
    paths = []
    for path in ROOT.rglob("*"):
        rel = path.relative_to(ROOT)
        if "__pycache__" in rel.parts or path.suffix in {".pyc", ".pyo"}:
            continue
        require(not path.is_symlink(), "Project symlinks are not archive inputs: " + str(rel))
        if path.is_file() and rel.as_posix() != MANIFEST:
            paths.append(rel.as_posix())
    return sorted(paths)


def build_certificate() -> dict:
    require(sp.__version__ == "1.13.3", "Install the pinned SymPy 1.13.3 dependency")
    x, y, lam, loss = sp.symbols("x y lambda loss", real=True)
    reactions = [
        ((1, 0), (0, 0), sp.Integer(7)),
        ((2, 0), (3, 0), sp.Integer(4)),
        ((3, 0), (2, 0), sp.Integer(1)),
        ((1, 1), (2, 1), sp.Integer(4)),
        ((1, 2), (0, 0), sp.Integer(1)),
        ((1, 2), (1, 4), sp.Integer(1)),
        ((0, 0), (0, 1), lam),
        ((0, 1), (0, 0), sp.Integer(1)),
    ]
    field = sp.zeros(2, 1)
    columns = []
    for source, product, constant in reactions:
        require(constant == lam or constant.is_positive is True, "Nonpositive fixed rate")
        column = sp.Matrix([product[i] - source[i] for i in range(2)])
        columns.append(column)
        propensity = constant * x ** source[0] * y ** source[1]
        field += column * propensity
        consumed = {i for i in range(2) if column[i] < 0}
        support = {i for i in range(2) if source[i] > 0}
        require(not consumed or support <= consumed, "Strict net-consumption support condition failed")
    require(sp.Matrix.hstack(*columns).rank() == 2, "Stoichiometric rank is not full")
    equal(field[0], x * (-7 + 4*x - x*x + 4*y - y*y), "X field reconstruction failed")
    equal(field[1], lam - y, "Y field reconstruction failed")
    equal(field[0].subs(x, 0), 0, "X=0 is not invariant")
    equal(field[1].subs(y, 0), lam, "Y=0 inward condition failed")

    siphons = []
    for mask in range(1, 4):
        species = {i for i in range(2) if mask & (1 << i)}
        if all(not any(product[i] > 0 for i in species)
               or any(source[i] > 0 for i in species)
               for source, product, _ in reactions):
            siphons.append(sorted(species))
    require(siphons == [[0]], "Exhaustive siphon list differs from {X}")

    equilibrium = sp.expand((field[0] / x).subs(y, lam))
    equal(equilibrium, 1-(x-2)**2-(lam-2)**2, "Circle elimination identity failed")
    # An exact rational parametrization covers the circle except one point.
    t = sp.symbols("t", real=True)
    circle_x = 2 + (1-t*t)/(1+t*t)
    circle_lam = 2 + 2*t/(1+t*t)
    require(sp.cancel(equilibrium.subs({x:circle_x, lam:circle_lam})) == 0,
            "Exact circle parametrization failed")
    equal(equilibrium.subs({x:1, lam:2}), 0, "Missing circle point check failed")

    jacobian = field.jacobian([x, y])
    boundary = jacobian.subs({x:0, y:lam})
    equal(boundary[0, 0], -3-(lam-2)**2, "Boundary spectral certificate failed")
    require(boundary[0, 1] == 0 and boundary[1, 0] == 0 and boundary[1, 1] == -1,
            "Boundary Jacobian is not the claimed diagonal matrix")
    equal(7+lam*lam-4*lam, 3+(lam-2)**2, "R<1 positive denominator certificate failed")
    equal(jacobian[0, 0].subs(y, lam)-equilibrium, 2*x*(2-x),
          "Positive-equilibrium eigenvalue identity failed")

    folds = []
    for parameter in [sp.Integer(1), sp.Integer(3)]:
        point = {x:sp.Integer(2), y:parameter, lam:parameter}
        matrix = jacobian.subs(point)
        right = sp.Matrix([1, 0])
        left = sp.Matrix([1, matrix[0, 1]])
        require(matrix*right == sp.zeros(2, 1), "Fold right kernel failed")
        require(left.T*matrix == sp.zeros(1, 2), "Fold left kernel failed")
        require((left.T*right)[0] == 1, "Fold kernel normalization failed")
        parameter_coefficient = (left.T*field.diff(lam).subs(point))[0]
        quadratic_coefficient = (left.T*sp.Matrix([sp.diff(f, x, 2).subs(point) for f in field]))[0]
        require(parameter_coefficient != 0 and quadratic_coefficient != 0,
                "Fold is not an ordinary saddle-node")
        equal(quadratic_coefficient, -4, "Fold quadratic coefficient differs")
        folds.append({"lambda":int(parameter), "state":[2,int(parameter)],
                      "jacobian":[[int(v) for v in matrix.row(i)] for i in range(2)],
                      "left_kernel":[int(v) for v in left],
                      "parameter_coefficient":int(parameter_coefficient),
                      "quadratic_coefficient":int(quadratic_coefficient)})

    equal(field[0].subs(x, 4), -12-4*(y-2)**2, "Compact rectangle X-edge identity failed")
    equal(field[1].subs(y, 4), lam-4, "Compact rectangle Y-edge identity failed")

    # Deliberately demonstrate the limit of the one-parameter assertion.
    enlarged_field = field + sp.Matrix([(7-loss)*x, 0])
    path_loss = 4+4*x-x*x
    path = {y:sp.Integer(2), lam:sp.Integer(2), loss:path_loss}
    for component in enlarged_field.subs(path):
        equal(component, 0, "Full-rate-space equilibrium path failed")
    equal(path_loss, 4+x*(4-x), "Path positivity identity failed")
    equal(path_loss.subs(x, 3), 7, "Path does not meet the original circle")
    equal(path_loss.subs(x, 0), 4, "Path boundary limit rate failed")
    boundary_extended = enlarged_field.jacobian([x,y]).subs({x:0,y:2,lam:2,loss:4})
    require(boundary_extended == sp.diag(0,-1), "Enlarged-space boundary crossing identity failed")

    return {
        "schema":"relay-isola-exact-v1", "arithmetic":"exact symbolic rational arithmetic",
        "sympy_version":sp.__version__,
        "parameter_domain":"lambda>0; only the emptyset->Y rate varies",
        "reactions":[{"source":list(u),"product":list(v),"constant":str(k)} for u,v,k in reactions],
        "stoichiometric_rank":2, "siphons_zero_based":siphons,
        "vector_field":[str(sp.expand(f)) for f in field],
        "positive_equilibrium_identity":str(equilibrium),
        "positive_component":"y=lambda; (x-2)^2+(lambda-2)^2=1; complete closed connected circle",
        "component_coordinate_bounds":[1,3],
        "boundary_equilibrium":"(0,lambda)",
        "boundary_eigenvalues":["-3-(lambda-2)^2","-1"],
        "invasion_number":"4*lambda/(7+lambda^2)<1",
        "positive_equilibrium_eigenvalues":["2*x*(2-x)","-1"],
        "count":"two for 1<lambda<3; one double at lambda=1,3; zero otherwise",
        "folds":folds,
        "compact_box":"[0,4]^2 invariant for lambda in [1/2,7/2]",
        "full_rate_space_path":{"domain":"0<x<=3; y=lambda=2",
                                "X_loss_rate":"4+4*x-x^2=4+x*(4-x)>0",
                                "start":"x=3, X_loss=7 lies on original isola",
                                "boundary_limit":"x=0, X_loss=4; all kinetic rates still positive",
                                "consequence":"This witness does not refute a full independent kinetic-rate-space version."},
        "limitations":["Isolas are classical; this is a source-specific given-family counterexample.",
                       "The local relay theorem and Feng-Gavish model are unaffected.",
                       "Fresh substantive AI review accepted the stated mathematics; no formal release-gate certification, external human peer review or priority guarantee."]
    }


def proposed_manifest(certificate: bytes) -> dict:
    paths = public_paths()
    require(CERTIFICATE in paths, "Archived certificate is missing")
    return {"schema":"relay-isola-archive-v1", "algorithm":"sha256",
            "files":{name:digest(certificate if name == CERTIFICATE else (ROOT/name).read_bytes())
                     for name in paths}}


def verify() -> dict:
    manifest = json.loads((ROOT/MANIFEST).read_bytes())
    require(manifest.get("schema") == "relay-isola-archive-v1", "Unknown manifest schema")
    require(manifest.get("algorithm") == "sha256", "Unknown manifest digest algorithm")
    require(sorted(manifest.get("files",{})) == public_paths(), "Manifest does not cover every project input")
    for name, expected in manifest["files"].items():
        require(digest((ROOT/name).read_bytes()) == expected, "Input hash mismatch: " + name)
    recomputed = encode(build_certificate())
    archived = (ROOT/CERTIFICATE).read_bytes()
    require(recomputed == archived, "Recomputed exact certificate differs from the archive")
    return {"ok":True, "mode":"verify", "archive_modified":False,
            "public_files_checked":len(manifest["files"]),
            "certificate_sha256":digest(archived),
            "checks":["eight-reaction reconstruction","strict chemicality support","full stoichiometric rank",
                      "exhaustive siphons","complete equilibrium circle","strict boundary stability",
                      "ordinary folds","compact rectangle","full-rate-space scope path"]}


def regenerate(directory: Path) -> dict:
    output = directory.expanduser().resolve()
    require(output != ROOT and ROOT not in output.parents, "Regeneration output must be outside the project")
    require(not output.exists(), "Regeneration output must not already exist")
    certificate = encode(build_certificate())
    manifest = encode(proposed_manifest(certificate))
    output.mkdir(parents=False, exist_ok=False)
    (output/"evidence").mkdir()
    (output/CERTIFICATE).write_bytes(certificate)
    (output/MANIFEST).write_bytes(manifest)
    return {"ok":True,"mode":"regenerate","archive_modified":False,
            "exclusive_output":str(output),"certificate_sha256":digest(certificate)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("verify")
    regen = sub.add_parser("regenerate")
    regen.add_argument("--output-dir",type=Path,required=True)
    args = parser.parse_args()
    try:
        result = verify() if args.command == "verify" else regenerate(args.output_dir)
    except (VerificationError, OSError, ValueError) as exc:
        print(json.dumps({"ok":False,"error":str(exc)},sort_keys=True))
        return 1
    print(json.dumps(result,sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
