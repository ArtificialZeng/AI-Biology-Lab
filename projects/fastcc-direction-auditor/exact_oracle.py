#!/usr/bin/env python3
"""Check 24 exact rational signed-path certificates; default stdout, no archive writes."""
import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify():
    # Columns: uptake, pathway, drain, dead end, fixed-zero parallel transfer.
    # Rows: A, B, C balances and the fixed-zero parallel-transfer equality.
    base = sp.Matrix([[1, -1, 0, -1, -1], [0, 1, -1, 0, 1],
                      [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])
    cases = []
    for mask, (minimum, close_path) in product(range(8), [(0, False), (0, True), (1, False)]):
        signs = [-1 if mask & (1 << j) else 1 for j in range(3)] + [1, 1]
        matrix = base * sp.diag(*signs)
        if close_path:
            matrix = matrix.col_join(sp.Matrix([[0, 1, 0, 0, 0]]))
        kernel = matrix.nullspace()
        require(matrix.rank() == (5 if close_path else 4), "Unexpected exact rank")
        witness = sp.zeros(5, 1) if close_path else sp.Matrix([10 * s for s in signs[:3]] + [0, 0])
        require(matrix * witness == sp.zeros(matrix.rows, 1), "Witness violates a balance")
        bounds = []
        for j, sign in enumerate(signs):
            lower, upper = (minimum, 10) if j < 3 else ((0, 10) if j == 3 else (0, 0))
            if close_path and j == 1:
                lower, upper = 0, 0
            if sign < 0:
                lower, upper = -upper, -lower
            bounds.append([lower, upper])
            require(lower <= witness[j] <= upper, "Witness violates a bound")
        require(not kernel if close_path else kernel[0][3] == kernel[0][4] == 0,
                "Unexpected blocked-reaction nullspace coordinates")
        # Verify the complete one-dimensional kernel is the signed pathway line.
        if not close_path:
            expected_line = sp.Matrix(signs[:3] + [0, 0])
            require(matrix * expected_line == sp.zeros(matrix.rows, 1), "Wrong pathway direction")
        cases.append({"orientation_mask": mask, "minimum_flux": minimum,
                      "closed_path": close_path, "matrix": matrix.tolist(),
                      "rank": matrix.rank(), "nullspace": [v.T.tolist()[0] for v in kernel],
                      "bounds": bounds, "witness": list(witness),
                      "consistent_reactions": [] if close_path else ["uptake", "path", "drain"]})
    return {"case_count": len(cases), "arithmetic": "exact integers and rationals",
            "scope": "Sv=0 and reaction interval bounds", "cases": cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="optional separate JSON output outside the project")
    args = parser.parse_args()
    manifest = json.loads((HERE / "archive_manifest.json").read_text())
    for relative, sha in manifest["sha256"].items():
        require(hashlib.sha256((HERE / relative).read_bytes()).hexdigest() == sha,
                "Archived file changed: " + relative)
    rendered = json.dumps(verify(), indent=2, default=str) + "\n"
    if args.output:
        path = args.output.resolve()
        require(path != HERE and HERE not in path.parents, "Preserve the project's archived evidence")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
