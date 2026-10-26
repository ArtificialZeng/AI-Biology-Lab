# Adjacent-port catenary singular-locus auditor

For a positive bidirected chain with one fixed leak, input at node1 and calibrated output at node2, this package supplies exact exceptional equations for the polynomial coefficient Jacobian. It distinguishes that map from the measured transfer-function differential at pole-zero cancellation.

The original four-node question covers every fixed single-leak location. A separate theorem proves the same characterization for every chain length n>=2. See [statement.md](statement.md), the complete [proof.md](proof.md), [source.md](source.md), and [prior comparison](evidence/prior-check.md). A fresh substantive AI mathematical review accepted the stated scope; see [review-summary.md](review-summary.md). Historical priority is unestablished.

Install Python3.12 or3.13 and the pinned dependency, then run:

```sh
python -m pip install -r requirements.txt
python -B verify.py verify
python -B verify_independent.py verify
python -O -B verify.py verify
python -O -B verify_independent.py verify
```

Each command prints JSON, hashes every public input listed by archive_manifest.json, recomputes its own certificate, compares exact bytes, and leaves the archive unchanged. No assertions are used for mathematical or I/O checks, so optimization preserves checking. The first builder computes all square coefficient Jacobians for n2 through4 and all leak locations. The second independently constructs the numerator by a cofactor and characteristic polynomials by a recurrence, checking19exact cases: every n4 leak position generic/singular, an input-leak cancellation example with full coefficient rank, and generic/singular examples for every n5 leak position. These finite checks do not prove the arbitrary-length theorem.

To emit deterministic certificate JSON without writing files:

```sh
python -B verify.py compute
python -B verify_independent.py compute
```

Explicit regeneration computes both certificates into a new absolute directory outside this project:

```sh
python -B verify.py regenerate --output-dir /absolute/new/adjacent-catenary-check
```

The target must not exist and its parent must already exist. Relative, existing, symlink, and project-internal output paths are rejected. Regeneration writes only the two certificates and a candidate manifest into that exclusive output directory. It never replaces the archived files. The regenerated directory is an evidence bundle rather than a full project copy; comparing and installing any changed archive is a separate review action.

All rates and the sole leak are strictly positive and independently unknown. Port gains are known. The polynomial coefficient map retains its original coordinates at cancellation points, whereas zero-initial-condition transfer data may have smaller reduced degree; proof.md gives both differential-rank contracts. No boundary, global fiber-count, or human-priority claim is made.

Original material is licensed under MIT. Linked works retain their rights; see THIRD_PARTY.md.
