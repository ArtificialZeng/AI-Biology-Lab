# Local damping in a six-species chemical clock

A small, source-specific code/proof note for the six-species mass-action clock in Hemery (2026). The nonhyperbolic uniform equilibrium is locally asymptotically stable. Exact algebra gives the center cubic coefficient `-3/14 + i*sqrt(3)/7`; an explicit polynomial Lyapunov function certifies a small positive trapping neighborhood.

The original H6 question—whether an attracting positive periodic orbit with six ordered peaks exists elsewhere—remains unresolved. This package does not certify the source initial trajectory, global convergence, the converter, or all clock sizes. The source reports numerical oscillations and does not assert an attracting cycle.

## Verify

Python 3.12 or later is recommended. Install the one pinned dependency, then run:

```sh
python -m pip install -r projects/six-clock-local-damping/requirements.txt
python projects/six-clock-local-damping/verify.py verify
python -O projects/six-clock-local-damping/verify.py verify
```

Both commands check the immutable manifest, recompute both certificates in fresh temporary directories, compare exact JSON objects, and leave existing evidence untouched. The wrapper re-executes with `PYTHONOPTIMIZE=0`, so scientific assertions stay enabled even when called with `-O`.

To retain newly generated certificates, choose a **new directory outside this checkout**:

```sh
python projects/six-clock-local-damping/verify.py regenerate --output /tmp/six-clock-new-certificates
```

The destination must not already exist. Regeneration never replaces the archived certificates or manifest. The proof is in [proof.md](proof.md); scope, attribution, and review limits are in [STATEMENT.md](STATEMENT.md), [SOURCES.md](SOURCES.md), and [REVIEW-SUMMARY.md](REVIEW-SUMMARY.md).

`local_algebra.py` is the original exact calculation. `independent_checker.py` reconstructs the physical reaction field and projects its real inverse DFT, instead of importing the original calculation. Its public-path adaptation and package integration are writing work, not a new independent mathematical review. No manuscript or PDF is included.
