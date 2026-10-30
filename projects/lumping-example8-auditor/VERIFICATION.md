# Verification and archive contract

`verify.py` differentiates the source's three coordinate equations. `independent_verify.py` separately reconstructs the six-species mass-action field from reactions (27)–(31), applies the full coordinate map, and differentiates that result. Neither imports or runs the other. The second implementation also checks the positive coefficients of the cubic stability polynomials: `A1`, `A2`, `A3`, and `A1*A2-A3` have respectively 8, 20, 16 and 96 positive-coefficient monomials, with strictly positive parts when `a=b=c=d=0` and `e,f,g,h>0`.

Both checkers use explicit exceptions, so checks remain active under `python -O`. Their default operation prints one JSON document and does not alter reports or archives. The only explicit report-writing interface is `--output PATH`, which uses exclusive creation and refuses overwrite. Run with `-B` to suppress Python bytecode creation as in the documented commands. Ordinary interpreter or dependency initialization is outside the application-level file-writing claim.

The four reports under `archive/` were generated locally on 26 October 2026 with Python 3.13 and pinned SymPy 1.13.3. Each implementation was run with and without `-O`; the respective reports agree byte for byte. Their exact hashes and hashes of the checking/proof inputs are recorded in `archive_manifest.json`.

From the repository root:

```sh
python -m pip install -r projects/lumping-example8-auditor/requirements.txt
python -B projects/lumping-example8-auditor/replay.py
```

The replay reads the manifest, checks every bound file, executes the four ordinary/optimized runs with bytecode writes disabled, and compares stdout bytes with the frozen reports. It rechecks hashes afterward and emits one JSON summary. It has no regeneration mode and does not rewrite the archive or manifest. A deliberate input change requires review, explicit regeneration of reports and manifest, and a corresponding explanation; updating hashes alone does not verify new mathematics.

CI is configured for Linux and macOS with Python 3.12 and 3.13, installing SymPy 1.13.3. This statement describes the workflow configuration; remote run success must be assessed from the actual CI results. Local verification here was performed on macOS/Python 3.13. No other Python or platform run is inferred from it.

The scripts check exact algebra. The generator spectral proof and the optional long-time limit argument are supplied in [PROOF.md](PROOF.md) and have received a separate AI mathematical review; they are not encoded in a proof assistant. Passing a replay or CI job does not establish novelty, priority, external peer review or upstream acceptance.
