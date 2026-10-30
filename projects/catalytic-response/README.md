# A sharp response threshold for a catalytic three-species model

For the exact six-rate conserved catalytic model in Supplementary Figure 14B of
Floyd, Dinner, Murugan and Vaikuntanathan (2025), every positive input has one
physical steady state. Its B response is strictly increasing exactly when

\[
ad\le (a+c)(b+c).
\]

Above this threshold it has exactly one nondegenerate global maximum. The proof
also gives the exact peak state, input and overshoot height. The threshold and
peak height are independent of the two input-edge rates e and f.

The source already demonstrates a nonmonotone response, and its author notebook
already computes symbolic quadratic steady-state branches. The contribution here
is the complete positive-rate classification and sharp threshold for this one
model. General nonlinear expressivity and the separate structural-compatibility
question remain unresolved. See source.md, problem.md, proof.md and literature.md.

This is an AI-assisted mathematical result with a fresh-context AI review,
subject to the exact receipt in provenance/review-receipt.json. That review is
not external human peer review, a formal proof-assistant certificate, an author
confirmation, or a priority claim. The frozen proof retains its historical
pre-review scope wording; the receipt records the later review decision.

## Reproduce

Use Python 3.12 or newer and the exact dependencies:

```sh
python -m pip install -r requirements.txt
python -B verify.py verify
python -B -O verify.py verify
```

The commands require no network, API or GPU after dependencies are installed.
Verification checks the exact manifest inventory and every input hash, runs two
separate checkers in external temporary directories, compares the public
certificates, and rechecks all project bytes. No cache or output is written
inside the project. Both workers run with Python optimization 0 even when the
outer wrapper uses `-O`, preserving scientific assertions. Python 3.12 is the CI
target; the local packaging check used Python 3.13.5. Local execution is not
represented as a Python 3.12 test.

```sh
python -B verify.py regenerate --output-dir /absolute/new/external-directory
```

Regeneration creates a complete replayable project and new manifest in a new
exclusive external directory. It cannot overwrite the archive or an existing
directory. Extra files, missing inputs, symlinks, byte tampering and malformed
manifest inventories cause verification to fail.

The owner checker is byte-identical to the frozen research script and executes
all ten symbolic identities, the exact rational witness, and 2,000 supplementary
floating checks. Its nonzero exit status fails verification. The wrapper also
explicitly requires finite, nonnegative residual/sensitivity observations below
the original strict tolerances. The canonical owner certificate retains all
exact identity and witness data. It records only PASS, seed, sample size, ranges
and tolerances for the floating checks. Raw floating maxima and the count of
noncancelling sign comparisons are omitted because their bytes need not agree
across platforms. This projection is documented in provenance/source-mapping.json.
The theorem relies on proof.md, not on floating tests or the hash manifest.

No third-party PDF, author notebook, private correspondence or copied author
code is included. MIT licensing covers this packet's original code and
exposition; the cited paper and its authors' code remain their authors' work.
