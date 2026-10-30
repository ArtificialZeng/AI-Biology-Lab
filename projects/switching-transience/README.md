# Shared-species branching under random environment switching

This original proof note addresses Example 44/Conjecture 46 of Cappelletti,
Howells and Xu, arXiv:2507.09421v2. For the exact three-species process,
it establishes transience for every alpha>0 and every 0<kappa<=1, with total
population tending to infinity almost surely from every finite initial state.
This sufficient switching interval is not claimed sharp. This proof concerns
Example 44 only. The separate [common-support counterexample](../switching-common-support/README.md)
refutes the literal general Conjecture 43 with a different six-reaction model.
The fixed scientific certificates here describe Example 44; their earlier scope
label about Conjecture 43 is superseded by that separate proof note.

The mechanism is a labelled subpopulation retaining only the common species S2
at environment switches. Conditional on IID holding times it is a one-type
branching process. A self-contained conditional second-moment argument proves
positive survival; within-interval population domination then yields escape of
the full process. Mean growth alone would not prove this conclusion.

Read proof.md for the probability argument, source.md for the exact original
model, and literature.md for prior results and the bounded novelty screen.
The algebra/path certificates are reproducibility checks, not formalizations of
the stochastic proof. Substantial AI assistance and a fresh-context AI review
are disclosed in provenance/review-receipt.json. No human author confirmation,
proof-assistant certificate, trusted-runner release gate or priority is claimed.
No manuscript, PDF, author code, source paper or private correspondence is included.

## Reproduce

Use Python 3.12 or newer with pinned requirements:

```sh
python -m pip install -r requirements.txt
python -B verify.py verify
python -B -O verify.py verify
```

Both modes verify the complete input inventory and hashes, run the frozen owner
checker and separate reviewer checker in external temporary directories, compare
exact output bytes, and confirm that the project remains unchanged. Workers
always run with Python optimization 0, including when the outer wrapper uses -O.
No network, API or GPU is needed after dependencies are installed.

```sh
python -B verify.py regenerate --output-dir /absolute/new/external-directory
```

This creates a full replayable copy in a new exclusive external directory.
Existing/source/inside/symlink destinations, extra or missing files, symlinks and
byte tampering are rejected. The manifest protects provenance, not mathematical
truth. The owner's 46 symbolic identities and 128 small reaction paths do not prove
universal transience; that requires the complete argument and its hypotheses.
Original code and exposition in this packet use the included MIT license.
