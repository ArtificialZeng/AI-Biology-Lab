# Periodic forcing changes the mean of a reversible adapting motif

For the two-state mass-action signalling model

```
x' = (a+c)u(t) - bx - dxy
y' = cu(t) - dxy
```

all a,b,c,d>0 and every positive continuous T-periodic input have one positive
T-periodic response. It attracts every nonnegative initial state exponentially.
For constant input, the adapted output is y0=bc/(ad). On the periodic response,

```
y0 - mean(y) = h * mean((y-y0)^2 / z),
h=(a+c)/c,   z=x-h(y-y0)>0.
```

Thus every nonconstant input produces a strictly lower mean output than the
constant-input baseline. Equality occurs exactly for constant inputs. The proof
also gives the second-order loss for a small sinusoidal input and the frequency
maximizing that leading coefficient. It does not optimize finite amplitudes.

This is the reversible motif in Franco and Velázquez's published Example3.8,
not a newly designed circuit. Constant-input adaptation, generic failures of
periodic mean adaptation, and contraction methods are prior work. The narrow
contribution is an exact signed all-waveform identity/equality classification
and a complete explicit analysis of this source example. Its scientific scale
is a modest theoretical result, with priority and human verification pending.
It does not resolve the source's general periodic-input question or establish
dependence on the average input alone.

Read [proof.md](proof.md), [source.md](source.md), [problem.md](problem.md) and
[literature.md](literature.md). The published reaction labels imply that source
species2 is the adapted output; the source's prose says species3. This indexing
correction alone is not counted as a research contribution.

## Reproduce

Use Python3.12 or newer and the pinned exact-algebra dependencies:

```sh
python -m pip install -r requirements.txt
python -B verify.py verify
python -B -O verify.py verify
python -B verify.py regenerate --output-dir /absolute/new/external-directory
```

Verification checks the complete archive inventory and hashes, runs both
checkers in external temporary directories, and compares full deterministic
certificate bytes. Workers run with optimization0 even under the optimized
outer wrapper. Regeneration creates a new complete replayable packet and
requires an exclusive external absolute directory. Inputs may be read-only.

The owner checker has16 symbolic identities,2025 exact reaction/state cases,
648 invariant-polygon checks and six rejected mutation controls. The separate
checker reconstructs stoichiometry and19 symbolic identities. These checks
do not formalize Brouwer's theorem, ODE existence, contraction or the mean
theorem; those steps are justified in the written proof and AI review.

Substantial AI assistance and fresh-context AI review are disclosed in
[provenance/review-receipt.json](provenance/review-receipt.json). No human peer
review, author endorsement, proof-assistant certificate, authenticated runner
release gate or priority is claimed. No third-party paper, private mail,
manuscript or PDF is distributed. No GPU, API or paid computation is needed.
