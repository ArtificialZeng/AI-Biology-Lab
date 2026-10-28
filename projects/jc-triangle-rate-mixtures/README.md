# Shared site rates in a JC triangle: a compatibility dichotomy

For a fixed three-taxon Jukes–Cantor triangle network, arbitrary positive scalar
rate variation across sites has a useful exact consequence: if the resulting
distribution fits any **strict homogeneous triangle orientation**, it also fits
the true orientation. Alternatively, it fits none of the three strict images.

This note extends the homogeneous image description in
[Currie et al., arXiv:2606.26673v1](https://arxiv.org/html/2606.26673v1),
Theorem 3.1. It is a new restricted robustness question, not a correction to that
theorem or a claim to solve an author-declared open conjecture.

## Result and scope

Keep the six generating edge lengths positive and finite, the JC root law
uniform, and the generating inheritance probability fixed in `(0,1)`. At a
site, multiply every edge length by the same positive, almost surely finite
random rate. Its probability law can be arbitrary; no moment assumption is
needed. Write the mixed Fourier coordinates as
`(X,Y,Z,T)=(q_ACC,q_CAC,q_CCA,q_CGT)` and set `D=X*Y*Z-T*T`.

| Condition | Exact homogeneous compatibility |
| --- | --- |
| `D > 0` | The true orientation is compatible; other orientations may also fit. |
| `D <= 0` | None of the three **strict** homogeneous triangle images fits. |

The fitted homogeneous edge lengths and inheritance probability may differ
from the generating ones. The result does not ensure a unique orientation or
accurate parameter estimates. Boundary models, finite-data inference,
misspecified likelihood projection, lineage-specific rate changes,
rate-dependent inheritance, coalescence, larger networks and other substitution
models are outside its scope.

The proof preserves the orientation-specific source inequalities by a standard
monotone covariance identity. Its key model-specific step is that
`q_CAC(r)+q_CCA(r)-q_CGT(r)` decreases with the shared rate. The only source
membership condition that can fail after mixing is the common polynomial `D`.
See [the complete proof](proof.md) and [sources and prior-work scope](sources.md).

## Exact checks

Run from this directory with Python 3.9 or later; only the standard library is
used:

```sh
python verify.py
python -O verify.py
```

The read-only checker validates archive hashes and five finite examples with
12 rate categories. It sums the original internal ancestral states for all
64 DNA patterns, checks normalization and positivity, compares every Fourier
and inverse coordinate, checks displayed-tree path contractions, evaluates all
three complete strict image predicates, and checks the covariance identities.
Checks remain active under `-O`. The exact inputs and outputs are in
[witnesses.json](witnesses.json) and [certificate.json](certificate.json).
These are finite implementation checks; the arbitrary-law theorem is the
analytic argument in `proof.md`. The 64 entries have JC symmetries and only four
independent nonconstant Fourier coordinates.

A nonclosure example takes all six edge decays `1/2`, inheritance `1/2`, and
equally weighted rates 1 and 2. It gives
`(X,Y,Z,T)=(9/128,53/1024,53/1024,33/2048)` and
`D=-9567/134217728`. Thus even lawful scalar rate variation can prevent a strict
homogeneous triangle fit. A separate three-rate example has `D>0` and fits only
the true orientation; another fits all three.

Developed with substantial AI assistance, including mathematical review and a
separate reconstruction of the probability calculations in other AI contexts.
Independent human confirmation and exhaustive priority assessment remain
pending. This is a modest model-specific proof/code note, with no claim about
general estimator robustness.

Shawn Zeng (Zijian Zeng), UCSI University.
