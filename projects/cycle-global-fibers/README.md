# Exact ambiguity counts for two-leak directed cycles

For a directed n-cycle, input 1, output p, and two fixed leaks at r,t, assume
`n>=3`, `2<=p<n`, and `1<=r<p<=t<=n`. The existing local-identifiability
classification is due to [Ahmed et al.](https://arxiv.org/html/2412.05283v3).
This note addresses the remaining global ambiguity in that family.

On the nonzero, simple-root data domain stated in [SCOPE.md](SCOPE.md), the
complete affine complex parameter fiber has

```
p! (n-p)!
```

distinct points, with explicit inverse branches. Complex leak coordinates may
be zero. On the all-nonzero parameter torus the generic degree is the same,
while a special fiber can lose zero-leak branches.

For data from strictly positive rates, every affine complex branch is real and
has positive edge rates. Positivity of the two leaks reduces to a single root
interval. With the recoverable invariants defined in [PROOF.md](PROOF.md),

```
m = number of negated U-roots in ( C/(kappa B), A/kappa )
positive parameter vectors = m (p-1)! (n-p)!
```

Every value `m=1,...,p` occurs on a nonempty open positive parameter region.
For example, `n=4,p=3` has complex degree 6 but can have 2, 4, or 6 positive
solutions. This explains how positivity filters algebraic ambiguity without
claiming that an entire model family becomes globally identifiable.

## Reproduce

From the lab repository root:

```sh
python -m pip install -r projects/cycle-global-fibers/requirements.txt
python projects/cycle-global-fibers/verify.py --self-test
python -O projects/cycle-global-fibers/verify.py --self-test
python projects/cycle-global-fibers/independent_verify.py
```

The first checker covers 392 rational fixtures and 80,576 inverse branches,
plus 25 symbolic state-space/cofactor models. Three deliberately wrong inverse
branches must be rejected. The independent standard-library implementation
checks 182 symbolic models through n=8, 182 rational fixtures and 8,864 branches,
including zero-leak endpoint and cross-factor repeated-root controls. These
finite checks supplement the all-n algebraic bijection and multiplicity proof.
Both scripts print receipts; `--output PATH` saves a separate receipt.

[Sources](SOURCE.md) · [scope](SCOPE.md) · [proof](PROOF.md) · [verification](VERIFICATION.md)

The mathematical scope passed a separate AI review. Novelty and priority remain
unconfirmed by the source authors or external human reviewers. The source's
local classification and earlier zero/one-leak degree results are credited,
not claimed as new. See [third-party credit](THIRD_PARTY.md).
