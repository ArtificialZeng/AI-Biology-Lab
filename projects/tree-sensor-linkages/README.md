# How many calibrated sensors remove parameter ambiguity in a leaking tree?

For any finite bidirected tree, let **B** be its nonempty set of known leak locations, with **L** independent positive leak rates. The single input and one calibrated observation are retained at **r**. Let **O** contain all observed vertices, including r. All directed edge rates are independent and positive.

Let **q** be the maximum number of pairwise **vertex-disjoint paths** joining distinct observations to distinct leaks. Endpoints count; a path can have length zero.

The [theorem and proof](PROOF.md) give:

- generic transfer-data rank **2n − 2 + q**;
- generic parameter ambiguity dimension **L − q**;
- generic local structural identifiability exactly when **q = L**;
- **L − 1 additional calibrated observations** are necessary and sufficient when their locations may be chosen freely.

The underlying harmonic evaluation matrix has rank q at **every positive** edge-rate assignment. If q<L, an explicit (L−q)-dimensional positive family preserves every observed transfer function through every positive parameter point. The full physical rank statement is **generic**; exceptional points can have larger fibers. Global uniqueness and statistical conditioning are outside scope.

This extends the lab's [two-leak sensor result](../tree-two-leak-sensors/README.md), which is now a corollary of the same contribution. The known root-only rank theorem and the classical all-minors matrix-tree identity are dependencies. The bounded [literature comparison](LITERATURE.md) found no equivalent complete arbitrary-tree placement theorem in the inspected sources; **historical priority and external human review remain pending**. The broader multiple-input/output problem is unresolved.

## Example

In the [seven-vertex branching-tree certificate](evidence/certificate.json), leaks are at vertices 0,2,4,6 and the input/root sensor is at 1. Observations {1}, {1,3}, {1,3,5}, {0,1,3,5} have q=1,2,3,4 respectively. Thus the generic ambiguity dimensions are 3,2,1,0. Sensor count alone is insufficient: other placements of the same count can leave fewer disjoint observation-to-leak paths.

## Reproduce

From the repository root:

```sh
python -m pip install -r projects/tree-sensor-linkages/requirements.txt
python projects/tree-sensor-linkages/verify.py
python -O projects/tree-sensor-linkages/verify.py
python projects/tree-sensor-linkages/independent_verify.py
python projects/tree-sensor-linkages/replay.py
```

Ordinary commands print JSON without changing archived evidence. Use `--output artifacts/tree-linkages` to save a new primary certificate explicitly. The manifest/replay script compares regenerated results to the archive and rejects changed mathematical inputs.

The primary checker tests 204,625 harmonic placements on all unlabelled trees through n=7, 8,497 independently constructed physical Markov Jacobians through n=5, and 100 held-out labelled trees through n=16, using exact modular arithmetic. Finite checks support the proof; they do not establish its infinite quantifier or historical novelty. See [verification details](VERIFICATION.md), [model and scope](SCOPE.md), and [primary sources](SOURCE.md).

No biological experiment, unknown sensor gain, noisy-data estimation guarantee, or physical virtual-cell validation is claimed. Separate AI review is recorded as an automated mathematical check, not human peer review.
