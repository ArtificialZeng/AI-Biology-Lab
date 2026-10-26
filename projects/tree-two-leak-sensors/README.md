# Where should one extra sensor go in a two-leak tree?

Consider any finite bidirected tree with independent positive edge rates and
two positive leaks at distinct known vertices a,b. The sole input and one
calibrated output are retained at r; add a calibrated output at `s != r`.

**The model is generically locally identifiable exactly when the path from r
to s shares at least one edge with the path from a to b.** Sharing only one
vertex is insufficient. The [proof](PROOF.md) covers arbitrary finite trees,
including coincidences between a leak, the input, and a sensor.

The root-only model has a known one-dimensional parameter ambiguity. An
explicit positive diagonal-similarity curve describes that ambiguity. An
outgoing-rate harmonic voltage varies strictly along the leak-to-leak path and
is constant on its hanging branches. The new sensor removes the missing
direction exactly when it observes a different voltage from the root.

With the root measurement retained and sensor placement free, **one additional
sensor is both necessary and sufficient** for this fixed two-leak model.
This is a structural, noiseless, generic local statement. It does not guarantee
global uniqueness, full rank at every positive point, or reliable estimation
from noisy data. Unknown sensor gains and unknown leak locations are outside
the stated model.

## Reproduce

From the lab repository root:

```sh
python -m pip install -r projects/tree-two-leak-sensors/requirements.txt
python projects/tree-two-leak-sensors/verify.py
python -O projects/tree-two-leak-sensors/verify.py
python projects/tree-two-leak-sensors/independent_verify.py
```

The first checker tests all 3,464 root/leak/sensor placements on nonisomorphic
trees through n=6 over an exact finite field. Those ranks give lower-bound
certificates, supplemented by rational rank/minor and positive-ambiguity
witnesses. The independent standard-library implementation uses rational
resolvent derivatives, checking 164 small-tree configurations and five
seven-vertex controls. Neither finite test proves the arbitrary-tree theorem;
that is the role of the written argument and its cited dependencies.

Both scripts print receipts; `--output PATH` saves a separate receipt.
[Sources](SOURCE.md) · [scope](SCOPE.md) · [proof](PROOF.md) · [verification](VERIFICATION.md)

The exact mathematical scope passed a separate AI review. Priority and external
human confirmation remain pending. The broader minimal-input/output questions
in the source review remain unresolved. [Credit and permissions](THIRD_PARTY.md).
