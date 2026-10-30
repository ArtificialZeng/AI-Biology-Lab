# A common support need not sustain a switching population

This fixed positive-rate reaction network is a counterexample to the literal
Conjecture 43 of Cappelletti, Howells and Xu, arXiv:2507.09421v2. Both environments
admit nonnegative left vectors with positive drift on their supports. The
supports share species B, and the initial state contains one B. Nevertheless,
for every switching rate kappa>0 the population becomes extinct in finite time
almost surely.

| Environment 1 | Environment 2 |
|---|---|
| A -> 2A, rate constant 1 | A -> 0, rate constant 2 |
| B -> A, rate constant 1 | B -> C, rate constant 1 |
| C -> 0, rate constant 2 | C -> 2C, rate constant 1 |

The independent environment switches symmetrically at rate kappa. All six
existing rate constants are positive; no inflow reactions are present. The
population chain is reducible, which the literal source assumptions permit.
The source requires the environment chain to be irreducible.

The mechanism is that B is consumed once and never regenerated. Its descendants
remain A or C. Conditional on the environment path, their logarithmic mean growth
rate is -1/2, and the absorbing zero population then gives almost-sure extinction.
The source defines evanescence through positive probability of escape to infinity;
merely observing that the initial B state is transient would not refute it.

The source's Example 32 already supplies the disjoint-support motivation and
related switching birth/death behavior. Those facts and classical branching
methods are credited. The specific new contribution is the shared transient
feeder that meets the additional intersection hypothesis and defeats Conjecture
43 as stated. No priority or new general branching theorem is claimed.

Read [proof.md](proof.md) for the complete argument and [source.md](source.md) for
the precise source contract. The separate
[Example 44 / Conjecture 46 proof](../switching-transience/README.md) concerns a
regenerative common species and is retained as a distinct result in the same
research PR. Neither model is silently substituted for the other.

## Reproduce

Python 3.12 or newer, standard library only:

```sh
python -B verify.py verify
python -B -O verify.py verify
python -B verify.py regenerate --output-dir /absolute/new/external-directory
```

Verification checks the complete archive inventory and hashes, runs the frozen
owner and independent checkers in external temporary directories, compares their
full deterministic JSON outputs, and rechecks unchanged inputs. Workers always
run with Python optimization 0, including when the outer wrapper uses -O. The
independent checker uses assertions and is invoked in that supported mode.
Regeneration requires an exclusive new external directory and makes a complete
replayable packet. Extra/missing files, symlinks and byte tampering fail.

The owner checker performs 250 exact state/environment cases and six mutation
controls. The separate checker has 432 exact generator-coordinate checks and
additional lineage/occupation identities. These finite computations do not
formalize or prove the probability theorem. The complete extinction proof was
reviewed in a fresh AI context; provenance/review-receipt.json records that scope.
Human author confirmation and priority remain pending. No formal runner gate,
proof-assistant certification, manuscript/PDF or biological experiment is claimed.
Original code and exposition use the included MIT license. No source paper,
private correspondence, paid service, API key or GPU is included or required.
