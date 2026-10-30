# Claim and limitations

The audited invariant concerns models defined by `Sv=0` and reaction interval bounds. Negating selected columns of `S` and simultaneously replacing `[l,u]` by `[-u,-l]` is a coordinate change preserving zero/nonzero flux status. If a model has additional user constraints involving those flux coordinates, those constraints must also be transformed before asserting equivalence.

For a fixed negative-only reaction `l <= v <= u <= 0`, `w=-v` is its feasible nonnegative detection quantity. The repaired auxiliary LP uses `w >= z >= 0`. Reversible and positive-only reaction expressions remain unchanged; fixed-zero fluxes remain zero. The source patch does not change the original constraints.

The independent five-test suite has one deterministic baseline failure: a strict-negative pathway produces a solver-infeasible warning promoted to a test failure. Four controls pass. The candidate passes all five tests without warnings. The positive/strict-negative cases include an exact `Fraction` witness. Dead-end and fixed-zero balances identify the blocked set independently of the optimization backend.

The separate 24-case exact oracle enumerates all eight pathway orientations, open/closed transfer, and positive-minimum path bounds. It checks integer/rational rank, nullspaces, mass balance and feasible witnesses. It verifies this finite family; it is not a theorem about every metabolic network.

The earlier zero-allowed all-negative example can sometimes pass accidentally because later LP solutions activate unrequested pathway reactions. Consequently, a previously observed eight-failure suite is not used as a required portable count. The archived independent 1-failure/4-control suite uses strict-negative bounds instead.

Tests preserve reaction stoichiometry and bounds, the objective, existing solver variables and constraint expressions/bounds, and nested context depth. This preservation check includes an unrelated dummy solver variable/constraint. It does not assert sign-coordinate equivalence for arbitrary untransformed user constraints.

Coverage and status:

- GLPK 5.0 through `swiglpk` 5.0.13 and Optlang 1.9.1; alternate backends were not checked here.
- Five independent regression tests and 24 exact finite cases; no large-model performance or benchmark claim.
- Existing FASTCC `LARGE_VALUE`, numerical tolerance and infeasible-input behavior remain outside this repair.
- Separate AI review accepted the bounded software scope. This is not external human peer review, maintainer acceptance, or proof of universal FASTCC correctness.
- No first-discovery, new-algorithm, new-theorem or new-biological-mechanism claim. Related PR1495 is credited.
