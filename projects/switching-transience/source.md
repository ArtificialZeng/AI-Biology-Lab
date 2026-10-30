# Exact primary model and question

Daniele Cappelletti, Aidan S. Howells and Chuang Xu, *Stability of randomly
switching stochastic reaction networks with asymptotically linear transition
rates*, [arXiv:2507.09421v2](https://arxiv.org/abs/2507.09421v2), revised 12 August 2025.
The current source was checked on 26October2026; printed pp. 29–30, Example 44
and Conjecture 46, are controlling. The reaction diagram was visually checked.

With species(S1,S2,S3), alpha>0 and mass-action propensities:

| Environment1 | Environment2 |
|---|---|
| S1 -> 4S1+S2, rate constant 1 | S3 -> 4S3+S2, rate constant 1 |
| S2 -> S1, rate constant 1 | S2 -> S3, rate constant 1 |
| S3 -> 0, rate constant alpha | S1 -> 0, rate constant alpha |
| 0 -> S3, constant rate 1 | 0 -> S1, constant rate 1 |

The independent environment jumps in either direction at rate kappa>0.
model.json records every arrow and inflow. The source's Conjecture 46 asks for
some alpha>0 such that this chain is transient for all sufficiently smallkappa.
The note proves this exact claim with the stronger uniform range stated in
proof.md. It retains the source's inflows and changing supports; no deterministic
mean approximation replaces the stochastic process.

SourceTheorem 16 handles identical supports of positive drift directions.
Example 44 has only the shared S2 support. Proposition 45 is conditional on the
broader Conjecture 43; Proposition 47 rules out one Lyapunov candidate, not the
conjecture itself. General Conjecture 43 is not resolved here.
