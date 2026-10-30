# Global exclusion in a density-triggered microbial transfer model

This original code/proof packet proposes an affirmative answer to the strict two-environment global-exclusion question in Schreiber (2025), Section 5.2, equation (26). The mathematical statement concerns an existing model of serially transferred competing microbial populations; it is not a new experiment or an empirical claim about antibiotic resistance.

## Result

Write a_j,b_j>0 for the two species' growth rates in environment j and 0<delta<1 for their common dilution. Each growth phase stops when the total population reaches one; immediately after dilution the total is delta.

For the source alternating environments, if

    (a_0/b_0+a_1/b_1)/2 < 1,

then species 2 excludes species 1 from **every positive post-transfer composition**, for either starting environment. At complete two-transfer periods,

    odds_1(2n) <= odds_1(0) exp(-n gamma),
    gamma=log(1/delta) [2-a_0/b_0-a_1/b_1]>0.

The argument uses a composition-uniform bound for each growth phase, not a linearization near a resident-only state. The source already derives the invasion rates and local exclusion; those are credited rather than claimed as new.

The same proof gives a coherent finite-period extension. With

    U=log(1/delta) sum_j(a_j/b_j-1),
    V=log(1/delta) sum_j(1-b_j/a_j),

all nonneutral finite periodic schedules have the following global alternatives:

| Condition | Interior dynamics |
|---|---|
| U<=0 | Species 2 excludes species 1 |
| V>=0 | Species 1 excludes species 2 |
| V<0<U | A unique interior phase cycle attracts every positive composition |

If every a_j=b_j, all compositions are preserved. Threshold equalities are included in the exclusion result, with no claimed uniform exponential rate at equality. With distinct phase labels the augmented cycle has the prescribed period; only its population-composition block can have a shorter period when phases repeat.

## Read and reproduce

Read [the precise question](problem.md), [full analytic proof](proof.md), and [primary-source comparison](source.md). The all-parameter conclusions rest on the proof, not on finite simulations.

```sh
python -m pip install -r projects/serial-transfer-global-dynamics/requirements.txt
python -B projects/serial-transfer-global-dynamics/verify.py verify
```

The checker compares archived certificates and hashes without writing project files. It includes symbolic derivative identities, 450 exact rational transfer witnesses, and separately labelled floating-point examples. A second implementation independently reconstructs the threshold equations, an exact nonneutral coexistence cycle and finite numerical attacks. CI compares normal and optimized **wrappers**; both deliberately execute the scientific children with assertions enabled (`PYTHONOPTIMIZE=0`). Direct optimized execution of the independent checker is not a supported verification route. To regenerate an owner's certificate deliberately, run `science/owner_checker.py` and redirect its JSON output outside the project first; changing evidence requires renewed comparison and review, not merely new hashes.

## Scope and verification status

The exact public proof received a separate AI mathematical check, including an independent derivation of the transfer derivative and source-model comparison. The checked scope includes the strict source question, equality cases, finite-period classification and continuous-time setwise convergence. [The receipt](provenance/review-receipt.json) binds the proof and checker hashes and records the limitations. This is no author confirmation, human peer review, priority certification or proof-assistant formalization. The trusted-runner gate did not pass because genuine runner job IDs are unavailable; no manuscript/PDF or formal runner release is included. This original code/proof packet is provided for technical review.

Exactly two species, positive constant growth within each phase, a common dilution and a fixed shared density threshold are essential assumptions. The environment schedule is prescribed. Multispecies competition, density-dependent growth, unequal dilution, demographic randomness and feedback-dependent environment switching are outside the theorem. For nonperiodic schedules the proof states only a cumulative odds envelope, not a general classification.

After every transfer the surviving density approaches delta in exclusion cases. Between transfers it continues to grow and be diluted: the continuous-time conclusion is approach to the boundary flow-kick cycle, not convergence to a fixed density or synchronization to an absolute clock.

Nearby literature already contains global alternatives for a different nutrient-triggered fermentation model. No claim is made to the first global classification of state-triggered microbial competition. See [source.md](source.md) for the inspected passages and the limits of the literature search.

AI assistance was used for derivation, programming and separate checking passes. The work is an independent computational project and does not claim institutional endorsement.
