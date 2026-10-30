# Extinction hypotheses and a corrected first-moment theorem

This packet studies the general one-type-size-dependent multitype branching model of Mehra–Taylor (2026). It contains two exact hypothesis witnesses and an analytic exponential-extinction-moment theorem, pending author verification.

* A two-type example satisfies literal A/B/C and has a Hurwitz limiting mean matrix, yet a three-state positive closed class prevents extinction.
* A separate example has accessible zero and almost-sure extinction, but infinite mean extinction time. Its limiting matrix is critical despite pointwise subcriticality above the carrying capacity.
* With uniformly bounded first offspring moments, bounded positive individual rates, uniformly positive controller-production probability, an explicitly Hurwitz limiting mean matrix, and accessibility of zero, there are theta>0 and C<infinity such that E_x exp(theta*T0)<=C*(1+|x|). Offspring support may be unbounded and its second moment infinite.

The last theorem is in [first-moment-proof.md](first-moment-proof.md). Its corollary covers source A and either reading of B after explicitly imposing the Hurwitz-limit and accessibility conditions, **without bounded-support Assumption C**. The source's uniform second moments convert the displayed mean version of B into the probability version by Cauchy–Schwarz. The original unrepaired A/B conjecture is not claimed proved unconditionally, and the examples do not refute the paper's particular immune model.

The earlier bounded-offspring proof and both witnesses remain unchanged in [proof.md](proof.md). This is one coherent updated research packet, not a separate paper or additional discovery count. The model-specific cutoff weight V=R^(-min(x1,Y))*(xi dot x) is retained. The extension uses expectation bounds in place of upper jump bounds and explicitly justifies nonexplosion, unbounded overshoots, stopped Dynkin, and finite-set absorption.

Run from the repository root with Python 3.12 or later and only the standard library:

```sh
python projects/controlled-branching-extinction-audit/reproduce.py
python -O projects/controlled-branching-extinction-audit/reproduce.py
```

The checker reconstructs the original closed-class generator, checks rational/polynomial certificates, and checks 100 finite occupation recurrences. It also verifies packet integrity. **It does not verify either infinite-state theorem, infinite mean, or the heavy-tailed illustration.** Those claims use the analytic arguments in the two proof notes. [sources.md](sources.md) records the source contract and classical methods.

Development and checking used substantial AI assistance. Separate fresh AI mathematical and source/value reviews are distinct checks, not human peer review, author endorsement, or priority certification. No experimental biological validation, journal manuscript, or submission is included.
