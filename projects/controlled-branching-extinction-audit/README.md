# Extinction hypotheses and a corrected exponential-moment theorem

This packet audits the general one-type-size-dependent multitype branching model of Mehra–Taylor (2026), and proves a bounded-offspring exponential-moment theorem under explicit additional hypotheses. It contains original code and a mathematical proof note, pending author verification.

* A two-type example satisfies the literal A/B/C assumptions and has a Hurwitz limiting mean matrix, yet a three-state positive closed class prevents extinction.
* A separate two-type example has accessible zero and almost-sure extinction, but infinite mean extinction time. Its limiting matrix is critical despite pointwise subcriticality above the carrying capacity.
* With bounded offspring, uniformly bounded positive individual rates, uniformly positive controller production by other types, an explicitly Hurwitz limit, and accessibility of zero, the proof establishes constants theta>0 and C<infinity with E_x exp(theta*T0)<=C*(1+|x|).

The examples address literal general hypotheses. They do not refute the intended class with both missing conditions imposed, or the source's particular immune model. The positive theorem retains bounded offspring; no unbounded-offspring or experimental claim is made.

Run from the repository root, using Python 3.12 or later and only the standard library:

```sh
python projects/controlled-branching-extinction-audit/reproduce.py
python -O projects/controlled-branching-extinction-audit/reproduce.py
```

The checker reconstructs the exact closed-class generator, checks rational/polynomial certificates, and checks 100 finite occupation recurrences. **These finite checks do not prove infinite-state extinction, infinite mean, or the general theorem.** Those conclusions require the analytic arguments in [proof.md](proof.md). [sources.md](sources.md) records the source contract and prior methods.

The model-specific ingredient is the cutoff weight V=R^(-min(x1,Y))*(xi dot x), which gives a geometric drift bound even when non-controller population is arbitrarily large. A finite-set return/absorption argument completes the corrected theorem.

Development and checking used substantial AI assistance. Separate fresh AI reviews are computational/mathematical checks, not human peer review, author endorsement, or a claim of priority. No journal manuscript or submission is included.
