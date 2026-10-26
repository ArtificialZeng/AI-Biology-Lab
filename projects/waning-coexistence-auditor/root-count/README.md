# Full-model root-count classification

**A fresh manual AI review accepted this mathematical scope on 26 October 2026.** It is not external human peer review, a priority certificate or a formal runner gate. This extends the earlier accepted full-waning existence result. The [precise statement](problem.md), [proof](evidence/proof.md), both mathematical algorithms and both new certificates are byte-preserved from the frozen research candidate. The [source target](source.md) and literature introduction have only process metadata removed for public presentation; their source facts and mathematical scope are unchanged. Its research snapshot digest is `4d08856e7e88c31d27c0d93ec4eea728d6b3314295bdabb6b412cf99b181cf64`; proof SHA-256 is `0fe7cfcb2e28d591cdcc7b5185d13b84a1e81ccf14ac68a987c9e4bbf6637a5a`.

For the full eight-compartment model with all finite positive biological parameters and finite nonnegative waning rates, relabel so `p<=q`. The proof establishes strict convexity of the reconstructed population on the feasible interval. A cubic obtained from the **same source definition** `PG=-Delta*FG` yields explicit semialgebraic tests for no root, one simple root, two simple roots, or one double root. It includes `p=q`, a linear Delta, P-degrees one/two/three and endpoint tangencies. A repeated feasible root is exactly double. Coincident critical thresholds have a separate endpoint interpretation. The broad abstract Conjecture 1 is not resolved.

The accepted mathematical scope gives a complete answer to actual-model Question 2. Historical candidate/pending wording inside the byte-preserved proof, statement and certificates records their pre-review state; this presentation and the parent review summary give the current status. Cubic identities and 24 exact held-out cases do not substitute for the analytic proof. Scalar multiplicity does not establish ODE stability, persistence or an eight-dimensional saddle-node. Prior permanent/no-ADE, one-way-zero and large-transmission uniqueness is credited in [the literature record](evidence/literature.md); final JTB full-text access and the latest detailed SIAM HTML refresh remain stated limitations. No historical priority is certified.

## Source normalization

The source Eq. (6) and Lemma 6 define `PG=-Delta*FG`. The population identity gives exactly `PG=mu*((s+r1+r2-s0)*Delta+A1*W1+A2*W2)`, so its degree is at most three and `P(d)=mu*(A1(d)*W1(d)+A2(d)*W2(d))>0`. Here d is the reconstruction pole, not a root of this defining polynomial. Later source prose describing a quartic or asserting `PG(d)=0` does not apply to that definition. The candidate keeps the source's actual polynomial and feasible interval; it does not remove an s-dependent factor or count the pole as an equilibrium.

## Reproduce from this directory

Use the parent's pinned SymPy 1.13.3 environment, or install with `python -m pip install -r ../requirements.txt`. Keep virtual environments and generated outputs outside the project.

```sh
python check.py verify
python -O check.py verify
python check_heldout.py verify
python -O check_heldout.py verify
```

The first wrapper invokes the unchanged generic symbolic algorithm in `verify_root_count.py`. The second invokes the unchanged independent algorithm in `evidence/heldout_check.py` and captures its JSON for comparison. These wrappers change archive I/O only; they do not change formulas, test parameters, isolation methods or results. Both verify modes hash every file in the complete parent project and byte-compare their [generic](evidence/exact-root-certificate.json) or [held-out](evidence/heldout-certificate.json) archived certificate. They emit JSON and do not write files. `compute` emits the recomputed certificate as JSON without saving it.

To deliberately propose fresh copies of all four project certificates:

```sh
python check.py regenerate --output-dir /absolute/outside/path/new-waning-proposal
```

The absolute output must be outside the complete parent project and must not exist. Its parent must exist. The proposal includes all four recomputed JSON files and a proposed parent manifest; no archive is replaced. A changed mathematical input requires mathematical review, not merely new hashes. Original material is covered by the parent's directory-scoped [MIT license](../LICENSE); cited sources retain their rights.
