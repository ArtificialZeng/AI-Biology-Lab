# A positive, regular CRN stress case for species-coordinate shape

This original proof/code packet gives a six-reaction mass-action family whose steady-state ideal has four regular complex zeros and a unique globally stable positive equilibrium, while **no reordering of the two species variables gives the proposed shape basis**. The conclusion holds for every strictly positive choice of its six independent rate constants.

The diagnostic concerns the species-coordinate shape program in Ferrari et al., [arXiv:2501.12233v2](https://arxiv.org/html/2501.12233v2), Sections 1.2 and 4. It distinguishes reordering existing species coordinates from taking a generic linear combination. The latter repair is already covered by classical shape-lemma theory and is calculated explicitly here for this family.

## Family and result

```text
0 -> X          rate a
0 -> Y          rate b
2X -> X         rate c
2Y -> Y         rate d
2X -> X+Y       rate p
2Y -> X+Y       rate q
```

All a,b,c,d,p,q>0. With A=c+p, B=d+q and Delta=AB-pq=cd+cq+pd, define

```text
u=(Ba+qb)/Delta >0,   v=(pa+Ab)/Delta >0.
I=(x^2-u, y^2-v).
```

The reaction stoichiometry has no nonzero conservation law and no nontrivial independent decomposition. The ideal is radical and zero-dimensional; all four sign choices `(±sqrt(u),±sqrt(v))` are simple torus zeros. Each species coordinate identifies two distinct points, so neither native coordinate can parameterize the other polynomially across the full complex variety. Saturating away the coordinate axes preserves the ideal.

The physically positive equilibrium `(sqrt(u),sqrt(v))` is unique and globally asymptotically stable on the nonnegative quadrant. The proof uses an explicit coercive Lyapunov function and a symmetric negative definite weighted reaction matrix. At unit rates the ideal is `(x^2-1,y^2-1)` and the positive-state eigenvalues are exactly -2,-6.

For `t=x+lambda*y`, the exceptional choices are exactly `lambda=0` or `lambda^2*v=u`. Outside them, the proof displays the quartic eliminant and polynomial formulas for both species. This is a family-specific illustration of the established mixed-coordinate repair.

## Scope and sources

The original broad conjecture already has a credited public community counterclaim based on a nonradical origin-only example. The paper's own Example 1D already demonstrates abstract radical ideals with coordinate-projection collisions. First-disproof priority and a new general algebraic obstruction are not claimed. The proposed increment is the explicit admissible, dynamically well-behaved, all-positive-rate CRN realization and its reproducible certificates.

The broader task of finding practical necessary-and-sufficient native-shape conditions remains open here. The source Theorem 2 allows generic linear changes of variables; this packet is consistent with it. The source's reported GRN/Wnt calculations and the use of Gröbner bases to find positive equilibria are not challenged by this example. The example is an abstract mass-action model, with no claimed measured biochemical pathway or new experimental effect.

The current author preprint and repository were inspected. The publisher preview verifies the journal citation [10.1007/s40574-025-00515-1](https://doi.org/10.1007/s40574-025-00515-1); the final published full-text wording could not be compared. All precise claim attribution here targets arXiv v2. See [source.md](source.md) and [problem.md](problem.md).

## Reproduce

From this packet's directory:

```sh
python -m pip install -r requirements.txt
python -B verify.py verify
```

The wrapper checks the complete project manifest and replays archived exact certificates without changing project files. Owner checks reconstruct the ODE, enumerate all 31 unordered nontrivial reaction partitions, verify symbolic ideal/derivative/saturation/repair identities, and run 66 supplementary rational cases. These finite cases do not replace the all-rate analytic proof.

Normal and optimized wrappers deliberately execute scientific children with assertions enabled (`PYTHONOPTIMIZE=0`). Direct optimized execution of an assert-based child is not a verification guarantee. To regenerate a certificate, redirect a scientific checker's output outside this project first, compare its full content and obtain renewed review for mathematical changes.

## Review status

The exact public proof received a separate AI mathematical check. Its contribution was judged narrowly sufficient for a short source-specific code/proof note, with no basis established for a substantial standalone algebra or reaction-network theory article. [The receipt](provenance/review-receipt.json) binds the proof and both checkers. There is no author confirmation, human peer review, priority certification or proof-assistant formalization. Genuine trusted-runner job IDs are unavailable; its structural freeze failed, so no manuscript/PDF or formal runner release is included.

The independent checker was first run with SymPy 1.14.0 and replayed with the pinned 1.13.3; all mathematical fields agree and the dependency-version field records the pinned run. Its certificate's `math.md` label refers to the original reviewer-side analytic report; the public analytical argument is in [proof.md](proof.md). The checker has no runtime dependency on that report or any private path.

AI assistance was used for derivation, code and separate checking passes. This is an independent computational project without claimed institutional endorsement.
