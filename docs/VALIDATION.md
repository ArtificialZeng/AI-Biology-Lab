# External validation

No external human review is reported for the first project. We distinguish four questions:

1. **Model fidelity:** do the reaction rates, conservation laws, and reduced equations match the cited model?
2. **Mathematical correctness:** do the rational equilibrium proof, interval cycle/sensitivity computations, and local persistence argument establish their stated quantified scopes?
3. **Novelty:** is an equivalent proof or stronger result already available, including unpublished follow-up work?
4. **Significance:** is this quantitative partial result useful enough for a specialist publication?

Reviewers can begin with the [current project statement](../projects/mixed-phosphorylation/README.md), [new interval proof and period correction](../projects/mixed-phosphorylation/validated-cycle/README.md), and the [earlier equilibrium-only note](../projects/mixed-phosphorylation/manuscript/main.pdf). `python scripts/reproduce.py verify` reruns the rational equilibrium checks; the interval integration has a separate C++/CAPD route.

The original authors Carsten Conradi, Maya Mincheva, and Anne Shiu are natural contacts for the status of Conjecture 6.2. Kim-João Grasse presented directly related robustness work at SBHD 2023 and with Conradi at BIOMATH 2025. Obtaining the full follow-up results is an unresolved priority task. The conference programme and abstract do not themselves settle whether the present theorem is new.

Useful specialist checks include the physical parameter interval, coverage of the entire unbounded K ray, signs of denominators, all Bernstein identities, the restriction to the six-dimensional compatibility class, and the Routh root count. Review the written argument as well as the executable checks.

Any verification issue should identify the claim version or commit, verifier command, outcome, and limitations. A GitHub Actions pass confirms that the encoded checks ran successfully; it does not certify a research contribution. External feedback will be recorded publicly only with permission and with its scope preserved.

For the new result, check the phase-fixed Poincare map, strict Krawczyk inclusion, transverse derivative bound, exact first slow-manifold jet, interval sensitivity solve and time derivative. Verify the finite-r parameter Fenichel/IFT hypotheses and fixed-compatibility-class scope. The candidate increment is the model-specific nonzero coefficient enclosure; generic sensitivity methods and qualitative robustness have prior work. No explicit K0 or remainder constant is claimed.
