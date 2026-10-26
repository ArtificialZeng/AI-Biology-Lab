# Source and precise discrepancy

Justin Eilertsen, Valery G. Romanovski, Santiago Schnell and Sebastian Walcher, [*Lumping of reaction networks: Generic and critical parameters*, arXiv:2606.28895v2](https://arxiv.org/abs/2606.28895v2), Example 8, printed pp.38–39. The arXiv revision is dated 2 July 2026; the title page is dated 7 July 2026. The inspected source was current v2 on 26 October 2026. No source PDF or page image is redistributed.

Its six-species model follows Desire T. Gijima and Enrique Peacock-López, [*A Dynamic Study of Biochemical Self-Replication*, Mathematics 8 (2020), 1042](https://doi.org/10.3390/math8061042). The source authors supply the reaction network, coordinates, conservation laws and local stability question.

At critical rates `k5=km5=0`, the printed p.38 matrix is

```text
[-k1*c1-k2*c2-km1, -k1*c1,             -k1*c1+km2;
 -k2*c2,           -k3*c2-k4*c1-km3,   -k3*c2+km4;
  k2*c2,            k4*c1,              -km2-km4].
```

Here `kmj` means `k_{-j}`. The preceding equation is

```text
y5' = k3*(-y1+c2-y5-y6)*(y1-y4-y5-y6)
      -k4*(-y1+c1-y4-y6)*y5-km3*y5+km4*y6.
```

Its derivative with respect to `y4` at zero is `-k3*c2`. Only the printed `(2,1)` entry differs from direct differentiation. Footnote 9 states: “Presumably this condition holds for all positive parameters, but the verification is not trivial.” The local text asserts stability in an open parameter region; the counterexample below concerns the universal extension of the literal printed matrix, not that open-region assertion.

The [canonical companion repository](https://github.com/santiago-schnell/lumping-of-reaction-networks) was inspected read-only. Its `gpl_replication_network` in `src/lumping_analysis/examples.py` correctly assigns `k3` to `X2+X3 -> X5`. The archived inspection listed 49 code/document files and separately checked `.gitignore` and `LICENSE`; none implements this Example 8 transverse Jacobian or its eigenvalue claim. The finding is a displayed-formula discrepancy, not evidence of incorrect production outputs. Bounded searches for the exact title or arXiv identifier with “Jacobian correction”, “Example 8” and “Hurwitz” found no overlapping primary correction on 26 October 2026; that does not establish priority or exclude an existing issue.

The source's printed pp.37 and 45 discuss finite-window approximation errors. Its qualified future observation about many systems over all positive times provides no universal hypotheses. The optional slow-drift example in [PROOF.md](PROOF.md) addresses only a separately specified stronger assertion.
