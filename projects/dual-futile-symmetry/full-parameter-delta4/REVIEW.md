# Verification and separate review

A separate substantive AI mathematical review accepted the precise partial statement `all sixteen independent positive parameters and det(G)>=0 imply Delta4>0` on 26October2026. This is not external human peer review, formal proof-assistant verification, or a claim of absolute priority.

The reviewer independently built twelve stoichiometric reaction columns and all sixteen reaction derivatives from the primary network and ODEs. With species order `(A,B,C,E,F,U,V,X,Y)`, this gives the nine-dimensional Jacobian `J9=S*Jrate`. Conservation reduction uses `(Q1,Q2,U,V,Y,X)`, where `Q1=A+U` and `Q2=C+X`. Exact projection/embedding and tangent-intertwining identities recover every entry of `matrix.py`.

The independent polynomial route uses Faddeev-LeVerrier trace recursion with every division checked by multiplication back and a zero Cayley-Hamilton residual. A generic 4x4 Hurwitz determinant is derived separately with SymPy and evaluated over FLINT arbitrary-precision integer polynomials. This route does not import the principal-minor checker or a numerical solver.

Both implementations verify: `c6` has50terms and9negative coefficients; `Delta4` has1,712,913terms and266negative coefficients; `Q` has266positive coefficients; and `R=Delta4-Q*c6` has1,712,869positive coefficients with coefficient sum75,014,864. The strictness witness coefficient is one. They also verify `negative_part(Delta4)=a*c*e4*f10*v6*y12*Q`.

The review found no mathematical gap in this finite certificate and inference. Its source comparison distinguishes source Proposition4.2 (general positivity through Delta3) from Proposition5.6 (Delta4 positivity under mass action). The new lemma has an explicit determinant-sign condition, so it does not globally replace the latter proposition. The original stability conjecture remains open because Delta5 is unproved.

The public `independent_reconstruction.py` is adapted only for self-contained paths, public input hashes and read-only stdout receipts. No research prompts, internal job records, or manufactured gate certifications are distributed. The package's normal and optimized local runs, input-immutability check and corruption-rejection probes are recorded in `evidence/local-tests.json`.

The publication bundle additionally applies a shared complete-inventory guard at both checker entry points. The manifest includes all archived evidence; its own file and bytecode are excluded to avoid recursive hashing. Essential guards use explicit exceptions under optimization. Exact receipts omit elapsed times and the current manifest digest, so normal and optimized outputs can be compared directly without a receipt/manifest cycle.
