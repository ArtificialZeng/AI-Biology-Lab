# Verification and review limits

On 26 October 2026, a separate AI review accepted the restricted proof note and its suitability as an upstream verification contribution. It did not establish novelty or independent publishability. This is not external human peer review or author endorsement.

The main checker uses standard-library rational polynomial arithmetic and direct Leibniz determinants. It verifies `GP=P diag(A,B)`, the complete six-dimensional characteristic polynomial, both cubics, positive Routh expressions (86 and 56 positive monomials), and three threshold regimes. Three deliberately corrupted claims are rejected, including under optimized Python.

The second checker uses SymPy on the six-dimensional formulas. The additional `independent_reconstruction.py` starts instead from the twelve stoichiometric reaction vectors and formal rate gradients of the nine-species network. It constructs `J=N R`, derives the compatibility-class injection and projection, and obtains the six-dimensional matrix anew. Solving the commutator equations with the involution gives rank eight and exactly the eight stated parameter pairings. An independently chosen parity basis produces the same cubic factors. The zero boundary is simple; the negative regime has exactly one positive real root.

Thus the family is the full fixed locus of this derivative-parameter involution. General graph symmetry alone does not imply these kinetic equalities. No asymmetric neighborhood or arbitrary mass-action realization is covered.

The source already proves the mass-action no-Hopf theorem. The review also compared primary results of Hell–Rendall, Conradi–Feliu–Mincheva, and Feliu–Kaihnsa at the level of their kinetic assumptions. No directly identical eight-parameter classification was located in the bounded search. This does not prove priority, and the symmetry reduction is elementary enough that the source author may already know the corollary. The general sixteen-parameter conjecture remains unresolved.
