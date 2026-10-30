# Primary source and exact reaction mapping

E. Franco and J.J.L. Velázquez, *The adaptation property in non-equilibrium
chemical systems*, Journal of Mathematical Biology92:32(2026), published13
February2026. [Version of record](https://doi.org/10.1007/s00285-026-02344-y).

Published Example3.8, PDF page32, displays source species1 <-> species2+species3
and species3 <-> species1. Set S=1, X=3, Y=2, and use positive forward/reverse
constants a,b for S<->X and c,d for S<->X+Y. Chemostatting S=u(t) gives exactly
the model in model.json and proof.md. The displayed network has conservation
vector(1,0,1); its adapted output is species2. The prose's species3 claim is an
indexing inconsistency, not the new mathematical contribution.

The underlying reversible network satisfies detailed balance for every positive
choice of these four constants. Its conservative completion S+Z<->X+Y,
S<->X, with Z fixed at1, is consistent with the source's completion mechanism.
No new thermodynamic realization or constant-input adaptation is claimed.

Section5 of the published article raises questions about periodic response and
mean behavior. It does not assert a universal mean-preservation conjecture.
The present theorem treats this one exact motif. The older
[arXiv2502.17375v2](https://arxiv.org/abs/2502.17375v2) does not contain the
published Example3.8 or published Section5; use the journal version for these
locators. Source PDF/text and indexing were checked against the published PDF.

No paper-specific author code repository was established in the bounded source
search. This additive independent code/proof note is suitable for the lab
repository; no unrelated upstream project is used as a substitute.
