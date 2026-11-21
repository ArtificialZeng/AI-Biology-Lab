# Primary sources and bounded novelty log

Literature verification date: 2026-11-21 (UTC+8). Searches are bounded to the listed APIs and repositories and do not constitute universal priority certification.

## Target conjecture source

* **David A. Brewster, Jakub Svoboda, Dylan Roscow, Krishnendu Chatterjee, Josef Tkadlec, and Martin A. Nowak**, *Maintaining diversity in structured populations*, **PNAS Nexus** 4(8):pgaf252 (published 2025-07-31; `arXiv:2503.09841v2`, updated 2025-09-16).
  * DOI: `10.1093/pnasnexus/pgaf252` — https://doi.org/10.1093/pnasnexus/pgaf252
  * arXiv: https://arxiv.org/abs/2503.09841 (v2 source inspected: `main.tex`, `si/proofs/star.tex`, `si/proofs/well_mixed.tex`)
  * Europe PMC / PMC: `PMC12363668` / `PMID:40838023` (full text verified)
  * OpenAlex ID: `W4413186912`
  * **Verified definitions**:
    * Star graph $S_n$ with $N = n + 1$ vertices ($n$ leaves of degree $1$, $1$ central vertex of degree $n$, simple connected undirected graph without self-loops, `main.tex` lines 536–538, SI Section 1).
    * Initial state: all $N$ vertices start with pairwise distinct types ($|\tau| = |V|$, $X_0(V) = \tau$, SI Section 1).
    * Dynamics: neutral discrete-time birth-death (`bd`) updating with one birth-death replacement per time step (`main.tex` Eq. (1), SI Section 1, Methods Eq. (15)).
    * Conjecture Eq. (9) (`main.tex` lines 556–560): $T_N = n^3 - n^2 + n \cdot H_n$, where $H_n = \sum_{k=1}^n 1/k$. Only $\Theta(N^3)$ (`bd`) and $\Theta(N \log N)$ (`db`) asymptotic bounds were proved in SI Section 5 (`si/proofs/star.tex`).

## OpenAlex citing works of `W4413186912` (all 6 retrieved on 2026-11-21)

1. `W7154604876` (2026-04-17) — Jakub Svoboda, Hossein Nemati, Josef Tkadlec, Kamran Kaveh, Krishnendu Chatterjee, *The effect of the fitness gradient on fixation probability*, *Nature Communications* 17(1):5325, DOI `10.1038/s41467-026-71777-2` (also `PMID:41997932`). Studies 2-type fixation probability under spatial fitness gradients; does not prove Eq. (9).
2. `W7163063073` (2026-06-01) — *A voter-model representation of the multi-allelic Moran process: Exact stationary distributions and diversity thresholds from well-mixed populations to complex networks*, *Journal of Theoretical Biology*, DOI `10.1016/j.jtbi.2026.112517` (`PMID:42229807`). Studies stationary allele distributions under parent-independent mutation; does not prove Eq. (9).
3. `W7115587857` (2026-09-30) — David A. Brewster, Yichen Huang, Michael Mitzenmacher, Martin A. Nowak, *Mixed updating in structured populations*, *PLOS Computational Biology* 22(9):e1014829, DOI `10.1371/journal.pcbi.1014829` (`PMC13645234`, `PMID:42814787`, related `arXiv:2511.18252v3`). Proves 2-type fixation probability formulas on stars/cycles and general fixation-time bounds under mixed Bd/dB updating; full text inspected (`PMC13645234.txt` and `2511.18252.txt`), does not address $N$-type consensus time Eq. (9).
4. `W4414010756` (2025-09-05) — Nikhil Sharma, Suman G. Das, Arne Traulsen, Joachim Krug, *Population structure amplifies mutation load*, bioRxiv `10.1101/2025.09.04.674237`. Studies mutation load under deleterious mutations; does not prove Eq. (9).
5. `W7162518942` (2026-05-26) — David A. Brewster, Gabor Lippner, Josef Tkadlec, Martin A. Nowak, *Fixation location in structured populations*, `arXiv:2605.26411v3`. Studies conditional fixation location of a single mutant; full text inspected (`2605.26411.txt`), does not address Eq. (9).
6. `W7202117459` (2026-08-06) — Michal Pecho, Jakub Svoboda, Lenka Kopfová, Josef Tkadlec, Krishnendu Chatterjee, *Replacers and their evolutionary stability in the Moran process on graphs*, `arXiv:2608.06552v1`. Studies neighbor-aware "replacer" phenotypes; full text inspected (`2608.06552.txt`), does not address Eq. (9).

## Additional API checks

* **Crossref REST API** (`https://api.crossref.org/works/10.1093/pnasnexus/pgaf252`): verified title, all 6 authors, journal (*PNAS Nexus*, Vol. 4, Issue 8, published 2025-07-31).
* **Europe PMC API**: `get_citations MED 40838023` returned 2 indexed citations (`PMID:42814787` and `PMID:41997932`), both already in the OpenAlex citing set above.
* **arXiv API** (`search_arxiv.py`): targeted searches for `(ti:"consensus time" OR ti:"absorption time" OR ti:"diversity time" OR abs:"star graph") AND (abs:"Moran" OR abs:"Brewster")` and author searches for `Brewster`, `Tkadlec`, `Svoboda` through 2026-11 returned no paper proving Eq. (9).
