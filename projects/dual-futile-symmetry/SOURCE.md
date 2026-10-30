# Source and frozen scope — 26 October 2026

Nicola Vassena, *Sequential and distributive dual futile cycle: Hopf bifurcation can occur under parameter-rich kinetics but cannot occur under mass action kinetics*, [arXiv:2608.27081v2](https://arxiv.org/html/2608.27081v2), revised 14 September 2026.

- Section 2 gives the reduced Jacobian in coordinates `(Q1,Q2,U,V,Y,X)`, with sixteen strictly positive partial-derivative parameters. Its graph involution alone does not impose equivariant kinetics.
- Section 3, Conjecture 3.1 asks whether `det(G_red)>0` implies that `G_red` is Hurwitz for **all sixteen** positive parameters. The checked version explicitly leaves this unresolved.
- Section 5 proves the separate mass-action no-Hopf theorem. That theorem is not an open target here.

Our subsidiary question imposes `a=c`, `b4=b10=b`, `e1=f7=e`, `e4=f10=f`, `u2=x8=u`, `u3=x9=w`, `v5=y11=v`, `v6=y12=z`. The eight variables `a,b,e,f,u,w,v,z` are strictly positive. These are Jacobian-parameter equalities; no assertion is made that arbitrary microscopic mass-action constants realize all such tuples.

The deliverable is an exact classification on that lower-dimensional slice, not a proof of the source's unrestricted conjecture. No statement about a full asymmetric neighborhood follows merely from this classification.

The author's [MATLAB repository](https://github.com/nvassena/dualfutile_hopfnohopf) was inspected at `02b3c865f0618ce0e104e5b8a7e2111848580802`. No directly equivalent threshold was found in the bounded primary-source search on the date above. Priority and usefulness remain questions for the original author; absence from a search is not proof of novelty.
