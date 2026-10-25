# Candidate computer-assisted lemma and analytic argument

Status: base-cycle scope accepted by a fresh AI-assisted mathematical audit after Linux x86/FILIB reproduction. [Linux evidence](evidence/linux-interval-output.txt), [independent exact audit](independent-audit/independent_check.json), and the [finite-K persistence argument](finite_kinase_persistence.md) are available. This is not external human peer review or priority confirmation.

## Frozen ODE and intended conclusion

Let u,v,w,a,b be the limiting concentrations, p=8-a-b, s=40-u-v-w-a-b. The polynomial ODE is

u'=100b-u, v'=u-v, w'=v-100wp+a,
a'=100wp-(19/10)a, b'=3sp-101b.

The intended conclusion is existence of a locally unique periodic orbit crossing u=13/4 positively in the displayed initial box, and orbital asymptotic stability. No global attraction, uniqueness among all cycles, minimal-period assertion, exact K-independence, or finite-K threshold is claimed.

## Algebraic physical-domain argument

Adjoin the dependent variables p,s. Differentiating their definitions gives

p'=-(100w+3s)p+(19/10)a+101b,
s'=(9/10)a+b-3ps.

Every coordinate of (u,v,w,a,b,p,s) has nonnegative derivative at its zero boundary when all coordinates are nonnegative. The affine constraints a+b+p=8 and u+v+w+a+b+s=40 are conserved. Hence the nonnegative physical polytope is forward invariant and compact. Strictly positive initial coordinates stay strictly positive for every finite time: each equation has the form gain-loss*coordinate, with nonnegative gain and bounded nonnegative loss on this polytope, so Gronwall's inequality supplies a positive exponential lower bound. A validated fixed point whose initial box is strictly inside this domain therefore yields a strictly positive whole periodic orbit without a sampled-minimum argument.

## Existence check

Use coordinates z=(v,w,a,b) on the transverse section u=13/4 and P the next positive-crossing return map. CAPD IPoincareMap encloses P and its derivative on the entire initial box X. The source code fixes the exact rational coefficient 19/10 by enclosing interval division, rather than treating the nearest binary 1.9 as exact.

Put F=P-Id. c is the explicit binary-double center; R is an explicit binary-double numerical inverse preconditioner. Its interval determinant excludes zero. The checker encloses

K(c,X)=c-R F(c)+(I-R DF(X))(X-c).

Strict K(c,X) inclusion in the interior of X gives existence. The additionally checked infinity norm of I-R DF(X) below one gives injectivity/uniqueness in the box. These are finite interval inequalities, not a claim based on a small floating residual.

The return map's output lies in a time interval bounded away from zero, and its enclosed normal velocity 100b-u is strictly positive at return. The phase direction is removed by the four-dimensional section.

## Attraction check

D encloses DP(X). Interval multiplication encloses the fourth power of every one of these matrices, including the derivative DP(z*) at the certified fixed point. If its infinity norm bound is <1, then rho(DP(z*))^4 <= ||DP(z*)^4||_inf <1. Thus all four transverse multipliers have modulus <1, giving orbital asymptotic stability. Numerically tiny multipliers are never declared exact zeros.

## Separately audited finite-K consequence

The exact six-dimensional singular reduction in the repository has epsilon=1/K and epsilon*y'=-y+epsilon*(u+v)*y+u+100b. At epsilon=0 its critical graph is y=u+100b with fast normal derivative -1. The strict-interior hyperbolic limiting orbit and the local compact-neighborhood hypotheses have been checked. The separate argument uses local Fenichel persistence and hyperbolic-cycle persistence to obtain an attracting orbit within each fixed six-dimensional compatibility class for every sufficiently large finite K, without a numerical K0. The existing instability threshold 4660166/209569 is not such a K0. The full argument is in finite_kinase_persistence.md; the complete original conjecture remains unresolved.

## Literature and method credit

The ODE is derived from Conradi-Mincheva-Shiu (2019), https://arxiv.org/html/1809.02886v2 . Its oscillations were already numerically displayed there. Direct related Grasse 2023/Grasse-Conradi 2025 conference records prevent an unqualified priority claim. The checker uses standard interval Poincare/Krawczyk methods; CAPD's published examples document this technique: https://capd.sourceforge.net/capdDynSys/docs/html/a05240.html . CAPD is GPL-3.0; our original source is MIT and no dependency sources/binaries should be republished without retaining their license terms.
