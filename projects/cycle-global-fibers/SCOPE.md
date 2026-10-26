# Frozen mathematical scope

Let n>=3, 2<=p<=n-1. Take the directed cycle 1->2->...->n->1, input {1}, output {p}, and exactly two leak vertices r,t with 1<=r<p<=t<=n. Write k_i for the edge leaving vertex i, e_r,e_t for leaks, and d_i=k_i+e_i at leaks and d_i=k_i elsewhere. All parameter labels are fixed; relabeling different vertices counts as a different parameter vector.

Ahmed et al. Theorem 3.6 already implies generic local identifiability. Conversely, all generically locally identifiable two-leak cycles with exactly one input and one output reduce by rotation to this family: p=1 cannot interlace; p=n would need a leak at n and is exceptional; otherwise interlacing requires one leak in [1,p-1] and one in [p,n]. We cite this local classification instead of claiming it.

The coefficient data are the monic denominator D(s)=prod_i(s+d_i)-C and numerator N(s)=kappa*T(s), with C=prod_i k_i, kappa=prod_{i<p}k_i and T(s)=prod_{i>p}(s+d_i).

Target theorem: On the data domain kappa*C!=0, U(s)=(D(s)+C)/T(s) and T(s) each have nonzero simple roots. Put A=U(0), B=T(0). Then the full complex coefficient fiber has exactly p!(n-p)! distinct points, explicitly parametrized by independent assignments of the roots -alpha of U to vertices 1..p and roots -beta of T to vertices p+1..n. At the leak vertices set

k_r=kappa/(prod_{i<p,i!=r}d_i),
k_t=C/(kappa*prod_{i>=p,i!=t}d_i),
e_r=d_r-k_r, e_t=d_t-k_t;

set k_i=d_i elsewhere. C is recovered uniquely as minus the constant remainder of D on division by T. For data from strictly positive parameters all alpha,beta,kappa,C are positive. With a=C/(kappa*B), b=A/kappa and m=#{alpha:a<alpha<b}, the number of strictly positive parameter preimages is exactly m*(p-1)!*(n-p)!. Every value m=1,...,p occurs on a nonempty open positive-parameter region. The count is allowed to change across regions, so we do not call the positive count a single generic algebraic degree.

Nearest prior result: exact local classification and the coefficient map (Ahmed et al.), plus zero/one-leak permutation-based degree analyses (Huneke–Keen). Proposed delta: exhaustive explicit inverse branches for all locally identifiable single-input/single-output two-leak placements, with exact positivity inequalities and all attainable positive counts. Verification route: a complete algebraic bijection, independent determinant/adjugate checks, and exhaustive rational enumeration of bounded fixtures. No mathematical acceptance or novelty certification is asserted by the researcher.
