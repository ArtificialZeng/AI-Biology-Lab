# Rooted two-sink images and a boundary obstruction after suppression

This is a proof candidate for independent checking, not external human validation.
The statements below use a uniform stationary root and the displayed-tree model.

## Domain and the two distinct source questions

Identify A,C,G,T with 0,1,2,3 in the Klein four group, with addition XOR.
Use the character table

```
H = [[1, 1, 1, 1],
     [1, 1,-1,-1],
     [1,-1, 1,-1],
     [1,-1,-1, 1]].
```

A group-based edge kernel has entries M(x,y)=p(x+y), where p is a nonnegative
probability vector. Write its Fourier eigenvalues as (1,a_C,a_G,a_T). Throughout,
for the rooted classification, all three nontrivial eigenvalues are strictly between zero and one. K2P imposes
a_C=a_T; JC imposes equality of all three. K3P imposes no extra equality. Every
inheritance probability and its complement are strictly positive. Identity edges
are excluded. Zero individual entries are allowed in this domain. Let D denote
the selected edge family and D+ its strictly positive-entry subset.

This is the explicit spectral reading of Brits et al., arXiv:2607.12919v3 §2.3,
Eq. (3)/(4). Its nonidentity/positive-definite wording alone is weaker than its
subsequent strict spectral assertion; we use the latter and do not identify it
with continuous-time embeddability.

Sullivant's rooted reduction retains two hidden attachment vertices b1,b2 and
replaces their ancestral subgraph by a new root and two edges to those vertices.
Brits' semi-directed reduction instead contracts the sub-blob and suppresses its
degree-two image to a single edge. These are different parameterized targets on
the boundary of D. The zero-entry example below does not refute Sullivant's
splittable-model hypotheses or his general question for arbitrary root laws.

## 1. Exact joint-image theorem for the rooted two-sink setting

Let B be any finite rooted binary DAG with exactly two distinguished sinks b1,b2,
all its other vertices hidden, and every vertex ancestral to a sink. The sinks
are distinct and neither is the root. For each reticulation independently keep
one parent, with any fixed inheritance probabilities in (0,1). Give every edge
a kernel from D and the root the uniform law. The possible joint laws of the
two sinks are exactly

    P(X_b1=x, X_b2=y) = (1/4) M(x,y),  M in D+.

The same statement permits arbitrary finite parent indegrees with all parent
weights positive; the binary case suffices for the cited networks. No blob or
level bound is required for this algebraic theorem. In a larger network, the
result applies to an ancestral hidden subgraph with these two exit sinks and
unchanged conditional descendants; no edge may enter this ancestral subgraph
from outside.

Proof of the forward direction. Positive spectrum gives p(0)>0. If a kernel
had support of size at most two, that support would lie in {0,h} for some h.
A nontrivial character trivial on h would then have eigenvalue one, contrary
to the domain. Thus every permitted kernel has at least three support elements.
Two subsets of a four-element group, each of size at least three, meet after
every translation. Consequently the convolution of any two permitted kernels
is strictly positive at every group element. Further convolution with stochastic
kernels preserves strict positivity.

For each displayed parent choice sigma, the unique path between the two sinks
has length ell_sigma>=2: the sinks are distinct and neither can be a parent of
the other. Its difference kernel is the convolution of its edge kernels. This
follows by marginalizing the uniform-root group-based tree, or by independent
edge increments. Its three nontrivial Fourier eigenvalues are the products of
the corresponding edge eigenvalues. Mixing the finitely many parent choices
therefore gives a strictly positive kernel, with each eigenvalue in (0,1), and
preserves the JC/K2P equality constraints. This proves containment in D+.

Proof of surjectivity, with an explicit construction. Fix M in D+, with
probability vector p and Fourier eigenvalues a_g. Fix the desired strictly
positive parent-choice weights w_sigma; they need not be fitted to M. Let
L=max ell_sigma>=2 and

    F=max(max_{g != A} a_g, 1-4 min_x p(x)),
    t=1-(1-F)/(2L),
    f=sum_sigma w_sigma * t**(ell_sigma-1).

Then 0<F<1 and 0<t<1. Bernoulli's inequality gives
t**(ell_sigma-1) >= 1-(ell_sigma-1)(1-t)>F, so F<f<1.
Let J_t be the JC kernel with nontrivial eigenvalues t, and define

    T = M J_f^{-1},
    p_T(x) = [p(x)-(1-f)/4]/f,
    a_g(T) = a_g(M)/f.

The choice f>1-4 min p makes every entry of T positive. Entries sum to one;
f>max a_g ensures every nontrivial eigenvalue of T is in (0,1). The selected
model's equality constraints are preserved, so T is in D+.

Assign T to every edge entering b1, and J_t to every other edge of B. Each
displayed b1-b2 path contains exactly one edge entering b1, so its kernel is
T J_{t**(ell_sigma-1)}. Since JC kernels satisfy the same affine mixture rule
as their single nontrivial eigenvalue, their mixture is exactly T J_f=M.
This realizes the prescribed joint law, using no identity edge, zero inheritance
probability, approximation, or limit in the final parameters. All formulas are
rational when the target and inheritance probabilities are rational. QED.

Corollary (rooted fork replacement). The two-edge fork from a new uniform root
to b1,b2 has exactly the same joint image D+: its effective kernel is PQ, which
lies in D+, and the preceding construction with path length two realizes every
member of D+. Hence replacing any such hidden ancestral two-sink subgraph by
that fork preserves the full descendant distribution image. This answers the
rooted equality question in Sullivant §5 for the explicitly specified uniform
stationary JC/K2P/K3P models. It is not a statement about arbitrary root laws,
all transition-matrix families, or every semi-directed 2-sub-blob.

## 2. Counterexample for the single-edge suppressed model

Take the rooted binary network

```
rho -> u,v; u -> a,b; v -> a,b;
a -> s; b -> t; s -> 1,2; t -> 3,4.
```

Suppressing rho produces the semi-directed diamond with undirected edge u-v
and reticulation edges u->a,v->a,u->b,v->b. The induced B={u,v,a,b} is a
root-trapping maximal 2-blob: its five edges lie in cycles, its only outside
attachments are a-s,b-t, and both attachments must be outgoing because a,b
already have two fixed incoming reticulation edges. A root outside B cannot
enter either attachment. The original root is the last stable ancestor.
Contracting B and suppressing its degree-two image yields the quartet 12|34.
See graph-proof.md for the detailed graph argument; its finite enumeration is
a supplementary check, not a claim to enumerate arbitrary hidden expansions.

On the collapsed quartet give pendant edge 1 the JC eigenvalue 1/4 and edges
2,3,4 the JC eigenvalue 1/2. Give its internal edge the K2P kernel

```
M = [[3/5,1/5,0,1/5],
     [1/5,3/5,1/5,0],
     [0,1/5,3/5,1/5],
     [1/5,0,1/5,3/5]].
```

Its eigenvalues in the stated H convention are (1,3/5,1/5,3/5), so it satisfies
the exact domain and K2P equality, and also the containing K3P model. This target
has an all-rational binary rooting: insert a root on pendant edge 1 and give
both resulting edges JC eigenvalue 1/2. Keep the other four rooted edges as
specified. Thus target membership does not assume a prohibited identity edge
or unsplittable boundary root edge.

The nonmembership argument below in fact permits the broader stochastic
positive-definite/nonidentity domain in which some nontrivial eigenvalues may
equal one. This avoids relying on the erroneous implication that nonidentity
and positive definiteness alone force every nontrivial eigenvalue below one.
It still uses a uniform root and excludes identity edge matrices.

Here is the additional support lemma. For any nonidentity positive-definite
Klein-four stochastic kernel, p(0)>0 and its support contains 0 and at least
one other element. If S,T are two such supports, then S+T is either the entire
group or a two-element subgroup. Indeed, a support of size at least three plus
one of size at least two fills the group; two two-element supports {0,h} and
{0,k} fill the group when h!=k and sum to {0,h} when h=k. Convolving their
product with any further probability kernel produces full support or a union
of cosets of that two-element subgroup. Its support therefore has cardinality
two or four, never three. Nonnegative convolution has no cancellation, so this
is a statement about exact support, not a lower bound or approximation.

For any parameters on the original diamond network, its four displayed trees
all have the quartet split 12|34. The internal s-t kernel C is the product of
the two common attachment kernels and the mixture of the a-b path kernels.
The two attachment edges are nonidentity positive-definite kernels, while the
core mixture is stochastic. By the additional support lemma, C cannot have
exactly three support elements, regardless of all other parameter choices.
All its Fourier eigenvalues c_g remain positive, including under the broader
domain. If d_i,g are the four pendant eigenvalues, the observable
pair correlations R_ij,g=E[chi_g(X_i+X_j)] satisfy

    R_12,g=d_1,g*d_2,g, R_34,g=d_3,g*d_4,g,
    R_13,g=d_1,g*c_g*d_3,g,
    R_24,g=d_2,g*c_g*d_4,g.

Therefore c_g**2=(R_13,g*R_24,g)/(R_12,g*R_34,g). Denominators and c_g are
positive, so the unique positive root determines each c_g even if every
pendant and attachment parameter is retuned. For our target this reconstructs
c_C=c_T=3/5 and c_G=1/5, which forces C=M by invertibility of H. But C cannot
have support of size three, whereas M does. This contradiction proves exact
nonmembership in the original network image.

Consequently the universal equality asked after Brits v3 Lemma 5.1(i) is false
on its explicitly stated spectral/nonnegative-entry domain for K2P and K3P,
and also on the broader nonidentity positive-definite stochastic reading.
This does not refute the reported inclusion or closure equality; it gives a
point of strict inclusion. It does not settle the arbitrary-graph JC question.
The boundary target is not a strictly positive-entry or continuous-time K2P
target, so no such counterexample to those additional domains is asserted.
Brits §4 also describes its parameters as continuously generated; that wording
is stronger than the §2.3 stochastic/positive-definite conditions. This note
does not claim to refute a question with an added continuous-time hypothesis.
The rooted fork
corollary and this unrooted counterexample coexist because a two-edge fork
cannot produce this boundary kernel under the stipulated strict spectrum.

The closest source claims, assumptions, exact graph, rational leaf probabilities,
observable pair correlations and construction checks must be independently
audited before any public claim. No scientific impact on empirical phylogenies,
absolute historical priority, human peer review, or formal proof assistant
verification is asserted.
