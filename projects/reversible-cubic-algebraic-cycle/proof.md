# A reversible cubic mass-action system with an algebraic limit cycle

## Statement and model

There is a reversible two-species mass-action reaction network in which every source and target complex has molecularity at most three, and whose positive-quadrant dynamics have the ellipse

\[
 h(x,y)=(x-2)^2+4(y-5)^2-1=0
\]

as a hyperbolic attracting limit cycle. Reversible here describes the reaction graph: each used directed reaction has its reverse. It does not assert detailed balance of the rate constants or an experimentally validated biochemical mechanism.

The complete reproducible witness is `witness.json`. Let
\(V=\{(i,j)\in\mathbb Z_{\ge0}^2:i+j\le3\}\).
For each distinct ordered pair \(v,w\in V\), include the reaction
\(v_1X+v_2Y\longrightarrow w_1X+w_2Y\), with rate one except for the ten rational replacements in `exception_rates`. All replacements are positive. Thus all 90 directed edges are present, all ten complexes are mutually reachable, and both species participate. In particular the network is weakly reversible and belongs to the source class \(W_3\), with no ghost complex or higher-molecularity target.

Its vector field is exactly
\[
 F(x,y)=\sum_{v\ne w} k_{vw}x^{v_1}y^{v_2}(w-v).
\]
All rates are fixed finite positive rational numbers; there is no limiting parameter. Coefficients and time are nondimensional as in the source definition.

## An invariant positive oval without equilibria

Direct rational polynomial expansion, reproduced by `check.py`, gives
\[
 F=hL+q(-h_y,h_x),\qquad \nabla h\cdot F=Kh,\qquad K=\nabla h\cdot L,
\]
where the two affine components of \(L\), the quadratic \(q\), and the quadratic \(K\) are explicitly stored in the witness. These are polynomial identities, not identities only on the oval. The checker reconstructs \(F\) from the reaction rates before checking them.

Parametrize the ellipse by \(x=2+\cos\theta\), \(y=5+\tfrac12\sin\theta\). It lies in \([1,3]\times[9/2,11/2]\), has nonvanishing gradient, and is therefore a smooth positive compact oval. Using
\[
 C(t)=\frac{1-t^2}{1+t^2},\qquad S(t)=\frac{2t}{1+t^2},
\]
the exact value of \(q(2+C(t),5+S(t)/2)\) is
\[
 \frac{19465008821284t^4+108647529059908t^3
 +174707253133865t^2+108020928482124t+155417707052014}
 {454003455620(1+t^2)^2}.
\]
Exact Sturm counting shows that its numerator has no real roots. Its value at zero and leading coefficient are positive; hence the numerator is positive on the entire real line. The missing parametrization point is \((1,5)\), where
\(q=374327092717/8730835685>0\).
Consequently \(q>0\) everywhere on the ellipse.

On this oval the vector field is \(q(-h_y,h_x)\), so it is never zero. In the stated angular coordinates it satisfies \(\dot\theta=4q>0\). Compactness bounds this speed away from zero and infinity. The invariant oval is therefore one periodic orbit with a finite positive period; periodicity has not been inferred from a plotted trajectory.

## Isolation and attraction

In a tubular neighborhood of this regular oval, \(h\) is a valid transverse coordinate. Since \(\dot h=Kh\), the nontrivial Poincare multiplier is
\[
 \mu=\exp I,\qquad
 I=\int_0^{2\pi}\frac{K(2+\cos\theta,5+\sin\theta/2)}
 {4q(2+\cos\theta,5+\sin\theta/2)}\,d\theta.
\]
The formula follows by linearizing the scalar transverse equation along the periodic solution. Variation of the section return time adds no first-order term because \(\dot h=0\) on the oval.

For \(0\le t\le1\), use the four rotations
\((u_j,w_j)=(C,S),(-S,C),(-C,-S),(S,-C)\).
Together they parametrize the four successive quarter arcs once with positive orientation. As \(d\theta/dt=2/(1+t^2)\),
\[
 I=\int_0^1G(t)\,dt,\qquad
 G(t)=\sum_{j=0}^{3}
 \frac{K(2+u_j,5+w_j/2)}
 {2q(2+u_j,5+w_j/2)(1+t^2)}.
\]
The witness stores the numerator and denominator of this exact rational function, of degrees 16 and 18. The checker independently rebuilds them from \(K\) and \(q\). Exact Sturm counts give zero numerator roots and zero denominator roots on \([0,1]\), while
\[
 G(0)=-\frac{264624322128746825582180724627900821454239466649810046218}
 {28462632417020934682407571540030115554849232562755262615}<0.
\]
Continuity therefore gives \(G(t)<0\) throughout the interval and \(I<0\). It follows that \(0<\mu<1\). The periodic orbit is hyperbolic attracting and isolated, hence a limit cycle. This proves the statement.

## Relationship to the source and limitations

Craciun and Erban, *Planar chemical reaction systems with algebraic and non-algebraic limit cycles*, Journal of Mathematical Biology 90, 64 (2025), arXiv:2406.05057v2, Definitions 2-3, Lemma 10 and the paragraph preceding Lemma 11, give a cubic merely kinetic algebraic-cycle example and a quartic weakly reversible example, while leaving existence in \(W_3\) open. The explicit construction above addresses that existence question. It makes no claim to classify all cubic networks, minimize the number of reactions, or determine a maximum cycle count. The source also proves \(W_2\) has no limit cycles, so the degree-three example, if the stated witness is accepted, is sharp with respect to the molecularity degree.

Classical ingredients are the invariant-conic representation, exact weakly reversible realization checks, Sturm's theorem, and the transverse multiplier formula. The proposed contribution is this explicit admissible cubic reversible witness, not invention of those tools. Literature priority and human expert verification remain unconfirmed. The work used substantial AI assistance in construction, coding and separate checking; such checks are not human peer review. It concerns a mathematical mass-action model, not experimental biology or a therapeutic claim.
