# Sharp peak cap for exclusive tri-ligand scaffold binding


## Statement

For every K_i>0 and alpha_i>0 (i=1,2,3), define the unique root U_i(X)>1 of

    U_i²+(K_i X-alpha_i-1)U_i-K_i X=0, X>0,

and

    y(X)=X(U_1-1)(U_2-1)/(U_1 U_2 U_3).

Then y has at most two strict local maxima on (0,infinity). The source's known exact two-peak baseline makes the cap sharp, but is not a new construction.

## 1. Source variables and a branch-safe stationary equation

Set q_i=K_i X, a_i=1+alpha_i, D_i=(q_i-a_i)²+4q_i,

    T_i=(1+(q_i-a_i)/sqrt(D_i))/2,
    mu_i=(U_i-1)T_i=(a_i+q_i)/(2sqrt(D_i))-1/2.

All square roots are positive. For finite positive X, 0<T_i<1 and mu_i>0. The quadratic root also gives

    mu_i(1+mu_i)=alpha_i T_i(1-T_i).

With the logarithmic derivative operator D=X d/dX, direct differentiation gives

    D T_i=T_i(1-T_i)(1+2mu_i),
    D mu_i=mu_i(1+mu_i)(1-2T_i).

Thus T_i increases strictly from 0 to 1. Put G=T_1+T_2-1. It increases strictly from -1 to 1 and has a unique zero X0. The source logarithmic derivative is

    F=X y'/y=mu_3-G.

Consequently F>0 whenever X<=X0, and every stationary point lies in X>X0, where 0<G<1.

On this interval, mu_3=G if and only if

    1/[G(1+G)] = D_3/(alpha_3 q_3)
                 = A X+B/X+C,

where

    A=K_3/alpha_3>0,
    B=a_3²/(alpha_3 K_3)>0,
    C=2/alpha_3-2.

Indeed mu_3(1+mu_3)=alpha_3 q_3/D_3 and u -> u(1+u) is injective for u>0. This equivalence involves no uncontrolled squaring or extraneous radical branches.

It follows that stationary points correspond exactly to solutions of

    L(X)=A X²+C X+B,  L(X)=X/[G(X)(1+G(X))], X>X0.

## 2. Universal derivative lemma

Claim: L'''(X)<0 whenever X>X0, for all positive parameters of ligands 1 and 2.

Fix such an X. Write t_i=T_i(X), m_i=mu_i(X), v_i=1+2m_i, b_i=t_i(1-t_i), h_i=1-2t_i. The logarithmic derivative identities yield

    P=DG=sum_i b_i v_i,
    Q=D²G=sum_i b_i h_i(3v_i²-1)/2,
    R=D³G=sum_i b_i v_i[(3v_i²-2)-b_i(15v_i²-9)].

For example D(b_i)=b_i h_i v_i and D(v_i)=(v_i²-1)h_i/2; these give the displayed formulas by the product rule.

Put H=G(1+G)>0 and A0=1+2G. Since L=X/H,

    X² L''' = D³(H^(-1))-D(H^(-1)).

Multiplying by H⁴ gives exactly

    N=A0 H²(P-R)+6H(1+3H)P Q-6A0(1+2H)P³.

To check its sign without parameter sampling, use positive coordinates

    x=(1-t_1)/G>0,  z=(1-t_2)/G>0,  S=1+x+z,
    G=1/S,  t_1=(1+z)/S,  t_2=(1+x)/S.

The script uses the letter y for z. The displayed map is a bijection from 0<t_1,t_2<1, t_1+t_2>1 to x,z>0: its stated inverse has G=1/S, and (1-t_1)/G=x, (1-t_2)/G=z. Thus no part of the positive G domain is removed. At any fixed X, the source parameters also allow arbitrary 0<t_i<1 and m_i>0, with alpha_i=m_i(1+m_i)/(t_i(1-t_i)) and q_i=(t_i+m_i)^2/[t_i(1-t_i)]; only the forward positivity is required for this sign proof. Let

    bbar_1=x(1+z), bbar_2=z(1+x),
    Hbar=S+1, Abar=S+2,
    Pbar=bbar_1 v_1+bbar_2 v_2,
    Qbar=[bbar_1(x-z-1)(3v_1²-1)
         +bbar_2(z-x-1)(3v_2²-1)]/2,
    Rbar=sum_i bbar_i v_i[(3v_i²-2)S²-bbar_i(15v_i²-9)].

Then b_i=bbar_i/S², H=Hbar/S², A0=Abar/S, P=Pbar/S², Q=Qbar/S³, R=Rbar/S⁴, and therefore

    S⁹ N = Abar Hbar²(Pbar S²-Rbar)
           +6Hbar(S²+3Hbar)Pbar Qbar
           -6Abar(S²+2Hbar)Pbar³.

Finally substitute v_1=1+2m_1 and v_2=1+2m_2. Exact polynomial expansion has 322 nonzero coefficients, all strictly negative, in the four positive variables x,z,m_1,m_2. The complete coefficient list is coefficient-certificate.json. A smaller display factors out -12(x+1)(z+1) and groups the remaining positive polynomial by m_1,m_2; see coefficient-table.md. The identity is reproducible by verify_peak_cap.py. Every monomial is positive on the stated domain, so N<0. Since X²H⁴>0, L'''<0, proving the lemma.

## 3. Rolle's theorem and the cap

The difference L(X)-(A X²+C X+B) has strictly negative third derivative on (X0,infinity). Four distinct zeros would yield a zero of its third derivative by three applications of Rolle's theorem, a contradiction. Hence y' has at most three distinct zeros on (0,infinity).

Three strict local maxima of the nonconstant real analytic y would require at least five distinct derivative zeros: one at each maximum and a minimum between each consecutive pair. Therefore y has at most two strict local maxima. This handles degenerate stationary points as well; no assumption of simple roots was made.

## Evidence and limitations

- The claimed theorem uses the full positive parameter domain; it does not restrict binding strengths, total concentrations, or equality cases.
- No floating-point search is used in the sign proof.
- The source baseline in source-baseline.json exactly certifies at least two strict maxima, showing sharpness. This checker artifact strengthens the source numerical example for use here; exactifying that example alone is not the proposed contribution.
- The symbolic expansion is a finite algebraic identity certificate, not a substitute for checking the derivative and stationary-point reductions.
- The helper's exploratory ratio and derivative probes are non-exhaustive floating evidence retained as research history. They are not premises of this proof.
