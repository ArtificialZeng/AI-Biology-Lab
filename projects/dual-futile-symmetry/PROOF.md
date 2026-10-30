# An exact stability classification on a symmetry slice

Status: exact restricted result with a separate AI mathematical check, 26 October 2026. Novelty has not been established. This is **not** a proof of the unrestricted conjecture.

The source model is Nicola Vassena, *Sequential and distributive dual futile cycle: Hopf bifurcation can occur under parameter-rich kinetics but cannot occur under mass action kinetics*, [arXiv:2608.27081v2](https://arxiv.org/html/2608.27081v2), especially Sections 2–3 and Conjecture 3.1. The upstream [MATLAB repository](https://github.com/nvassena/dualfutile_hopfnohopf) is credited for the model and published computations. V2 was revised on 14 September 2026. The source proves the separate mass-action no-Hopf result; that statement must not be advertised as still open.

## Domain and scope

Write the reduced coordinates as `(Q1,Q2,U,V,Y,X)`. In the source's 16-parameter positive Jacobian family impose

\[
a=c,\quad b_4=b_{10}=b,\quad e_1=f_7=e,\quad e_4=f_{10}=f,
\quad u_2=x_8=u,\quad u_3=x_9=w,\quad v_5=y_{11}=v,
\quad v_6=y_{12}=z.
\]

All eight remaining variables are strictly positive. These are paired **Jacobian parameters**, not a claim that any eight microscopic rate constants can be assigned arbitrarily in a mass-action realization. Their pairing makes the matrix equivariant under the source graph's involution `(Q1 Q2)(U X)(V Y)`. A graph automorphism alone does not imply this kinetic equivariance. We impose eight additional equalities, so the result is on a lower-dimensional slice, with no assertion about a neighborhood of that slice or about a general asymmetric model.

Put

\[
h=a+e+u+w,\quad j=2b+f+v+z,\quad k=f+v+z,
\quad \Delta=w(f+v+z)-fz.
\]

The specialized source matrix is

\[
G=\begin{pmatrix}
0&0&-w&0&z&0\\
0&0&0&z&0&-w\\
a&0&-h&-e&0&0\\
-b&-b&-f&-b-f-v-z&-b&0\\
-b&-b&0&-b&-b-f-v-z&-f\\
0&a&0&0&-e&-h
\end{pmatrix}.
\]

**Proposition.** Counting algebraic multiplicities, the spectrum of this matrix has:

| Condition | Right half-plane | Left half-plane | Imaginary axis |
|---|---:|---:|---:|
| \(\Delta>0\), equivalently \(\det G>0\) | 0 | 6 | 0 |
| \(\Delta=0\), equivalently \(\det G=0\) | 0 | 5 | one simple zero |
| \(\Delta<0\), equivalently \(\det G<0\) | 1 (a simple positive real eigenvalue) | 5 | 0 |

Thus Conjecture 3.1 holds on this slice. No member of the slice has a nonzero purely imaginary eigenvalue. This excludes a Hopf bifurcation while staying inside this family; it does not assert absence of periodic orbits or extend the no-Hopf statement to arbitrary parameter-rich kinetics.

## Proof

Use the even basis `Q1+Q2, U+X, V+Y` and odd basis `Q1-Q2, U-X, V-Y`. The change-of-basis matrix has nonzero determinant (the checker computes it exactly), and direct multiplication gives

\[
P^{-1}GP=\operatorname{diag}(A,B),\qquad
A=\begin{pmatrix}0&-w&z\\a&-h&-e\\-2b&-f&-j\end{pmatrix},
\quad B=\begin{pmatrix}0&-w&-z\\a&-h&-e\\0&-f&-k\end{pmatrix}.
\]

For a monic cubic \(t^3+c_1t^2+c_2t+c_3\), the strict Routh–Hurwitz criterion is \(c_1,c_2,c_3>0\) and \(c_1c_2-c_3>0\).

The coefficients for the even block are

\[
c_{1A}=h+j,\qquad c_{2A}=hj-ef+aw+2bz,
\]
\[
c_{3A}=awj+2bew+afz+2bzh>0.
\]

In particular \(hj-ef>0\): writing \(h=\alpha+e\), \(j=\beta+f\), where \(\alpha=a+u+w>0\) and \(\beta=2b+v+z>0\), gives \(hj-ef=\alpha j+e\beta\).

For completeness, the following exact identity exhibits positivity of the even Routh expression:

\[
\begin{split}
H_A={}&c_{1A}c_{2A}-c_{3A}\\
={}&\alpha^2j+e\alpha j+e^2\beta+ej\beta+awh+2bzj\\
&+e\big[(a+u)\beta+w(v+z)\big]
+(u+w)j^2+a(j^2-fz)>0.
\end{split}
\]

Every summand is positive: \(j^2-fz\) has strictly positive coefficients when expanded because \(j=2b+f+v+z\). Thus \(A\) is always Hurwitz.

The odd coefficients are

\[
c_{1B}=h+k,\quad c_{2B}=hk-ef+aw>0,\quad c_{3B}=a\Delta.
\]

Here \(hk-ef=\alpha k+e(v+z)>0\), and

\[
H_B=c_{1B}c_{2B}-c_{3B}
=(h+k)(hk-ef)+awh+afz>0.
\]

If \(\Delta>0\), the cubic criterion makes \(B\) Hurwitz. If \(\Delta=0\), its polynomial is \(t(t^2+c_{1B}t+c_{2B})\). Its zero is simple because \(c_{2B}>0\), and the other roots have strictly negative real parts.

If \(\Delta<0\), the coefficient signs of the odd cubic are \((+,+,+,-)\). It has exactly one positive root \(r\); indeed its derivative is strictly positive on \([0,\infty)\), and its values at zero and infinity have opposite signs. The other two roots have sum \(-c_{1B}-r<0\) and product \((-c_{3B})/r>0\). They therefore have negative real parts, whether real or complex conjugates. The positive root is simple. Finally,

\[
\det G=c_{3A}c_{3B}=c_{3A}a\Delta,
\]

since both characteristic factors have odd degree and \(G\) has even dimension. The asserted classification follows. □

## Exact reproduction

Run from this project directory with Python 3. No third-party package is required:

```sh
python3 verify_symmetry_slice.py --self-test
python3 -O verify_symmetry_slice.py --self-test
```

The checker uses its own sparse polynomial arithmetic over rational coefficients and a direct Leibniz determinant. It verifies the basis identity, the complete 6×6 characteristic-polynomial factorization, all coefficient formulas, both grouped positivity identities, and three rational examples. The expanded Routh expressions contain 86 and 56 monomials respectively, all with positive coefficients. Its three deliberately corrupted inputs must be rejected; checks do not use Python `assert`. The companion SymPy derivation provides a different algebra implementation for comparison. Receipts are computational verification, not external expert review or a formal proof-assistant certificate.

## Priority check and remaining work

The source latest version expressly states the general conjecture remains unproved. A bounded search on 26 October 2026 for symmetric dual-futile-cycle stability, parameter-rich symmetry and prior Hopf work did not locate this exact eight-parameter threshold. This does **not** establish novelty. Earlier stability and bistability literature, including Hell–Rendall's [*A proof of bistability for the dual futile cycle*](https://arxiv.org/abs/1404.0394) and Conradi–Feliu–Mincheva's [2020 Hopf paper](https://doi.org/10.3934/mbe.2020027), must be compared at the level of their kinetic assumptions before claiming a new theorem. Ask the source author whether this elementary symmetry corollary is already known before describing it as original research.
