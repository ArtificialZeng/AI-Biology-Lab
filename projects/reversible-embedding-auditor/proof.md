# Exact four-state counterexample to Theorem 3.16 (1) ⇒ (3)

This is a counterexample to the stated implication in arXiv:2511.22800v1, not
a claim that every result in that paper is false. It uses the source's row
convention. The input matrix is rational, strictly positive, symmetric and
therefore reversible for p=(1/4,1/4,1/4,1/4). Its spectrum is strictly positive.

Put h=9 log 2 and define

    A = (1/36) [-63  33  33  -3
                33 -79  29  17
                33  29 -79  17
                -3  17  17 -31],

    B = (1/12) [ 9 -21   9   3
                 9  -5  -7   3
               -21  33  -5  -7
                 3  -7   3   1].

Then Q=hA+πB is a real Markov generator. Both A and B have zero row and
column sums. Its twelve off-diagonal entries, with repetitions suppressed, are

    (11h−21π)/12, (11h+9π)/12, (−h+3π)/12,
    (29h−21π)/36, (17h+9π)/36,
    (29h+99π)/36, (17h−21π)/36.

Every expression is strictly positive. Indeed log 2>2/3 and log 2<1 give
6<h<9, while 3<π<22/7 gives 11h>66>21π and 3π>9>h. The other
inequalities follow immediately. The inequality log 2>2/3 follows, for
example, from log 2=2[1/3+(1/3)^3/3+…]. No rounded constant is used.

For direct verification of its exponential, take the orthonormal vectors

    e=(3,3,3,3)^T/6,
    v=(−3,−1,5,−1)^T/6,
    w=(−3,5,−1,−1)^T/6,
    u=(−3,−1,−1,5)^T/6.

Set E=ee^T, P=vv^T+ww^T, R=uu^T and J=2vw^T−(1/2)wv^T.
These are mutually orthogonal projectors E,P,R with E+P+R=I, and

    A=−R−3P, B=2J, J²=−P, JE=JR=0, JP=PJ=J.

Consequently A commutes with J. On ran(P), exp(2πJ)=I because J²=−I;
on its complement J=0. Therefore

    exp(Q) = exp(hA) exp(2πJ)
           = E+exp(−h)R+exp(−3h)P
           = E+2^−9 R+2^−27 P =: M.

Explicitly the rational input matrix is

    M = 1/268435456 *
      [67239937 67152554 67152554 66890411
       67152554 67123429 67123427 67036046
       67152554 67123427 67123429 67036046
       66890411 67036046 67036046 67472953].

The displayed entries and row sums independently verify strict positivity and
stochasticity. Symmetry verifies detailed balance with strictly positive p.
The orthogonal decomposition gives eigenvalues 1,2^−9,2^−27,2^−27, all
positive. Thus the principal matrix logarithm is

    Log(M) = −hR−3hP = hA,

whose entry (1,4) is −h/12=−3 log(2)/4<0. It is not a Markov generator.
All the hypotheses of Theorem 3.16(1) hold, while (3) fails.

Since M is symmetric with positive spectrum, its only reversible logarithm
for the uniform p is its principal logarithm: any such logarithm is symmetric
and orthogonally diagonalizable with real eigenvalues; real scalar
exponentiation is injective. Hence this also separates arbitrary embeddability
from reversible embeddability for this M. This uniqueness step is standard
and is already part of the source's Lemma 3.10 / Jia's 2016 result.

The construction actually works for h in the open interval
21π/11<h<3π, replacing the dyadic eigenvalues by exp(−h),exp(−3h).
This is only explanatory context; the frozen decisive claim is the one
rational matrix above. No minimum-dimension, open-set, or new classification
claim is made.

Proof gap list: none known in the counterexample argument. Historical novelty
is a separate bounded literature check and remains subject to fresh review.
