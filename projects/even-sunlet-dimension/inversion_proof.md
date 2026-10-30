# Cyclic-inversion 3-sunlet dimensions

## Statement and scope

For every integer q>=5, let G=C_q and let B={id,inversion} act on G.
Put m=floor(q/2), k=m+1. For the literal group-based 3-sunlet model with
edge parameters identified on B-orbits, the complex affine cone has dimension

    dim V_N^(G,B) = 5k-4 = 5 floor(q/2)+1.

This is a theorem for this natural nontrivial-B subclass. It leaves the
broader nontrivial-B dimension problem unresolved, and does not assert
topology identifiability or a biological experimental conclusion. The source
model, dimension bounds and earlier general-model results are credited below.

## Literal coordinates and dimension bounds

Represent G by integers modulo q, and write

    d(x)=min(x mod q,(-x) mod q),  0<=d(x)<=m.

The B-orbits have labels0,...,m. There are6k independent parameters w_i,a.
Retain all q^2 consistent leaf-label coordinates (g,h,-g-h), as in the
Gross–Krone–Martin2024 model definition. Their two literal summands are

    p_g,h = w1_a w2_b w3_c (w4_c w5_a + w4_b w6_a),
    a=d(g), b=d(h), c=d(g+h).

The image closure over C is the variety in question. Equal output polynomials
are retained as duplicate coordinates; they do not affect dimension. We
never assume that two orbit labels determine a third orbit label or that
the codomain has dimension k^2.

For a direct upper bound, absorb w1_a into alpha_a=w1_a*w5_a and
beta_a=w1_a*w6_a. At a=0, g=0 and b=c, so only gamma0=alpha0+beta0 occurs.
On the dense open set where w2_0,w3_0,w4_0 are nonzero, divide the three
entire edge blocks by those zero coordinates and absorb their product into
alpha,beta,gamma0. The resulting normalized parameter count is

    3(k-1)+2(k-1)+1 = 5k-4.

The dense-open image has the same closure as the full image. Thus the
dimension is at most5k-4; this also agrees with the source upper bound.

For the lower bound we use a strict selected exponent matrix A. One can
check the source tropical criterion directly: put w_i,a=t^lambda_i,a in
the logarithmic Jacobian. Factoring the unique smallest power from each
output column leaves its selected exponent vector as the constant term.
If a minor of A is nonzero, that determinant is the nonzero leading
coefficient of the corresponding Jacobian minor. Hence the polynomial
Jacobian has rank at least rank A at some nonzero parameter point. In
characteristic zero its generic rank equals the dimension of the image
closure. Therefore dim V is at least rank A. Strictness is essential.

## A strict orbit-constant chamber

For orbit labels a=0,...,m choose actual positive integer weights

    lambda1_a=lambda2_a=lambda3_a=1,
    lambda4_a=4a+1,
    lambda5_a=2,
    lambda6_a=2+mu_a,
    mu0=-1, mu_a=1 for a>0.

Every weight is constant on its B-orbit. The first-minus-second literal
monomial weight is exactly

    4(c-b)-mu_a.

At g=0, a=0 and c=b, giving+1, so every such column selects the second
summand. At g!=0, mu_a=1, so the difference is a multiple of4 minus1.
There is no tie: the first term is selected precisely when c<=b, and
the second precisely when c>b. Denote the6k-by-q^2 selected exponent
matrix by A. All comparisons above hold uniformly in q.

## Reduce the entire left kernel

Let (A_a,B_a,C_a,D_a,E_a,F_a), a=0,...,m, be an arbitrary row functional
annihilating every column of A. These are coefficients of a row functional,
not values of the original parameters. From g=0 second-summand columns,

    A0+B_b+C_b+D_b+F0=0  for every b.

Put

    K=A0+F0, U_b=C_b, W_b=C_b+D_b,
    r_a=A_a+E_a-K, s_a=A_a+F_a-K.

Then B_b=-W_b-K. A first-summand column is annihilated exactly when

    r_a+W_c-W_b=0,                         (I)

and a second-summand column exactly when

    s_a+U_c-U_b=0.                         (II)

Both statements follow by substituting the five literal factors. No
equation or sign condition is discarded.

## Step1 makes both functions affine in orbit distance

Let rho=r1 and tau=-s1. For each j=1,...,m choose g=1 and h=-j modulo q.
Then b=j,c=j-1, so the column selects the first term. Equation(I) gives

    W_j-W_(j-1)=rho.

Choose instead g=1,h=j-1. Now b=j-1,c=j, so the column selects the
second term and(II) gives U_j-U_(j-1)=tau. Consequently

    W_j=W0+j rho, U_j=U0+j tau,  0<=j<=m.   (A)

This uses only actual labels and applies to both parities, including the
largest orbit label.

## Odd q

Let q=2m+1, so q>=5 means m>=2. At g=1,h=m, one has b=c=m. This tie
in orbit distance still selects the first summand, because its weight
difference is-1. Equation(I) forces r1=0. Thus rho=0 and W is constant.

For U choose g=2. Its orbit label is2, also when q=5. At h=0 one has
b=0,c=2, so(II) and(A) give s2+2tau=0. At h=m-1 one has b=m-1 and
c=d(m+1)=m, so this is also a second-summand column, giving s2+tau=0.
Subtracting yields tau=0. Hence U is constant. All named labels exist
already at the smallest odd case m=2.

## Even q

Let q=2m, so q>=5 means q>=6 and m>=3. At g=2,h=q-1 one has b=c=1,
so equation(I) gives r2=0. At g=2,h=q-2 one has b=2,c=0, also a
first-summand column. Equations(I) and(A) give r2-2rho=0, hence rho=0
in characteristic zero. Therefore W is constant.

For U choose g=3. Its orbit label is3 because m>=3, including the
self-inverse element3 at q=6. At h=0 one has b=0,c=3, giving
s3+3tau=0 from(II). At h=q-1 one has b=1,c=2, giving s3+tau=0.
Both select the second summand. Their difference is2tau=0, so tau=0
in characteristic zero. Thus U is constant, including q=6.

## Full kernel converse and rank

Write C_a=C and D_a=D for the resulting constants. Then B_a=B is
constant and K=-B-C-D. For any nonzero orbit a, choose g=a. The column
h=0 has b=0,c=a>0 and selects the second term. The column h=-g has
b=a,c=0 and selects the first. Therefore(II) and(I) respectively force

    E_a=F_a=K-A_a,  a>0.

At zero, F0=K-A0 by definition, while E0 is unused because every g=0
column selected the second summand. Conversely choose arbitrary values
of all k coefficients A_a, arbitrary constants B,C,D and arbitrary E0.
Set K=-B-C-D, define F0 and all nonzero E_a,F_a by these formulas, and
take each B_a=B,C_a=C,D_a=D. Direct substitution annihilates every
selected column. Each of these k+4 choices is recoverable individually
from the row functional, so they are independent and describe the
entire left kernel. Thus

    rank A = 6k-(k+4) = 5k-4.

The matching upper and lower bounds prove the stated theorem for every
q>=5. The proof does not infer an infinite quantifier from finite tests.

## Diagnostics, boundaries and contribution

The exact checker constructs both literal monomials for all q^2 labels,
checks actual positive weight sums and all branch decisions, verifies
an explicit k+4 integer kernel and its independence, and computes exact
rational selected ranks. It also evaluates the ordinary Jacobian of the
sum of both terms at three serialized integer points per case. Held-out
q11,12,17,18 give ranks26,31,41,46. q5,...,10 reproduce the triage ranks
11,16,16,21,21,26. These computations are finite diagnostics only.

The q3/4 controls have only5/10 distinct output polynomials and selected
ranks5/10, below the formal5k-4 counts6/11. They are outside the theorem
and not new small-order dimension claims. The odd slope argument needs
m>=2; the even argument needs m>=3, matching the stated domain.

The new component is a strict orbit-distance chamber and a uniform
kernel theorem for cyclic inversion constraints, a canonical nontrivial-B
class. It does not follow merely by imposing those constraints on the
completed B={id} theorem. The source all-B triangle-free results exclude
three-leaf cycles, and the source2025 discussion leaves nontrivial B
unexplored. Known general-model finite examples and the R11/R12 proofs
are credited prior work. Bounded primary comparison found no equivalent
uniform cyclic-inversion theorem, without certifying priority.

No gap is currently known in the stated characteristic-zero argument.
Fresh mathematical/contribution review is required. No trusted-runner,
proof-assistant, external human or priority certification is asserted.
