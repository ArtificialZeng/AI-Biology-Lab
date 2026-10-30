# Local nonlinear damping in the six-species clock

Owner analytic candidate. The original H6 stable-periodic-orbit question remains unresolved. This note proves a strictly identified local theorem, not absence of global periodic orbits or failure of the converter.

## Model, normalization and theorem

Use exactly the six-species clock in Hemery2605.23745v1 Eq.(4), with common catalytic rate alpha>0 and annihilation rate kappa>0. Every unordered annihilating pair is included once; at n=6 the six i,i+2 pairs are distinct. The deterministic mass-action ODE is

    T_i'=alpha*T_(i-1)-kappa*T_i*(T_(i-2)+T_(i+2)), i mod6. (1)

Let tau=alpha*t, z_i=2*kappa*T_i/alpha. Then

    dz_i/dtau=z_(i-1)-z_i*(z_(i-2)+z_(i+2))/2.             (2)

The uniform positive equilibrium is z_i=1, or T_i=alpha/(2*kappa).

**Local theorem.** For every alpha,kappa>0 this equilibrium is locally asymptotically stable, although its Jacobian has one purely imaginary eigenvalue pair. In the normalized Fourier coordinates below, a quadratic correction exposes the exact center cubic coefficient

    g=-3/14+i*sqrt(3)/7.

An explicit polynomial Lyapunov function certifies asymptotic stability on a nonempty positive neighborhood. Consequently there is a neighborhood containing no nonconstant periodic orbit. No source-figure initial condition is asserted to lie in that neighborhood. The theorem supplies a local damping result only.

## Fourier equations and exact cubic coefficient

Put w_i=z_i-1 and zeta=exp(i*pi/3). Use the DFT convention

    w_i=sum_{k=0}^5 c_k*zeta^(ki),
    c_k=(1/6)sum_i w_i*zeta^(-ki).

Reality gives c_0,c_3 real, c_5=conj(c_1), c_4=conj(c_2). Equation(2) gives exactly

    c_k'=lambda_k*c_k - sum_{j=0}^5 cos(2*pi*j/3)*c_j*c_(k-j mod6),
    lambda_k=zeta^(-k)-1-cos(2*pi*k/3).                     (3)

Writing omega=sqrt(3)/2, these eigenvalues in order are

    -1, -i*omega, -1-i*omega, -3, -1+i*omega, i*omega.

Hence the equilibrium is not hyperbolic and a conclusion based only on linearization would be insufficient.

Set A=c_1 and

    H=(1/2)/(2*lambda_1-lambda_2)=2/7+i*sqrt(3)/7,
    y_0=c_0-|A|^2, y_2=c_2-H*A^2, y_3=c_3.              (4)

Let Y=(y_0,y_3,Re(y_2),Im(y_2)), with its ordinary real Euclidean norm. This is an invertible real polynomial change of coordinates near zero (indeed triangular with an explicit polynomial inverse). Direct substitution in(3) gives

    A'=lambda_1*A+g*A*|A|^2+R_1(A,Y), g=H-1/2,
    Y'=L_s*Y+R_2(A,Y),                                  (5)

where the symmetric part of L_s has eigenvalues -1,-3,-1,-1. All pure-center terms of stable-coordinate degree zero in R_2 have center degree at least three. Every other R_2 term either has both center and stable degree at least one, or stable degree at least two. R_1 has only these latter two types. There is no unremoved pure-center quadratic or cubic remainder in R_1. For example its exact complex remainder is

    R_1=-y_0*A/2+y_2*conj(A)
        -y_3*(conj(H)*conj(A)^2+conj(y_2))/2.            (6)

This quadratic normalization does not assert that Y=0 is an invariant center manifold. In particular its generally nonzero stable remainder is retained in the Lyapunov proof.

## Explicit remainder bound

On r=|A|<=1 and s=||Y||<=1, the finite polynomials in(5) satisfy

    |R_1|<=C*(r*s+s^2),
    ||R_2||<=C*(r^3+r*s+s^2), with C=61.                 (7)

The attached exact certificate supplies all six real remainder polynomials. Here is its complete finite-bound rule: expand each real component as a polynomial in ReA,ImA,y_0,y_3,Rey_2,Imy_2. A monomial of center degree a and stable degree b has absolute value at most r^a*s^b. The degree categories following(5) bound each by the relevant expression in(7) on the unit box. Every coefficient is q_0+q_1*sqrt(3) with rational q_0,q_1; use |q_0|+2|q_1| as a rational upper bound. The sum of these coefficient bounds over all components is at most61, so it bounds each component-vector norm via its l1norm. This is exact coefficient enumeration, not a numerical fit. The certificate checks each term's degree and rational-coefficient decomposition.

## Strict Lyapunov certificate

Choose K=34730, so K>=28*C^2/3, and define

    V=|A|^2+K*||Y||^2.

From(5)–(7) and Re(g)=-3/14,

    V' <= -3*r^4/7-2*K*s^2
          +2*C*r^2*s+2*C*r*s^2
          +2*K*C*r^3*s+2*K*C*r*s^2+2*K*C*s^3.          (8)

Use

    2*C*r^2*s <= 3*r^4/28+(28*C^2/3)*s^2,
    2*K*C*r^3*s <= K*s^2/2+2*K*C^2*r^6.

Whenever r,s<=delta and

    delta<=1/16,
    delta<=K/[8*C*(1+2*K)],
    delta^2<=3/[56*K*C^2],                              (9)

the remaining r*s^2,s^3 terms total at most K*s^2/4, while2*K*C^2*r^6<=3*r^4/28. Therefore

    V'<=-3*r^4/14-K*s^2/4<0 unless A=Y=0.               (10)

The rational choice delta=1/100000 satisfies all three inequalities exactly. Initial V<delta^2 keeps r,s<delta, so the strict derivative prevents exit. Inverse coordinate transformation makes each closed smaller V-sublevel compact. It also gives |w_i|<=6*delta+delta^2<1 there, so all z_i remain strictly positive. Polynomial local existence plus boundedness ensures forward existence inside the sublevel. Strict Lyapunov decrease and invariance yield convergence to the only zero-derivative point, A=Y=0, and Lyapunov stability. This proves local asymptotic stability directly; no unproved center-manifold closure or numerical Floquet criterion is used.

Changing to physical time multiplies V' by alpha>0. Thus the theorem holds for all common positive alpha,kappa, including the exact source rates1,1000. The rescaled neighborhood is small and is not claimed to include the source initial state with one species concentration1 and the rest zero.

## What remains open and what the pilot means

H6 asks for a nonconstant positive orbitally asymptotically stable periodic orbit anywhere, with six cyclically ordered peaks. Such an orbit could lie outside the local stable neighborhood, so H6 is unresolved. Global convergence, the full source-figure trajectory's limiting behavior, all n-clock claims, module loading, converter correctness and the halter are not proved here. The source itself asserts oscillatory numerical behavior, not asymptotic attraction to a cycle; this theorem is not a counterexample to an author stability claim.

The bounded numerical pilot used30.66CPU seconds and observed slowly changing peaks rather than a closed stable-return certificate. Its approximate transverse multiplier above one was computed on a nonclosed trajectory and cannot establish periodic-orbit instability. A separate exact rational point refutes one proposed global entropy-Lyapunov inequality; that is a method-specific failed route, not a disproof of H6. Those records are retained separately and are not used as proof of the local theorem.

The 2017 four-species signed sine/cosine compiler is a known different construction and must be credited. Current source/overlap checking found no inspected six-species nonlinear stability theorem, but priority remains unverified. Fresh mathematical and contribution review is required; exact algebra, separate AI review and a successful build are not external human confirmation or formal proof-assistant verification.
