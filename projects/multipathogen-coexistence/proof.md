# Counterexamples to the simultaneous-endemic-equilibrium assertion

## 1. Full source model and admissible parameters

The source ODE is, for p=1,…,P and v=1,…,V,

dS_p/dt=b_p−d_p S_p−S_p sum_{q,v} beta_pqv I_qv,
dL_pv/dt=S_p sum_q beta_pqv I_qv−(epsilon_pv+d_p)L_pv,
dI_pv/dt=epsilon_pv L_pv−(gamma_pv+d_p)I_pv,
dR_p/dt=sum_v gamma_pv I_pv−d_p R_p.

Use P=2,V=3,k=1 and b_p=d_p=epsilon_pv=gamma_pv=1. All states are nonnegative. Every specified contact beta_pqv is strictly positive. Hence both host populations genuinely interact and every strain can transmit in either direction. Define N_p=S_p+sum_v(L_pv+I_pv)+R_p. Direct summation gives Ndot_p=1−N_p; therefore N_p=1 and 0<S_p<=1 at every nonnegative equilibrium. Strict positivity of S_p follows because S_p=0 would give Sdot_p=1.

## 2. Simple exact witness

Set beta_pqv=lambda_v for every p,q, with lambda=(1,3,4). At the disease-free state S_p=1, source(9) gives B_v=(lambda_v/4)11^T, with eigenvalues0 andlambda_v/2. Thus R01=1/2<1,R02=3/2>1,R03=2>1, exactly the strict hypotheses of Conjecture5.

Put J_v=I_1v+I_2v and T=S_1+S_2. For ANY state of the full16-variable ODE,

sum_p (Ldot_pv+2 Idot_pv)=J_v(lambda_v T−4).

At an equilibrium, a strain endemic in at least one host has J_v>0, so lambda_v T=4. Simultaneous endemicity of strains2 and3 would require T=4/3 and T=1. These are incompatible. Equivalently4(3T−4)−3(4T−4)=−4 gives a constant exact contradiction. No such mixed equilibrium exists, stable or unstable. This refutes the universal original conjecture, not merely its stability adjective.

For this witness all nonnegative equilibrium supports are completely classified: the empty support is the DFE; strain1 cannot be endemic because it would require T=4>2; every support with two different lambdas is impossible by the same identity. The two remaining equilibria are singleton strain2 or singleton strain3. Indeed all hosts have a common S=s=2/lambda_v because their S-equations have identical birth/death and force of infection. The I-equations and L-equations imply I_pv=t,L_pv=2t,R_p=t, and N_p=1 then gives t=(1−s)/4. Thus strain2 has (s,t)=(2/3,1/12), and strain3 has (s,t)=(1/2,1/8), in every host; all other infected components vanish. No continuum or omitted asymmetric equilibrium escapes this classification.

## 3. Full-rank positive witness

The simple witness already meets the source assumptions. To remove reliance on a rank-one contact matrix, replace the three matrices beta_v by

beta_1=[[1,1/2],[1/2,1]],
beta_2=[[4,1],[1,4]],
beta_3=[[5,2],[2,5]].

All entries are positive and all determinants are nonzero. Source reproduction numbers are respectively3/8,5/4,7/4 because B_v=beta_v/4. At an equilibrium, the I-equations give L_v=2I_v and the L-equations give I_v=(1/4)diag(S)beta_v I_v. A nonzero vector I_v therefore requires det(I−T_v)=0, where T_v=diag(S)beta_v/4.

Writing U=S_1S_2,T=S_1+S_2 gives

d2=det(I−T2)=1−T+(15/16)U,
d3=det(I−T3)=1−(5/4)T+(21/16)U.

If both strains2 and3 are endemic, d2=d3=0. The identity d3−(5/4)d2=−1/4+(9/64)U gives U=16/9; then d2=0 gives T=8/3. This contradicts T<=2 at every nonnegative equilibrium. This second exact witness again satisfies all strict original hypotheses while lacking the required mixed equilibrium.

## 4. Standard dominance explanation and robustness

This paragraph applies the classical Perron–Frobenius theorem; it is not asserted as a new biological principle. For arbitrary positive source rates, put c_pv=epsilon_pv/[(epsilon_pv+d_p)(gamma_pv+d_p)] and D_v=diag(c_v)beta_v. At any nonnegative equilibrium S_p>0 and I_v=diag(S)D_v I_v. If beta_v is strictly positive and I_v is nonzero/nonnegative, that equation forces I_v>0; Perron–Frobenius implies rho(diag(S)D_v)=1.

If D_w>D_v entrywise, then diag(S)D_w>diag(S)D_v entrywise, and strict Perron monotonicity makes their spectral radii unequal. Hence v,w cannot both be endemic. One direct proof of strict monotonicity is to take a positive right Perron vector x of the smaller positive matrix A and a positive left Perron vector y of the larger C: rho(C)y^Tx=y^TCx>y^TAx=rho(A)y^Tx.

In the full-rank witness c_pv=1/4 and beta3−beta2=11^T, so strict dominance holds. All inequalities of positivity, dominance and R01<1<R02,R03 are strict and continuous in the full collection of model parameters. They consequently hold on some open neighborhood of this witness. The nonexistence is not an isolated rank or exact-symmetry artifact. No global convergence or trajectory-persistence theorem is inferred.

## 5. Local outcome and controls for the simple witness

The full16D Jacobian at the strain3-only point has characteristic polynomial

(z+1)^3(z+2)^7(z+3)(z^2+4z+1)(z^3+6z^2+8z+4).

The linear factors have negative roots. The quadratic is Hurwitz because its two coefficients are positive. The cubic is Hurwitz because6,8,4>0 and6*8>4. Thus all16 eigenvalues have negative real part and the point is locally asymptotically stable. Invader reproduction numbers for strains1 and2 there are1/4 and3/4. These are evaluated at depleted susceptible levels, not at the DFE.

The checker also verifies an equal-supercritical-contact control: lambda=(1,4,4) permits positive coexistence, with each host S=1/2,L_2=L_3=1/8,I_2=I_3=1/16,R=1/8. This shows that the proof does not incorrectly forbid all multi-strain endemic equilibria. It checks all six strain relabelings and additional host-size singleton controls. These finite controls validate implementation scope; the displayed algebraic contradiction is the proof of nonexistence.

## 6. Scope and gaps

No mathematical gap is known in this argument, pending fresh independent review. The source Conjecture5, read as its literal universal simultaneous-endemic-equilibrium claim, is disproved. We do not classify general heterogeneous coexistence, claim competitive exclusion is new, assert global convergence, identify an author-confirmed error, or infer anything about experimental fish infections. Theorem4's DFE/endemic wording is not needed for the main disproof. A bounded current literature search found no correction to the exact conjecture in inspected sources; this is not universal priority certification.
