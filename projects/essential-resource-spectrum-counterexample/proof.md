# Exact construction and scope

The target is the real-negative-B stability conjecture stated before Fig. 6 in Gibbs et al. (2022), not their separate symmetric-matrix sufficient theorem. We retain their equal-resource/equal-consumer analytical setting and both conversion pathways of published Eq. (1). Matrices and rates are dimensionless.

For resource index i and consumer index k, write

```
g_k(R) = min_{j:C_jk != 0} epsilon_jk C_jk R_j
dR_i/dt = rho_i - R_i sum_k C_ik N_k
          + sum_j P_ji sum_k (C_jk R_j - g_k/epsilon_jk) N_k
          + theta sum_j Pt_ji sum_k (g_k/epsilon_jk) N_k
dN_k/dt = N_k [(1-theta)g_k - eta_k].
```

P and Pt are nonnegative ROW-stochastic matrices; their transposes occur in the resource equations. There is no resource washout term. Take diagonal efficiency e=1/63, every off-diagonal efficiency 1, theta=1/13, and

```
C  = [[7,13/2,1], [3/5,7,13/3], [1/6,1/3,7]]
P  = [[0,8/11,3/11], [1/2,0,1/2], [0,1,0]]
Pt = [[0,1/3,2/3], [1,0,0], [0,1,0]].
```

For every n>0 set R*=1, N*=n*1, eta_k=4/39 and
rho=n*(1043/90,16697/2970,4142/1485). Every inflow and consumption coefficient is strictly positive, while permitted conversion pathways include structural zeros. Both production matrices have zero diagonal and rows summing one. At R* each consumer's diagonal growth value is e*7=1/9, and every other candidate is larger by at least 1/18. Thus the minimum is unique in a neighborhood and g_k=e*7*R_k there. The consumer equations vanish since (1-theta)/9=4/39. Direct substitution in the original resource sums gives zero as well. Original biomass accounting provides the independent identity

```
d(sum_i R_i + 65 sum_k N_k)/dt = sum_i rho_i - 65 eta sum_k N_k.
```

## Local Jacobian and real spectrum of B

Let Cd=7, D=diag(C*1), Etilde have diagonal 1 and off-diagonal e, and beta=(1-theta)eCd=4/39. For the general uniform equilibrium R*=r*1,N*=n*1, differentiation gives

```
B = -C + P^T(C-Cd Etilde) + theta Cd Pt^T Etilde
A = -D + P^T D - Cd P^T Etilde + theta Cd Pt^T Etilde
J = [[n A, r B], [n beta I, 0]].
```

Our r=1 makes the supplement's r-scaled B and the main-text B identical. The positive inflow condition is rho=-rn B*1. Row stochasticity implies

```
1^T A = 1^T B = -kappa 1^T,
kappa=(1-theta)Cd(1+2e)=20/3.
```

Hence total resource and consumer perturbations form a quotient block with characteristic polynomial x²+(20/3)nx+(80/117)n, Hurwitz for every n>0. The sum-zero composition subspace is invariant; quotient and restriction characteristic factors multiply even without diagonalizability. Direct exact differentiation gives the A and B matrices in certificate.json, and

```
det(xI-B)=(x+20/3)[x²+(1334/195)x+105544993/9034740].
```

The quadratic has positive coefficients and discriminant 803167/11293425>0. Its roots are real and strictly negative, as is -20/3. B is independent of n.

## Complete open-interval classification

Computing det(xI-J), or independently det(x²I-nxA-n beta r B), gives the stable aggregate quadratic above times

```
q_n(x)=x⁴+a1 x³+a2 x²+a3 x+a4,
a1=5377n/195,
a2=n(793372103n+3169584)/4517370,
a3=1010399042n²/88088715,
a4=422179972n²/3435459885.
```

All a_i are positive for n>0. The Routh determinants are

```
Delta2=a1*a2-a3
      =n²(4265961797831n+6938862748)/880887150 > 0,
Delta3=a1*a2*a3-a3²-a1²*a4
      =73n⁴(29522765162582466287n-1639836184859336)/38798108551756125.
nc=1639836184859336/29522765162582466287.
```

The quartic Routh first column is (1,a1,Delta2/a1,Delta3/Delta2,a4). For 0<n<nc its signs are (+,+,+,-,+), so precisely two roots lie in the open right half-plane. For n>nc all are positive, so q_n is Hurwitz. No zero pivot or imaginary-axis root occurs on these open intervals: the coefficients, Delta2 and Delta3 are nonzero there. The aggregate factor is always Hurwitz. Therefore the full positive equilibrium is unstable below nc and locally asymptotically stable above nc. Physical feasibility and the strict limiter hold throughout both intervals, while B retains its same real negative spectrum. The counterexample is thus to the stated universal protection conjecture in the source's analytical setting.

The point n=nc is excluded from the strict classifications. Neither a nonlinear Hopf bifurcation nor a periodic orbit is proved. Decreasing n decreases all inflows in the fixed ratio displayed above; other network parameters and mortality are fixed. Instability lies in composition modes, while the aggregate block stays stable.

## Two consumers within the canonical class

This lemma permits arbitrary nonnegative row-stochastic P and Pt, including self recycling, but requires common diagonal consumption Cd>0, common diagonal efficiency e>0, off-diagonal efficiencies 1, 0<=theta<1, strict diagonal limiting, and positive uniform R*=r*1,N*=n*1. It makes no claim about unequal efficiencies or unequal equilibria.

For S=2 put kappa=(1-theta)Cd(1+e)>0 and beta=(1-theta)eCd>0. The aggregate quadratic x²+n kappa x+nr beta kappa is Hurwitz. On the invariant sum-zero line w=(1,-1), B acts by b. If its spectrum is real strictly negative, b<0. Write

```
A=B-L, L=(I-P^T)(diag(C*1)-C).
P=[[1-u,u],[v,1-v]], 0<=u,v<=1.
Lw=ell*w, ell=2(u*C12+v*C21)>=0.
```

Thus A acts on that line by a=b-ell<0. The remaining quadratic x²-na x-nr beta b has strictly positive nonleading coefficients and is Hurwitz. This reasoning permits repeated B eigenvalues and requires no diagonalizability. The one-consumer aggregate case is immediate. Three consumers are therefore the first possible size for this failure within this canonical class, not within the entire family of possible resource models.

## Interpretation and evidence limits

The source effective matrix B alone omits the resource-relaxation information carried by A. Here the full composition feedback depends on both. This illustrates established block-spectrum and quadratic-pencil mathematics, rather than establishing a new general spectral principle. Prior conceptual and ecological results, including Meacock et al. (2026 v2) and Liu et al. (2025), are credited in sources.md. The bounded new object is this physically feasible native-model counterexample with an exact feed interval and the scoped two-consumer obstruction.

No result about the separate symmetric sufficient theorem, nonlinear cycles, extinction, global behavior, empirical consortia, or general minimality is claimed. Independent human confirmation and exhaustive priority clearance remain pending. The accompanying exact checker is computational evidence and an inspectable derivation, not formal proof-assistant certification or human peer review.
