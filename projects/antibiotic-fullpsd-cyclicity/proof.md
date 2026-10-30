# Fixed-role obstruction in three-strain full-PSD antibiotic communities

Shawn Zeng (Zijian Zeng), AI Biology Lab, UCSI University.

A proposed bounded proof/code note, 28 October 2026. External author review and priority confirmation are pending.

## Statement and source scope

Use the three-strain well-mixed mixed-inhibition map of Athreya, Gokhale and Verma (2025), DOI10.1371/journal.pcbi.1013763, Eqs1,2,4–6. For any finite M>=1 every channel column contains exactly one P, one S, one D. There is no R. Let K_P,K_D>0 be common channel strengths and r_0,r_1,r_2>0 constant effective growth factors. On the interior of the simplex,

p_ij=S_ij exp(-K_D sum_k D_kj x_k)(1-exp(-K_P sum_k P_kj x_k)),
f_i=r_i product_j(1-p_ij),
F_i=x_i f_i/(sum_k x_k f_k).

The positive r_i condition includes the original positive-growth parameterization with common basal growth and ordered channel phenotype costs. The stronger impossibility statement does not rely on choosing nonphysical rates.

A cyclic pair means two distinct column words related by one or two cyclic rotations of the three strain entries. Repeated identical columns do not count. Absence of a cyclic pair is invariant under relabeling strains and channels.

**Claim.** If the phenotype matrix has no cyclic pair, F has no locally asymptotically stable interior fixed point. In the two-column-type cases with a common producer or common degrader, every interior fixed point has a real multiplier greater than1.

This resolves only a natural restricted full-PSD class. It does not resolve the source conjecture for other columns, resistance, more strains, periodic coexistence, or spatial dynamics. It does not claim cyclic-pair sufficiency. A topology containing an SSP column, for example, is outside this theorem.

## 1. Complete support reduction, independent of M

Identify the six column permutations of (P,S,D) with S_3. Multiplication by a three-cycle partitions them into two classes of three, equivalently the two permutation parities. Any two distinct columns in the same class are nontrivial cyclic rotations. Therefore no cyclic pair permits at most one distinct column from each class: at most two distinct types, however many times they are repeated.

If there is just one type, the producer and degrader are immune to every channel. If there are two types, their relative permutation is odd and therefore a transposition. It fixes exactly one entry position, so one strain has the same phenotype in both types and hence in every channel. Its common phenotype is P, S, or D. This exhausts all supports: six singleton supports and nine cross-parity two-type supports. No finite M enumeration is being used as an infinite proof.

## 2. Two immune strains rule out attraction

If two strains a,b are never sensitive, f_a=r_a and f_b=r_b everywhere. An interior fixed point requires all fitnesses to coincide. Thus r_a!=r_b excludes an interior fixed point. If r_a=r_b, the ratio F_a/F_b=x_a/x_b is exactly invariant at every iterate. Every neighborhood of an interior fixed point contains points with a different ratio, which cannot converge to that fixed point. Hence it is not locally asymptotically attracting.

This handles a singleton support and the common-S two-type case, whose other two strains alternate P and D. The invariant argument is essential: a unit linear multiplier alone would not exclude nonlinear attraction.

## 3. Exact coordinate test when strain0 is immune

In common-P and common-D cases strain0 is immune, so f_0=r_0. Put y=(x_1,x_2), x_0=1-y_1-y_2 and z_i=log(x_i/x_0). The log-ratio map is exactly

z_i^+=z_i+log f_i(x)-log r_0.

Write G=partial(log f_1,log f_2)/partial(y_1,y_2). The inverse logistic chart has derivative

B=partial y/partial z=diag(y)-y y^T,
det B=y_1 y_2 x_0>0.

Thus at a fixed point, the simplex derivative is conjugate to J_z=I+GB, and

det(J_z-I)=det G det B.

If detG<0, GB is a real2x2 matrix with negative determinant. Its two eigenvalues are real and have opposite signs. A positive eigenvalue of GB gives an eigenvalue of J_z strictly greater than1. Since the map is analytic on the interior, the fixed point is unstable and cannot be locally asymptotically stable.

## 4. Common producer

After relabeling, the two types are (P,S,D) repeated m_1>=1 times and (P,D,S) repeated m_2>=1 times. Write a=K_P,b=K_D and h=exp(-a x_0). Let

u_1=1-exp(-b y_2)(1-h),
u_2=1-exp(-b y_1)(1-h).

Then f_1=r_1 u_1^(m_1), f_2=r_2 u_2^(m_2). All u_i are strictly positive. Define k_1=exp(-b y_2),k_2=exp(-b y_1). Direct differentiation with x_0=1-y_1-y_2 yields

partial_1 u_1=a k_1 h,
partial_2 u_1=a k_1 h+b k_1(1-h),
partial_2 u_2=a k_2 h,
partial_1 u_2=a k_2 h+b k_2(1-h).

Consequently G has the form

G=[[A,A+B],[D+C,D]],

where A=m_1 a k_1 h/u_1, B=m_1 b k_1(1-h)/u_1,
D=m_2 a k_2 h/u_2, C=m_2 b k_2(1-h)/u_2 are all strictly positive. Hence

detG=A D-(A+B)(D+C)=-A C-B D-B C<0.

## 5. Common degrader

Now the two types are (D,S,P) repeated m_1>=1 times and (D,P,S) repeated m_2>=1 times. Put k=exp(-b x_0),h_1=exp(-a y_2),h_2=exp(-a y_1),

u_1=1-k(1-h_1),
u_2=1-k(1-h_2).

The same f_i=r_i u_i^(m_i) representation holds. Derivatives are

partial_1 u_1=-b k(1-h_1),
partial_2 u_1=-b k(1-h_1)-a k h_1,
partial_2 u_2=-b k(1-h_2),
partial_1 u_2=-b k(1-h_2)-a k h_2.

Thus G=[[-A,-(A+B)],[-(D+C),-D]], with all A,B,C,D>0, giving the same strictly negative determinant. Section3 excludes stable interior fixed points.

The support reduction and the three common-role cases prove the claim for every finite M and all stated positive parameters. No numerical approximation, finite scan, cost fitting, or source catalogue-membership assumption is needed for this implication.

## Meaning and remaining unknowns

The mechanism is a fixed-role obstruction: without two co-oriented PSD channels, either two immune strains retain an incompatible or invariant abundance ratio, or cross-feedback dominates the diagonal response and forces an unstable multiplier. This supports a bounded necessity statement in the original population map, not a clinical or experimentally established microbial-design rule.

The broader source conjecture remains unresolved. A bounded primary-source screen found no equivalent proof in the inspected sources; it does not establish priority. Separate fresh AI reconstruction checked the argument, the native derivatives and supporting examples. This is not external human peer review, a proof-assistant certificate, or a clinical/experimental result. The accompanying script checks exact symbolic identities and finite examples; it does not mechanically prove the analytic theorem.

## Two exact illustrative interior fixed points

These demonstrate that the unstable-multiplier argument applies to nonempty parts of the source's equal-basal-growth, strictly ordered-cost domain. They do not demonstrate stable coexistence. In both cases x=(1/3,1/3,1/3), g=1 and c_R=1/20.

* Common P, columns PSD and PDS: c_P=1/5, c_D=1/10, K_P=3 log(9/5), K_D=3 log(4/3). The effective rates are (3/5,9/10,9/10). Both sensitive survival factors equal 2/3, so all three fitnesses equal 3/5.
* Common D, columns DSP and DPS: c_P=3/20, c_D=1/10, K_P=3 log(17/15), K_D=3 log 2. The rates are (4/5,17/20,17/20), and sensitive survival is 16/17. All fitnesses equal 4/5.

Their numerical simplex multipliers are approximately (1.292858,0.856159) and (1.058670,0.951562), respectively. Those approximations are illustrations; instability follows from the analytic determinant argument. Other fixed points covered by the theorem need not be saddles: one multiplier greater than1 suffices for instability, regardless of the other multiplier.
