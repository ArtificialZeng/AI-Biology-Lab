# Model and current definition

Source status: **new restricted question**, not an explicit named open conjecture and not a claim that the source theorem is wrong. The source establishes generic nonzero circulation but does not supply the global N=3 sign/zero locus in the inspected text.

Fixed objects: sites I={0,1,2}; exactly N=3 particle slots. Each particle at i moves to j!=i at rate q(i,j)=1. It is killed at rate k_i>0 and immediately chooses one of the other two slots uniformly, copying both site and ancestral label. Copies at the same site remain part of the labelled dynamics.

Population counts c in Omega_3={c in N0^3:sum c=3}, ten states, use the row-generator convention

    L(c,c-e_i+e_j)=c_i + c_i*c_j*k_i/2  (i!=j),

with diagonal minus outgoing sum. pi is the unique stationary probability, pi L=0.

Founder h(i,d) is the probability the label of one specified founder at site i ultimately fixes, when the other two slots have counts d, sum d=2. This gives 3*6=18 unknowns. Equation A.4, valid at arbitrary positive k, defines the linear system T h=0:

    (Th)(i,d)=sum_(j!=i)[h(j,d)-h(i,d)]
       +sum_a d_a sum_(b!=a)[h(i,d-e_a+e_b)-h(i,d)] - k_i h(i,d)
       +sum_a d_a k_a [h(i,d-e_a+e_i)
          + (1/2) sum_b (d_b-1_(a=b)) h(i,d-e_a+e_b) - h(i,d)].

Terms with zero multiplying count are omitted before arguments are evaluated. Anchor h(0,(2,0,0))=1/3. Source Proposition3.2 guarantees the anchored solution is unique. The ten population normalizations sum_i c_i h(i,c-e_i)=1 can be checked independently.

Set

    F_ij=sum_(c:c_i>0) pi(c)*c_i*h(j,c-e_i),
    Delta_3(k0,k1,k2)=F_02-F_20.

Orientation is fixed: Delta>0 means 0-to-2 ancestral current exceeds 2-to-0. On the symmetric triangle the other oriented cycle currents equal it by conservation.

The target is the complete sign and zero locus on this entire positive domain. Primary sources and current review status are documented separately.

## Four-particle extension

Retain the same graph, current orientation, and same-site copying convention; set N=4, so the donor divisor is three. The candidate path is exactly k=(1,t,t+1), t>1. `four_particle_model.py` implements the source's general count/founder equations for this path, while `verify_four_particle_independent.py` independently enumerates victim/donor events on six colour/site categories. The 15 count states and 30 founder states, polynomial solutions, current formula, and phase polynomial are under `evidence/four-particle/`. See `four-particle-proof.md` for the exact threshold, minimality, and local geometry.

Only killing-rate order is fixed along the path; ratios change. The original broader questions about full N4 phase structure, fixed-ratio scale reversal, general graphs, and arbitrary populations remain unresolved.
