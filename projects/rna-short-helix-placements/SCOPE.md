# Model, scope, and review status

Model: pseudoknot-free secondary structures on the full sequence, alphabet A,C,G,U; only A-U/C-G pairing; zero minimum span; energy minus number of pairs. Competing folds are all compatible noncrossing partial matchings. A maximal helix is a maximal uninterrupted stack of consecutive nested pairs. An isolated stack is a maximal helix of length 2.

The interval tree has one node per pair, singleton nodes for unpaired positions, and a virtual root. Paired degree includes the incoming paired incidence except at the root. Exclude m5 (paired degree >4) and m3dot (unpaired child and paired degree >2). A proper coloring uses B,W,G, with exposed incidences containing at most one B, one W, two G, where the incoming incidence is color-inverted (B/W swap, G fixed). Levels are root sums B=+1,W=-1,G=0. Modulo-2 separation means gray and unpaired level residues are disjoint.

Frozen broad original target: classify the source's three-or-more-short-helix structures. Source status: `open-supported` as of the inspected August 2026 preprint; follow-up searches on 26 October found no equivalent classification. The absence of a search hit is not proof of priority.

Admitted subsidiary claims:

1. Complete classification of motif-free targets whose entire maximal-helix set consists of exactly three length-2 helices, with arbitrary unpaired gap lengths. They fail modulo-2 separation exactly for two geometries described in `PROOF.md`. Under the gap-presence normalization there are 288 classes, 4 failures, 284 constructive successes.
2. As a short corollary of the source F/Q construction, every motif-free target with no length-1 helix, for which all length-2 helices are linearly ordered by ancestry in the contracted helix tree, admits a uniform modulo-2 separated coloring, regardless of their number. Other helices can have arbitrary lengths >=3 and branch arbitrarily.

These are newly formulated restricted questions, not a claim to have solved the full source question. Prior broad finite-state decision capability is already known. The main explicit refinement is the small complete success/failure placement classification. The unlimited-count comparable-short statement is a short corollary of the published F/Q construction. A separate AI review found no mathematical gap after requiring these attribution and status qualifications. Priority and independent human expert review remain pending. Neither the new statements nor the Python checks have been Lean formalized. No runner/release-gate acceptance is claimed.
