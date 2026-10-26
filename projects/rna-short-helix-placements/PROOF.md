# Structural placement results in the Watson–Crick pair-count model

A separate AI mathematical review found no gap in this stated scope; independent human expert review and priority confirmation remain pending. The comparable-short result below is a short corollary of the F/Q construction in Jogalekar, arXiv:2608.25194v1, §§5–8. The new statements and Python checks have not been Lean formalized, and this is not an accepted runner/release-gate record.

## Common local lemma

Work with relative residues: all unpaired nodes must have residue 0 and gray pairs residue 1. An actual root level can be translated to either relative residue by choosing which absolute parity is assigned to unpaired nodes.

For a helix subtree define an entry state by its parent's residue and the head color. Reflection B<->W preserves feasibility. The only possible states are

- X: entry residue 0 and head B or W;
- Y: entry residue 1 and head B or W;
- Z: entry residue 1 and head G.

When a letter occurs in a state set, both B/W choices are feasible in X or Y. Let F={X,Y,Z} and Q={X,Z}. Empty head-state sets are permitted. These statements concern existence, so larger feasible sets can always be restricted to F or Q for induction.

For a single nonterminal pair with one paired child and no unpaired child, adjoining a pair outside a subtree maps its possible-state mask by

| S | empty | X | Y | XY | Z | XZ | YZ | XYZ |
|---|---|---|---|---|---|---|---|---|
| D(S) | empty | Y | XZ | XYZ | XZ | XYZ | XZ | XYZ |

This is checked directly: a B parent allows a B or G child, W allows W or G, G allows B,W,or G (one incoming plus one outgoing G). A B/W adds 1 modulo 2; G adds 0 and is permitted only at residue 1. Thus D(Y)=Q, D(Q)=F, D(F)=F.

At a helix terminus, motif avoidance and maximality give exactly: E, no child; L, unpaired child(ren), with zero or one child helix; M, no unpaired child and two or three child helices. Assume each outgoing helix supports F except possibly one supporting Q. The exact terminal-pair state sets are:

| terminal | E | L | M |
|---|---|---|---|
| state set before adding outer stack pairs | F | {Y} | Q |

For E, only the gray-residue constraint applies, so F. For L the loop residue must be 0 and its closing color non-gray, giving Y; if a paired child is present, it can start B/W at residue 0, using the same color as the closing pair. Both F and Q permit that.

For M, the loop residue must be 1: the at least three exposed incidences cannot all be B/W, so one is gray. At residue 1, a closing B allows child colors B,G (two children) or B,G,G (three), and similarly for W. A closing G permits B,W when there are two children or B,W,G when there are three. If one designated child only supports Q, assign it G; for a G closing pair and two children, use G and B (or W), which is also proper. This exhausts all closing colors at residue 1, exactly Q. At residue 0 no gray incidence is possible and there are too many incidences for B/W.

A length-2 helix adds one application of D, so every case supports Q. A length>=3 helix adds at least two applications: E gives F, L gives D²(Y)=F, M gives D²(Q)=F, and further additions keep F fixed. This proves the lemma for every length, rather than testing a bounded list of lengths.

At the virtual root, zero or one Q child and the rest F can always be joined. If unpaired positions exist, there are at most two children; choose relative root residue 0 and assign children B,W. If there are no unpaired positions and at most two children, use the same choice. For three or four children choose root residue 1: assign B,W,G or B,W,G,G, assigning G to the Q child if there is one. This respects both capacity and all entry-state sets.

## Corollary 1: arbitrary comparable isolated stacks

Let T avoid m5 and m3dot, have no length-1 maximal helix, and have every length-2 helix comparable to every other by ancestry in its contracted helix tree. Then T admits a proper uniform modulo-2 separated coloring. Consequently T is uniquely designable in the strict Watson-Crick maximum-pair model.

Proof. In any helix subtree, all short helices lie on one ancestor chain. Hence at most one outgoing child subtree contains a short helix. By bottom-up induction, a subtree without a short helix supports F (its own helix is long, and all children support F). A subtree with a short helix supports at least Q: its children are all F except at most one Q, and the common local lemma applies for either its own short or long helix. At the virtual root at most one child subtree contains short helices, so the root lemma joins the subtrees. Each induction step preserves the residue and properness conditions and leaves the target unchanged. There is no bound on the number of short helices. The implication from separated coloring to a unique maximum-pair sequence is the prior Haleš theorem, not a new claim. QED.

## Theorem 2: complete three-helix classification

Let T be motif-free and have exactly three maximal helices, all of length 2. Unpaired runs may have any finite nonnegative length consistent with the three maximal helices. T fails modulo-2 separation if and only if it has one of these two forms:

(A) Three root-child helices, each ending in a nonempty unpaired hairpin. In dot-bracket notation: `((.+))((.+))((.+))`, where `.+` means one or more dots.

(B) One root-child length-2 helix ends in a two-child multiloop. Both children are length-2 helices with nonempty unpaired hairpins. The exterior has at least one unpaired nucleotide before or after the outer helix. In notation: `.*((((.+))((.+)))).*`, with the two exterior dot counts having positive sum.

All other targets in this class are uniquely designable by a uniform modulo-2 separated coloring. The theorem does not say that failures of this certificate are undesignable.

Proof. The ordered contracted helix forest has exactly three nodes. With three root children it consists of three leaves. A leaf helix has exact state F for an empty hairpin, Q for a nonempty one. If not all three are Q, choose relative root residue 1. With two Q children, assign them G,G and the F child B. With one Q child, assign it G and the F children B,W. With zero Q children, use B,W,G. Each choice is proper. If all three are Q, at root residue 0 all must start B/W, exceeding its two non-gray incidences; at residue 1 all must start G, exceeding its two gray incidences. Thus exactly (A) fails.

With two root children, one child subtree contains two helix nodes and the other one. The two-node subtree is a chain: its intervening loop is L, since without unpaired positions the helices would merge. The local lemma gives its state Q. The root with up to two Q children can use residue 0 and B,W, even if exterior unpaired positions exist. So all such targets succeed.

With one root child, the contracted tree is a three-node chain or a fork. The chain succeeds by Corollary 1. For a fork the two leaf child helices have state F or Q. If at most one is Q, the common local lemma shows the outer length-2 helix supports F, and every exterior pattern succeeds. If both are Q, the multiloop can only occur at residue 1 and both outgoing heads are G. The closing color must therefore be B/W. Consequently the exact terminal state is X, and the length-2 outer helix has state D(X)=Y. A root with no exterior unpaired positions may choose relative residue 1 and join this child. A root with exterior unpaired positions must have relative residue 0 and cannot join Y. This is exactly (B). These are all ordered forests on three nodes. QED.

For both failure geometries a gray node is forced by a degree-three set of incidences, and unpaired nodes are present. Thus any modulo-2 separated coloring would necessarily have the uniform two-residue form tested above; vacuous separation cannot evade either argument.

## Finite domain and independent computation

A run of unpaired positions can be compressed to one position without affecting any coloring constraint: they have the same parent and level. Therefore gap-presence normalization is exhaustive for the two coloring theorems, and arbitrary positive dot-run lengths are covered by the proofs. It does not identify general thermodynamic folding models.

`verify_designs.py` generates every three-node ordered helix forest. Zero-child loops have 0/1 dot; one-child terminal loops have a nonempty choice among their two gaps; two/three-child loops have no dots; root gaps are arbitrary when it has <=2 children and all empty otherwise. Each helix is expanded to exactly two pairs. The resulting 288 distinct strings each have exactly six pairs and exactly three length-2 maximal helices. Every one of their 729 possible B/W/G colorings is checked directly on the interval tree, without importing the compressed-state checker. Exactly four normalized strings fail, corresponding to (A) once and (B) with exterior unpaired positions on the left, right, or both. Each other string has an explicit sequence whose exact Nussinov optimum is 6 with count 1 (284 sequences). No floating-point computations occur.

A separately chosen nested family includes 3,4,5,7,9,12 short helices, long side branches of alternating lengths 3/4, unpaired vertices at several depths, and an exterior unpaired vertex. Direct interval-tree recursion and exact Nussinov checks give unique target optima up to length 170. This is supplementary evidence for Corollary 1, not the proof of its infinite quantifier.

## Scope and gap list

A separate AI review found no mathematical gap in this precise scope. External human expert review and priority confirmation remain pending. The full source problem with arbitrary additional long helices and non-nested short placements remains unresolved. The algorithmic decision problem was already solved by Boury et al.; the contribution, if priority survives review, is only the explicit structural classification and unlimited nested-stack success theorem. The known three-root-branch example is credited to prior work. No thermodynamic, kinetic, wet-lab, or biological-realizability conclusion follows.
