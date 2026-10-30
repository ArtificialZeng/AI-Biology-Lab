# Explicit root-trapping diamond

Use leaves 1,2,3,4 and rooted DAG edges

```
rho -> u, v
u   -> a, b
v   -> a, b
a   -> s
b   -> t
s   -> 1, 2
t   -> 3, 4
```

Vertices a,b are the two reticulations. The other nonleaf vertices except rho are tree vertices. Every edge is unique; the displayed order is topological. The two root children have routes to both reticulations, so neither u nor v lies on every root-to-leaf path. Vertex a can be avoided on paths to 3,4, and b can be avoided on paths to 1,2; s,t and each leaf are likewise avoidable. Thus rho alone lies on all root-to-leaf paths, as required by source Definition 2.1.

Undirect nonreticulation edges and suppress rho. The resulting semi-directed core has undirected edge u--v and directed edges u->a, v->a, u->b, v->b. Set B={u,v,a,b}. The only neighbors of B outside B are s adjacent to a and t adjacent to b. The induced core is connected and consists of two triangles sharing u--v; every one of its five edges lies in a cycle. None is a cut edge of the whole graph. B therefore is a 2-sub-blob and, in fact, the unique nontrivial maximal blob.

The cut edges a--s and b--t separate the two outside components. A root located in the s,1,2 component would need a directed path into B across the unique attachment a--s, directed toward a. But a already has the two fixed incoming reticulation edges from u,v and exactly one outgoing edge. Its outside attachment is necessarily outgoing, so no such entrance is possible. The same argument excludes a root in the t,3,4 component or on either attachment edge. Hence the root is trapped inside B, including the source convention that a hidden larger blob may be suppressed into B. In an expanded rooting, a root outside the preimage of B would still have to cross one of these outward attachment cuts, giving the same obstruction. This direction/cut argument is the topological reason; the finite orientation enumeration is a separate check on simple edge-insertion rootings, not an enumeration of arbitrary hidden expansions.

Contracting B produces a degree-two vertex between s and t. Suppressing that vertex leaves the quartet tree with internal edge s--t and split 12|34. No outside reticulation remains.

For the displayed-tree choices in the displayed rooted DAG, if a,b both retain their edge from u, their unique attachment path is a-u-b; if both choose v it is a-v-b. If they choose distinct parents, the paths are a-u-rho-v-b and a-v-rho-u-b. Thus the respective original-edge counts are 2,2,4,4. All four choices have strictly positive mixing weight when both mixing parameters lie in (0,1).

More generally, the checker exhausts simple root insertion on all eleven edges and all orientations of the remaining undirected edges, with fixed incoming reticulation directions, binary degrees, acyclicity, reachability and last-stable-ancestor conditions. Exactly five root positions are valid, all on B's edges; the six outside positions are invalid. For each valid rooting it checks all four reticulation choices, giving twenty attachment paths, all of length at least two (and at most four). These paths count original kernels before suppressing vertices. Since the semi-directed graph has no edge a--b, re-rooting cannot create a one-edge original attachment path.

This establishes a lawful graph for the proposed boundary test. It does not establish any transition-kernel claim, image nonmembership, failure of splittability or solution of the published equality question.
