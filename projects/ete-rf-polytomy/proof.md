# Why the common-leaf resolution count is correct

Let C be the common set of labels accepted by the existing API's duplicate
checks, and n=|C|. Retain the minimal tree spanning C, deleting unretained
terminal branches and suppressing unary rooted paths or unrooted degree-two
vertices. An unrooted master root of degree two is also suppressed. These
operations do not change the represented nontrivial partitions of C.

For rooted trees, let m be the number of distinct proper nontrivial clades.
For n>=2, if I is the number of internal vertices after suppression, m=I-1.
Since the tree has n+I-1 edges, the excess-child sum is

    sum_internal (children(v)-2) = n-I-1 = n-2-m.

Each extra binary resolution decreases this sum by one; it is zero exactly
for a fully binary induced rooted tree. For n<=1 there is no internal resolution
and both the clade count and the nonnegative maximum count are zero.

For unrooted trees with n>=3, every internal vertex has degree at least three.
Writing I for their number, the handshake identity gives

    sum_internal (degree(v)-3) = n-I-2.

The number of nontrivial edges, or distinct nontrivial splits, is m=I-1.
Consequently the missing-resolution count is n-3-m. For n<=2 there are no
nontrivial splits or unresolved internal vertices, so the count is zero.

The replacement helper projects every existing cached leaf clade to C.
For rooted trees it stores proper nontrivial clades once. For unrooted trees
it stores unordered complementary pairs once. Thus duplicates from unary
paths, either half of a degree-two master root, and clades made equal by
removing unshared leaves are counted once. A degree-three unrooted master
root contributes no missing resolution. The helper returns

    max(0,n-2)-m  (rooted), or max(0,n-3)-m  (unrooted).

The max acts on the attainable number of splits for tiny taxon sets, not on
the returned RF score. This is exactly the degree-excess count above, hence
never negative for valid induced trees and independent of representation.
The old exception remains when both induced counts are positive.

When support filtering is disabled, all trivial common partitions cancel and
the raw RF difference is m1+m2-2c, where c is the common nontrivial-part count.
If both deficits are zero the correction is zero. If one deficit d is positive
and the other tree is fully resolved with r possible parts, then m1=r-d,
m2=r and subtraction of d gives 2(m1-c)>=0. This explains why the old negative
quartet score is an accounting error and why the corrected accounting cannot
produce such a negative result in this defined topological setting.

This argument proves the topology count, not correct functioning of existing
support filtering or an optimal-refinement interpretation. The patch leaves
that filtering code intact. Its known failures are separately recorded. Actual
compiled implementation was compared with an independently written graph
pruning/suppression/degree oracle on complete finite domains; this finite audit
is implementation evidence and not an infinite-domain execution proof.
