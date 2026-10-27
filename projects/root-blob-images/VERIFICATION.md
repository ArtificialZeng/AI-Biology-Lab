# Verification coverage

`check_exact.py` uses exact `fractions.Fraction` arithmetic. Its certificate
contains all 256 rational leaf-pattern probabilities for the collapsed target,
their sum and strict positivity, three observable internal-edge reconstructions,
the explicit character table and edge spectrum, and a rational split of the
target's root pendant. Twenty-five support-set pairs check the finite
combinatorial support lemma. A further 735 complete finite support combinations
check the broader two-attachment obstruction, including kernels with unit
nontrivial eigenvalues. Nine compensating constructions cover JC, K2P and
K3P targets at three sets of fixed inheritance weights; each is compared with
direct summation over the actual ancestral random variables.

`check_graph.py` independently checks the explicit diamond's graph contract,
cut edges, attachment vertices and collapsed quartet. It exhausts simple root
insertion on eleven edges and all remaining orientations, finding five valid
root positions and twenty parent-choice paths. This finite domain excludes
arbitrary hidden expansions; graph-proof.md states the separate cut argument.

`verify.py` binds every regular project file except the manifest itself,
rejects extra files and symlinks, runs both checkers in temporary working
directories, compares their complete archived outputs and rechecks the input
hashes. It does not regenerate evidence or install dependencies. Both normal
and optimized wrapper invocations preserve the archive. Scientific subprocesses
explicitly set `PYTHONOPTIMIZE=0` and disable bytecode writing. The scientific
checker also uses explicit exceptions rather than removable `assert` statements.

The infinite graph theorem, sink conditions and all-parameter nonmembership are
analytical claims in proof.md. No finite certificate proves them alone. The
archive is not a trusted-runner release record, formal proof, human review or
biological validation. Changes to mathematical inputs require renewed review;
updating hashes alone does not verify a theorem.
