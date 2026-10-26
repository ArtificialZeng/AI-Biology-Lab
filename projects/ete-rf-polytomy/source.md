# Primary sources and scope

The exact source is [etetoolkit/ete](https://github.com/etetoolkit/ete) default
branch ete4, commit9562dfb6a02795dfda3975be5025f83b8dc884b1.
The source archive and production file checksums are in provenance/source-mapping.json.
[Issue745](https://github.com/etetoolkit/ete/issues/745) is an existing open
report about unrooted root-degree correction; this packet supplies an additive
repair and negative-score/common-leaf regressions.

[ETE's tutorial](https://etetoolkit.github.io/ete/tutorial/tutorial_trees.html)
describes conceptual unrooted trees, master-root representations and comparisons
on shared labels. [The current method reference](https://etetoolkit.github.io/ete/reference/reference_tree.html)
describes the RF return contract. The inspected source's production correction
uses full-tree raw child counts while its RF partition sets use shared labels.
The existing support-filter tests in tests/test_tree.py are commented out as TODO.

CONTRIBUTING.md requests tests and ./run_tests.py. Both native default fast runs
and their common failures are disclosed in README and evidence/validation-summary.json.
The inspected complete current tree and contribution guide had no AI-specific
restriction. No maintenance acceptance or exhaustive novelty search is implied.
