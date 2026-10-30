# Review and validation scope

The source research used substantial AI assistance. A separate fresh AI context reconstructed the LWL85 aggregation, checked both public APIs and both sequence orders, added a held-out alignment with L2<L4 and all transition/transversion classes nonzero, and independently recalculated 32 existing API expectations. This is an internal scientific-software check, not human peer review or maintainer acceptance.

The six-file patch changes the two numerator weights in each API, adds three manual-count regression methods, and updates dependent existing LWL85 expectations and one tutorial output. The initial full offline suite exposed the tutorial discrepancy; an unpatched-copy comparison passed, confirming that the expected output needed to change with the fix. That failure and the later repair were retained rather than described as an unrelated pre-existing failure.

All specified pre-commit hooks passed locally for the six changed files, without changing their frozen contents. Three related test modules and the 89-example tutorial doctest were also checked. The full offline-suite outcome and upstream CI are reported separately in the pull request; this packet's small release audit does not claim to execute the complete upstream suite. Network tests and unavailable optional dependencies are outside an offline run.

The public certificate uses four deliberately bounded fixtures. Its rational site counts are exact; evaluated logarithms and rate comparisons are floating point. Broad biological impact, historical priority, universal codon conventions and other dN/dS methods remain unassessed. No third-party reply or private contact record is included.
