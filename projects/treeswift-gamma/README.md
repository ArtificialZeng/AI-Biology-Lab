# TreeSwift gamma calculation audit

TreeSwift v1.1.51 computes an incorrect Pybus–Harvey diversification summary by using absolute branching depths as waiting intervals. The local correction constructs intervals from all branching events, independently of node labels, ignores incoming crown stems and unary prefixes, and retains arbitrary zero-edge resolution of polytomies.

For the positive-length ultrametric tree `((A:1,B:1)I:9,C:10)R;`, the three tips are at time 10 and the internal split is at time 9. Thus g2=9, g3=1 and total lineage time T=21. The established formula gives **5√3/7**, whereas the original implementation returns **−√3**. This can reverse the sign of the diversification summary. We do not claim that any observed study's conclusion changed.

This is an engineering repair, not a new biological result, statistical method or mathematical discovery. The statistic assumes an extant-tip ultrametric tree. Fewer than three tips retain the existing RuntimeError; zero total lineage time remains undefined. No broad validation rejecting other previously callable input classes is added. Finite checks do not establish universal floating-point accuracy.

The archived original and modified runtimes derive from author commit `e8ff53deabf57fa0dc8e8fba3daeaf90b96e413c`, identical to latest tag v1.1.51 at the audit. Only `gamma_statistic` and its docstring change, with one added stdlib unittest module. The full original three-module Python runtime and existing test fixtures are retained byte-for-byte; the minimal patch is included.

Run with Python 3.12 (standard library only):

```
python3 -B verify.py verify
python3 -O -B verify.py verify
```

Verification is read-only, rejects missing/modified/extra files and symlinks, and reruns scientific workers in an external disposable copy with assertions enabled (`PYTHONOPTIMIZE=0`). Existing upstream tests write only inside an external temporary working directory. Both invocations produce identical deterministic stdout.

Explicit regeneration writes the two certificates and manifest to a new, absolute directory outside this archive; an existing directory is rejected:

```
python3 -B verify.py regenerate --output-dir /absolute/new-directory
```

The original Python Fraction/event-count oracle and 714 deterministic fixtures cover generated 3–8-tip binary histories, labels, scales, root stems, unary nodes, zero-gamma controls, tied times, zero edges and true polytomies. Candidate: 714 matches, no exceptions/mutations. Baseline nonoverlapping outcomes: 24 matches, 622 numerical mismatches, 68 exceptions. Actual floating outputs are compared at absolute/relative tolerance 2e-13; exact rational coefficients, radicals and squared expectations are retained.

A separate fresh checker builds 1774 adversarial event histories, including very small/large scales and polytomies, without using the owner or candidate timing helpers. It retains its reviewed absolute/relative tolerance 2e-12; all 1774 candidate cases match without mutations. Its only portable source adaptation is the runtime directory binding. These finite generated checks are not exhaustive over all trees or edge assignments, and fresh AI review is not external human peer review.

Formula credit: [Pybus and Harvey (2000)](https://doi.org/10.1098/rspb.2000.1278), and Emmanuel Paradis's independently authored [ape::gammaStat](https://raw.githubusercontent.com/cran/ape/master/R/gammaStat.R). The archived independent Python implementations are original; no R source or paper PDF is redistributed. Original scanned paper PDF access was unavailable during this audit; no inspection of that PDF is claimed.

Repository and public contract: [TreeSwift](https://github.com/niemasd/TreeSwift), [author documentation](https://niema.net/TreeSwift/), [v1.1.51](https://github.com/niemasd/TreeSwift/releases/tag/v1.1.51).

## License and assistance

**This entire subproject is GPL-3.0 and overrides any enclosing repository's MIT default.** Upstream full GPL-3.0 LICENSE and attribution are preserved. TreeSwift author: Niema Moshiri; citation: Moshiri (2020), SoftwareX 11:100436, [doi:10.1016/j.softx.2020.100436](https://doi.org/10.1016/j.softx.2020.100436). AI assistance was used to audit, repair, generate cases and perform a separate review; evidence and limitations are provided for independent inspection.
