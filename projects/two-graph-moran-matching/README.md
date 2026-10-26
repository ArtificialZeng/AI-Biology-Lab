# Complete-graph fixation bounds: regular graphs and matching deletions

This project studies neutral **type-specific Moran birth–death**: a parent is chosen uniformly; its type determines the graph on which its offspring replaces a uniformly chosen neighbor. All graphs here are simple, connected, undirected and unweighted.

For every noncomplete d-regular graph G on N vertices and every starting vertex v, put q=N−1−d>0. The proposed exact bounds are

    h_v(G,K_N) ≤ 1/N − q/(d N³) < 1/N,
    h_v(K_N,G) ≥ 1/N + 2q/[d N⁴(N−1)] > 1/N.

The first argument gives the mutant graph and the second gives the resident graph. No vertex-transitivity or asymptotic limit is assumed. These bounds come from an exact squared-drift identity on every labelled configuration, finite absorption and one- or two-step comparison. They do not assert that all regular graphs have the same fixation probability or give sharp fixation formulas. See the [regular-graph proof](regular-graphs/proof.md), [scope](regular-graphs/problem.md), [primary sources](regular-graphs/source.md), and [bounded overlap review](regular-graphs/literature.md).

The earlier matching-deletion theorem remains a complementary subclass. For all N≥3 and 1≤m≤⌊N/2⌋, every singleton mutant on K_N minus m disjoint edges, against complete-graph residents, satisfies

    h_v ≤ 1/N − 2m/[N²(N−1)(N(N−2)+2m−2)] < 1/N.

See [all-matching proof](all-matchings/proof.md). Except for perfect matchings, these complements are generally irregular, so the regular-graph theorem does not replace the full matching theorem. The historical m=2 [proof](proof.md) and all earlier scientific workers and certificates are preserved byte-for-byte. Historical wording that calls larger matching cases unresolved is superseded by the all-matching proof. This is a single augmentation of existing lab PR18.

The broader conjectures for arbitrary irregular graphs remain unresolved. The biological model, known single-edge exact formulas and numerical regular-graph examples are due to Tkadlec, Kaveh, Chatterjee and Nowak, [Evolutionary dynamics of mutants that modify population structure](https://doi.org/10.1098/rsif.2023.0355) (2023). Their [official Figshare deposit](https://doi.org/10.6084/m9.figshare.16910170) is credited. The earlier type-specific model is due to Melissourgos et al., [An extension of the Moran process using type-specific connection graphs](https://doi.org/10.1016/j.jcss.2021.07.007). No appropriate author GitHub submission endpoint was verified, and no author paper or software is redistributed.

The derivation and checks used AI. A separate fresh AI context independently reconstructs the argument, checks held-out labelled graphs and solves selected finite fixation systems. This is computational and mathematical AI review, not human peer review or proof-assistant certification. Trusted runner provenance was unavailable; the formal workflow gate did not certify the work. Neither passing checks nor a bounded literature search establishes universal priority. Frozen proof documents retain their pre-review status wording; current review scope and reproduction mappings are recorded under `provenance/`.

## Reproduce

Use Python3.12 and pinned SymPy1.13.3/mpmath1.3.0:

    python -m pip install -r requirements.txt
    python -B verify.py verify
    python -B -O verify.py verify

Both wrapper modes run eight scientific workers with PYTHONOPTIMIZE=0 in external disposable directories, retaining assertions. Verification checks the complete manifest before and after recomputing all eight exact certificates; it rejects altered, missing, unexpected or symlink inputs. Project files are never output targets. The finite graph fixtures are explicitly fixed checks, not exhaustive classification of all regular graphs. The infinite statements depend on the written proof.

To regenerate all certificates and the entire manifest, choose a new absolute directory outside the project, whose parent already exists:

    python -B verify.py regenerate --output-dir /absolute/new-directory

The wrapper refuses existing and in-project destinations. The resulting certificates and manifest must match the archive byte-for-byte. The manifest is a consistency inventory, not a signature or mathematical proof. Direct worker execution in the archive is unsupported; use the wrapper for output isolation. No network, GPU or author software is needed after installing dependencies. Original code is MIT licensed. Linux/macOS CI replays the read-only normal and optimized wrappers.
