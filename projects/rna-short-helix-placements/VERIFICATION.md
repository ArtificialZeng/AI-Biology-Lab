# Exact verification and archive replay

Run from any directory using the path to `replay.py`:

```sh
python projects/rna-short-helix-placements/replay.py
```

Python 3.12/3.13 and its standard library suffice. No network requests, package installation, numerical library, GPU, model weights, or API key are used. The CI matrix requests Ubuntu and macOS with both Python versions. A local run does not certify that all hosted CI jobs have completed.

## What is checked

1. `archive_manifest.json` binds the proof, scope, source/credit documents, checker programs, replay program, and deterministic JSON archives by SHA-256. Missing files, changed bytes, or malformed manifest entries cause an explicit error.
2. `verify_states.py` enumerates the finite entry-state transition algebra and permitted terminal-loop joins. Its F/Q invariants support the proof's induction. Testing helix lengths 2–8 is supplementary: the identity D(F)=F is what extends the stated argument to all longer helices.
3. `verify_designs.py` builds the 288 normalized three-helix target classes and checks each of their 729 B/W/G assignments directly on the interval tree. It constructs one sequence for each of the 284 positive targets. An exact integer Nussinov recurrence decomposes at the leftmost nucleotide and returns optimum/count `(6,1)` for every such sequence. It also checks the six larger examples and their actual short-helix comparability.
4. `independent_verify.py` independently regenerates the domain from every six-pair Dyck skeleton with motif-legal binary dot gaps, checks ordinary and uniform separation separately, verifies the serialized sequence/color witnesses, and uses a rightmost-nucleotide optimum/count recurrence. Its negative controls expose ties and sequences without pairs. It imports neither of the other checker programs.
5. Replay runs each checker in both ordinary and optimized Python. Each complete parsed result must equal its archived JSON; the raw stdout must also agree between modes. The three checker programs do not import each other. The subprocesses use isolated mode and disable bytecode writes. Failure checks use explicit exceptions, not `assert` statements removed by `-O`.

The replay output records hashes and counts, but passing it is a computational reproduction, not a proof of priority, human peer review, or formal theorem-prover certification. The full source placement question remains unresolved.

## Archived receipts

- `evidence/states.json`: complete finite-state table and local checks.
- `evidence/designs.json`: all 288 normalized targets, coloring counts, 284 constructive sequences and optimum counts, plus six larger examples.
- `evidence/independent.json`: independent generation, both separation counts, rightmost-fold checks, local algebra, and negative controls.

The normalized domain has five ordered helix-forest shapes, contributing `8 + 96 + 96 + 72 + 16 = 288` classes. The three-root-leaf shape contributes one modulo-2 failure; the fork contributes three, corresponding to exterior dots before, after, or on both sides. The other shapes contribute none. Arbitrary positive dot-run lengths are covered by the proof's equal-parent/equal-level argument, not by enumerating an infinite list of sequences.

Ordinary modulo-2 separation can hold vacuously when gray or unpaired nodes are absent; the archived checker counts the stronger uniform form. The two failure geometries force both kinds of nodes to exist, so the proof establishes failure of ordinary modulo-2 separation there as well. The positive uniform witnesses are sufficient everywhere else. Uniform and ordinary coloring **counts** must not be equated.

## Archive policy

The default commands read the archive and print to stdout. To save a run, redirect stdout to an ignored `artifacts/` directory of your choice. The checkers never overwrite the archived receipts, and replay provides no automatic manifest-repair option. A changed mathematical input requires a justified proof/certificate update and fresh review before replacing the archive. Merely updating hashes does not verify a theorem.

A separate AI reviewer reconstructed the local algebra, regenerated the finite domain from six-pair Dyck words, and used a rightmost-nucleotide optimum/count recurrence to check the supplied sequences. That review found no mathematical gap and requested the attribution and status qualifications retained in the public documents. It is not external human review, Lean verification, or release-gate acceptance.
