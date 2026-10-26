# Exact audit of a multipathogen coexistence conjecture

Two admissible parameter examples contradict the simultaneous-endemic-equilibrium assertion in Conjecture 5 of Djuikem–Arino (2025). Each has two mutually transmitting hosts, one subcritical pathogen and two supercritical pathogens, but the specified joint endemic equilibrium cannot exist. A second example uses full-rank positive contact matrices. The biological explanation uses classical competition for a shared susceptible pool.

This is an AI-assisted, paper-specific mathematical correction awaiting author/human review. It is not a new experimental discovery or a general classification of pathogen coexistence. See source.md, proof.md and provenance/review-receipt.json for exact scope and review limits.

## Reproduce

Use Python 3.12 and install the exact dependencies before verification:

```sh
python -m pip install -r requirements.txt
python -B verify.py verify
python -B -O verify.py verify
```

The verification command uses no network, API or GPU. Dependency installation requires access to the package index. Verification reads the exact manifest inventory, runs two separate exact checkers in external temporary directories, compares every certificate byte, and rechecks the archive. It does not write caches or outputs into the project. Both workers run with optimization0 even when the outer wrapper uses `-O`, retaining their internal scientific checks.

Regenerate both certificates and the full manifest into a new absolute external directory:

```sh
python -B verify.py regenerate --output-dir /absolute/new/external-directory
```

The output directory must not already exist or lie inside the project. Regeneration produces a complete replayable copy, with recomputed certificates and manifest. It never overwrites the original archive.

The hash manifest detects inventory and byte inconsistencies; it is not a cryptographic author signature or proof of mathematical novelty. No third-party paper, private correspondence or unpublished source material is included. MIT licensing applies to the original code and exposition in this packet; the cited article remains its authors' work.
