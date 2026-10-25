# Working in AI Biology Lab

## Purpose

Pursue original, externally checkable results in mathematical and computational biology. The lab aims for discoveries with a clear scientific and public story. Prefer one useful, verified contribution to a large collection of unsupported claims.

## Research work

- Start with a precise question and primary sources. Reproduce the published model before changing it.
- Distinguish prior results, exact theorems, conditional arguments, floating-point evidence, and future targets.
- Check all parameter assumptions and quantifiers. Equilibrium instability does not imply periodic-orbit existence in these dimensions.
- Review directly overlapping work before claiming novelty. Do not infer an unresolved problem from an old paper alone.
- Seek a proof and a counterexample symmetrically. Record failure or uncertainty honestly.
- AI review, passing scripts, and GitHub Actions are computational checks, not external human peer review.
- Keep new projects self-contained under `projects/`, each with a statement, sources, reproduction route, evidence, and status.

## Reproduction and publication

- Use `python scripts/reproduce.py verify` for the current exact result. Use `regenerate` after changes to certificate generation.
- Numerical reproduction uses `python scripts/reproduce.py numerics`; figures can be redrawn using `figure`.
- Keep archived evidence unchanged during ordinary checks. Outputs belong under ignored `artifacts/`.
- If mathematical inputs change, justify the change and update certificates, proofs, and `archive_manifest.json` together. Updating hashes alone does not verify a changed theorem.
- Public documentation must retain limitations and credit the source model.
- Do not publish credentials, internal prompts, private correspondence, patient information, or third-party papers without permission.
- External outreach and submissions need authorization from the human user. Repository publication does not authorize contacting researchers or journalists.

## Practical scope

Favor exact algebra, validated low-dimensional computation, and carefully designed benchmarks feasible on a 128GB Apple Silicon Mac. Renting GPUs or paid services requires an established budget. Current reproduction has no API or GPU dependency.
