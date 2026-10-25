# AI use and provenance

Zijian Zeng initiated this project with an AI coding/research assistant. AI-assisted passes contributed model transcription, symbolic derivation, numerical implementation, proof search, literature comparison, exact checks, and manuscript preparation. Separate contexts implemented additional mathematical and numerical checks.

This process is AI-assisted research with inspectable outputs. A separately generated implementation is not external human peer review, and multiple agreeing AI reviews do not establish priority or biological validity.

The first public package was prepared on 2026-10-26 from the local mixed-phosphorylation result. It preserves the mathematical certificate JSON files and the five-page manuscript. Source overviews, reproduction commands, and provenance-dependent checks were adapted for a portable public repository. Internal job identifiers, execution transcripts, system prompts, local paths, build logs, and downloaded third-party papers are not part of the public package.

The second-reduction verifier is adapted from an AI review pass. Its checks reconstruct conservation laws, equilibria, the tangent-space reduction, Hurwitz determinants, and Bernstein identities. The public version hashes its mathematical inputs at the start of each run and checks they remain unchanged during that run; it does not claim to authenticate the private review process.

The root reproducibility script checks the public archive's mathematical-input checksums before running the verifiers in a temporary copy. `projects/mixed-phosphorylation/archive_manifest.json` records the public artifact hashes. This is integrity metadata, not a proof of novelty or external approval.

Future releases should describe substantive AI and human contributions, unsuccessful approaches that affect interpretation, and the evidence supporting each promoted claim. Public model names, cost claims, and autonomy claims should be backed by records before they are reported.

## Validated-cycle and evaluation follow-up

The new period-correction work used AI for exact reduction, construction of interval verifiers, and separate skeptical mathematical audits. CAPD performs the validated integration; independent rational implementations reconstruct equations and postprocess deliberately widened interval output. They do not replace or independently implement CAPD's integrator. Linux CI provides a separate platform run. These checks do not constitute external human peer review or establish priority.

Protein Metric Auditor operates on frozen public numeric scores without model inference or sequence design. MEMOTE direction auditing tests a small exact software invariant and real library entry points; it does not establish genome-scale biological impact. Public descriptions retain these limits and no new third-party messages were sent for these studies.
