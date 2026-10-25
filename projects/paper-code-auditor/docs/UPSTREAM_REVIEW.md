# Review packet: execution-route corrections

Submitted with the lab owner's authorization on 2026-10-26: [upstream PR #1](https://github.com/ziul-bio/DMS-SiteEffect-PLM/pull/1). Head commit: `9e924cc7ad5fd52a1e8dca02470b30448e58a665`. Open at submission; maintainer acceptance and training validation remain pending.

## Version and scope

- Repository: https://github.com/ziul-bio/DMS-SiteEffect-PLM
- Current version checked: `e168a2c825be0c3e28c300b8f6277dda1d999c4d`
- Public historical recovery source: `e0c2d6ae06903cb356abb7e1b579058115e5623a`
- Purpose: repair specific public reproduction entrypoints without substituting a different model or training method.

## Evidence

The README-recommended cellular wrapper defines one dataset array but uses a different, undefined array in its inner loop. The isolated control flow schedules zero calls and returns success. Correcting the array schedules the intended 198 calls in the counter test; this does not launch training.

Both wrappers reference a partial trainer and matching model module deleted from the public tree. Their public historical versions can be recovered together. The restoration retains the partial model source. Module invocation supplies the documented repository-root import context. The original wrapper/trainer also contain an inconsistent one-visible-GPU/logical-index-1 configuration; the patch uses logical index 0 and allows an existing visibility setting to override the default.

Separately, the current LoRA data module passes two keywords that its split helper does not accept. The proposed one-line correction uses its existing `test_size` argument. The split helper itself is unchanged.

## Validation performed

`python projects/paper-code-auditor/audit.py` runs offline. It checks source integrity, actual extracted data-setup behavior, three seeded split comparisons, patch applicability, reconstructed loop counts, restored source identity, shell syntax, and checkpoint-argument compatibility. The input is a minimal identifier-only projection of an upstream assay, not synthetic biological observations.

## What remains to establish

The original authors should assess whether the historical partial route is the intended public reproduction route and whether additional unpublished configuration is required. Compatible GPU dependencies and end-to-end training have not been tested. No conclusion about the validity of the paper's reported metrics follows from these software checks alone.

The original Lasso and LoRA helper splits agree in the checked cases. We explicitly withdrew a preliminary different-fold hypothesis. This packet should not be represented as evidence of an unfair-comparison result.
