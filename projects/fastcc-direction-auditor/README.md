# FASTCC direction auditor

A negative-only irreversible reaction can cause the existing COBRApy FASTCC implementation to build an infeasible auxiliary LP for an otherwise feasible metabolic model. This project supplies an exact small-network oracle, a production-only patch, and a portable reproduction at upstream commit `a9c5ced29670bed5035fb19e87e7b0b588adb62d`.

For a source–transfer–drain pathway represented with bounds `[-10,-1]`, mass balance permits flux `(-7,-7,-7)`. FASTCC nevertheless imposes `v >= z >= 0` on its negative-only reactions. The patch uses `-v` as the nonnegative detection quantity. It does not change the original stoichiometry, bounds, objective or preexisting solver constraints.

The independent review suite gives **1 failed / 4 passed before the patch, and 5 passed after it**. The exact rational verifier covers 24 direction/bound combinations. Tests cover positive and strict-negative paths, blocked and fixed-zero controls, a reversible cycle, and preservation of nested contexts and solver state. Archived runs used GLPK; see [scope and limitations](SCOPE.md).

```sh
git clone https://github.com/opencobra/cobrapy.git artifacts/fastcc-upstream
git -C artifacts/fastcc-upstream checkout a9c5ced29670bed5035fb19e87e7b0b588adb62d
python3.11 -m venv .venv
. .venv/bin/activate
python -m pip install -r projects/fastcc-direction-auditor/requirements.txt
python -B projects/fastcc-direction-auditor/reproduce.py \
  --checkout artifacts/fastcc-upstream
python -B projects/fastcc-direction-auditor/exact_oracle.py
```

Run these commands from the lab repository root. Python 3.11 and Git are required. Dependency installation requires network access; the reproduction itself is local. Both scripts print JSON by default and leave archived evidence unchanged. To retain a separate run:

```sh
python -B projects/fastcc-direction-auditor/reproduce.py \
  --checkout artifacts/fastcc-upstream \
  --output artifacts/fastcc-direction/receipt.json \
  --logs artifacts/fastcc-direction/logs
```

The runner checks 80 frozen source-file hashes and the actual COBRApy/FASTCC imports, copies the complete source into temporary baseline/candidate directories, and applies only [upstream.patch](upstream.patch) to the candidate. Both pytest processes use an explicit seven-key environment allowlist and `--tb=short`. Requested logs replace local source, project, temporary and environment paths with labels. The input checkout is never patched. Outputs must be outside this project directory.

This is an implementation repair, not a new biological mechanism, new theorem, or proof of full FASTCC correctness. The original algorithm assumes nonnegative irreversible representations; the change supplies the corresponding normalization inside the existing API. [PR1495](https://github.com/opencobra/cobrapy/pull/1495) already covers a related negative-only pathway in the new SPECTRA sibling and discusses FASTCC false negatives; it does not patch the existing FASTCC source. Credit and inspected versions are recorded in [SOURCE.md](SOURCE.md).

The existing-API repair is submitted as [COBRApy PR1500](https://github.com/opencobra/cobrapy/pull/1500), with maintainer acceptance pending. The repair received a separate AI review and independent checks. These do not constitute external human peer review or maintainer acceptance. See [the review](evidence/independent-review.md), [reproduction receipt](evidence/reproduction.json), [archive hashes](archive_manifest.json), and [third-party licensing](THIRD_PARTY.md).
