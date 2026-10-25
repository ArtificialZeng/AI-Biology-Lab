# Reproduction

Use Python 3.12 or 3.13. Commands below run from the repository root. No model API, GPU, or paid service is required.

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-exact.txt
python scripts/reproduce.py verify
```

The exact route checks the archive's mathematical-input hashes, copies the project to a temporary directory, and runs:

1. `manuscript/verify_certificate.py`: reconstruct the manuscript's six-by-six matrix, Hurwitz determinants, 125 coefficients, and table bounds.
2. `src/verify_finite_kinase.py`: reconstruct the nine-dimensional characteristic polynomial and independently check Bernstein expansions with rational evaluations.
3. `audit/verify_reduction.py`: rederive conservation laws and equilibrium identities, verify the tangent-space reduction and exact Routh root count, and check all Bernstein polynomial identities.

Outputs and logs go to `artifacts/mixed-phosphorylation/verify/`. A successful `result.json` records `status: passed`. Each command exits on a failed assertion. Do not use Python's `-O` flag or `PYTHONOPTIMIZE`: these scripts rely on assertions.

## Regenerate the certificate

```sh
python scripts/reproduce.py regenerate
```

This independently rebuilds the algebra and Bernstein files, compares their full JSON content with the archived inputs, and then runs all three exact verifiers. Generated outputs remain under ignored `artifacts/`. Regeneration is more expensive than verification.

## Reproduce numerical observations

```sh
python -m pip install -r requirements.txt
python scripts/reproduce.py numerics
```

This runs `baseline.py 100 1000 10000 inf`, the independently implemented Radau/BDF and periodic-shooting checks, first-order response computations, and figure generation. It saves results under `artifacts/mixed-phosphorylation/numerics/`. Floating-point outputs can vary across platforms; the recorded shooting periods are approximately 4.77828075, 4.69985145, 4.69286392, and 4.69209600 for K=100, 1000, 10000, and the limit, respectively. These are numerical reference values, not interval enclosures.

To redraw the figure from archived numerical inputs without rerunning the full integrations:

```sh
python scripts/reproduce.py figure
```

## Build the manuscript

A TeX distribution providing pdfLaTeX and latexmk is needed only for rebuilding the PDF. From the project directory:

```sh
cd projects/mixed-phosphorylation
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build manuscript/main.tex
```

The archived `manuscript/main.pdf` is the original five-page note. The mathematical model, coefficients, and assumptions are described in it and in the proof. Original third-party articles are linked in the bibliography rather than bundled.

## New interval orbit and period correction

The equilibrium-only Python command above does not rerun the new periodic-orbit integration. For that claim, use the [separate validated-cycle instructions](../projects/mixed-phosphorylation/validated-cycle/README.md): pinned CAPD, Linux x86/FILIB, then the parameter-sensitivity checker. Independent exact postprocessors verify the archived widened interval bounds; they complement the fresh integration. No explicit finite-K threshold or remainder constant is claimed.
