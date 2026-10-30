# Measured model failures are not missing observations

This audit reproduces and proposes a correction of the Gaussian hierarchical
likelihood defect already reported by **Fabian Fröhlich (FFroehlich)** in
[pyPESTO issue1794](https://github.com/ICB-DCM/pyPESTO/issues/1794).
It is a scientific-software repair, not a new biological discovery.
Upstream base: `053079aa71710a62c1db17f3d0342e96499a3768`.
The issue and candidate patch still require maintainer review.

A successful ODE integration can produce an undefined measured observable, for
example x/(x+y)=0/0 at time zero. The baseline's NaN-skipping sums discard that
residual while still counting the measurement when estimating noise. In the
small analytic model here, it returns a finite objective −0.9607627105.
The candidate rejects that point and still accepts an undefined output when its
measurement really is absent, as well as genuinely valid zero predictions.

`upstream.patch` contains the exact two-source-file correction and native
regression. Its SHA256 is
`6ca3c00f863e9d6b9cf57b6c88850b64e35b5381e337890914a467052a81a920`.
Four Hessian-containing sensitivity requests already rejected the bad point in
the baseline; this project does not present those as newly fixed failures.

## Fast archive and utility checks

From this project directory, with Python3.12:

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-light.txt
python reproduce.py verify
python reproduce.py light
```

The fast route verifies the manifest, extracts **two unchanged function ASTs**
from the exact archived utility source, and compares them with an independent
closed-form Gaussian oracle. It also checks archived native floating-point
observations by tolerance. It does **not** import the full package, compile
AMICI, or execute the public collector. A light CI job may run these commands;
passing it is not evidence that native integration ran on CI.

## Full native replay

This route fetches the pinned upstream source, so it requires network access.
It builds real AMICI C++ models and requires a working C/C++ compiler and AMICI
build dependencies; installation can take several minutes. Full upstream
sources remain outside this public project.

```sh
python reproduce.py bootstrap --work ../../artifacts/pypesto-nan-replay
python -m pip install -r requirements-native.txt
python -m pip install -e '../../artifacts/pypesto-nan-replay/baseline[petab]'
python reproduce.py native --work ../../artifacts/pypesto-nan-replay
```

Use a new work directory. `bootstrap` verifies the exact commit, applies the
exact patch, and checks both changed source files and the native test against
archived hashes. Optional `--from-local PATH` supplies an existing Git checkout
instead of a network fetch; the commit and content checks remain mandatory.
The replay freezes the chosen package import before pytest can alter paths,
rejects skips/errors, and expects **baseline6failed/33passed** and
**candidate39passed**. Baseline failures are expected reproduction evidence;
the wrapper exits successfully only if those exact scoped outcomes occur.
Results are written in the external work directory, not into this archive.

The original research executed this matrix on macOS Apple Silicon with Python
3.12.13, AMICI1.1.0 and PEtab0.10.0. This is not a claim of Linux execution,
offline reproducibility, complete upstream CI coverage, or every PEtab-v1
configuration. `PACKAGING_STATUS.json` records the separate packaging replay.

## Independent check and limits

For x=x0+rate*t, y=y0+t and observed ratio h=x/(x+y), profile variance
v=mean((h-d)^2) over the declared observed set. The Gaussian negative
log-likelihood is n/2*(1+log(2*pi*v)), plus the separate ordinary quantitative
measurement contribution. `gaussian_oracle` evaluates this formula directly.
The envelope gradient uses dh/drate=t*(y0+t)/(x0+y0+(rate+1)*t)^2.
The configured second-order return is a **conditional Fisher information**, not
the exact Hessian of the profiled objective.

The 39-case native suite covers three non-finite utility values at three
positions, missing/empty masks, and four real compiled models across all seven
nonempty subsets of sensitivity orders{0,1,2}. The archived supplemental probe
also records empty requests, supported adjoint requests, explicitly unsupported
adjoint Hessians and real failed integration status−1. Supplemental archived
observations are evidence from the original run; the native replay command
executes the neutral 39-case regression, not that supplemental script.

Perfect-fit zero-noise MLE degeneracy, noise-bound policy, log-normal support,
all optimizers and novel biological interpretations are outside scope.
There is an explicit low-level boundary: calling `compute_nllh` manually with
zero sigma and a perfect fit changes its legacy return from −inf to NaN because
the new reductions preserve undefined terms. In the actual profiled perfect-fit
path, both versions still raise `AssertionError` while estimating sigma. This
repair does not resolve that degeneracy or claim every low-level edge return
remains identical.
The original method is [Efficient parameterization of large-scale dynamic models based on relative measurements; DOI10.1093/bioinformatics/btz581](https://pmc.ncbi.nlm.nih.gov/articles/PMC9883733/).
No claim is made that the original Boehm biological paper is wrong.

Early validation errors and the import-contamination correction are summarized
in `VALIDATION_HISTORY.md`. AI checks and passing tests are computational
evidence; no external human endorsement is implied.

Original audit code is MIT licensed. Archived upstream source and patch context
retain pyPESTO's BSD3-Clause notice under `licenses/`.
