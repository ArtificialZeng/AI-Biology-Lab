# Replay a GillesPy2 NumPy ODE failure-propagation candidate

**Patch status: accepted by a fresh independent AI correctness review for upstream consideration. Package review: accepted by a separate fresh AI packaging review.** This package demonstrates a scientific-software correctness defect and the reviewed correction. It does not claim a biological discovery, universal solver accuracy, trusted-runner certification, or external human peer review.

For the first-order reaction A → B with A(0)=1, B(0)=0 and rate one, the exact solution is A(t)=exp(−t), B(t)=1−exp(−t), and A+B=1. When SciPy cannot reach a requested output time, current GillesPy2's NumPy ODE wrapper can return an early state in a future-labelled row with status Success/0. For example, VODE step exhaustion gives A(2) near one instead of exp(−2)≈0.135335. This illustrative numerical value is runtime-dependent; it is not an archived numerical certificate.

The proposed two-file patch checks the backend success flag before recording that row, raises through the existing SimulationError path with requested/actual time and backend code in the cause, and resets per-run state so retries work. It does not change backend configuration or numerical defaults, or relabel partial failed data.

## Reproduce

Python 3.12 is the intended Linux/macOS CI runtime; local checks also use Python 3.13.5. Install the pinned runtime dependencies into an environment outside this project. SCons, a C++ compiler, libSBML and C++ solver builds are not needed for this NumPy replay. The complete vendored library retains its optional-backend import warning when dependencies are absent.

```sh
python -m pip install -r requirements.txt
python -B verify.py verify
python -O -B verify.py verify
```

Each command runs actual byte-preserved baseline and candidate libraries in separate processes, the independently written real-SciPy analytic spy, the seven new regression methods against both trees, and the complete 13-method candidate ODE module. The baseline must reproduce the false success and fail the new tests; the candidate must raise on failed steps and pass the controls. This is not a mock-only replay. Temporary logs and plotting caches are outside the project, and verification changes no archive input.

The deterministic certificate records checked properties and counts. It contains no raw backend times, warning text, temporary paths, OS/Python platform, or machine-specific numerical values. Floating observations are evaluated with explicit tolerances, and genuine failed-step codes/call counts are checked. See [verification scope](VERIFICATION.md), [source provenance](source-provenance.json), [rights and attribution](THIRD_PARTY.md), and [candidate patch](minimal.patch).

## Archive operations

`compute` emits the recomputed property certificate. Explicit regeneration recomputes it and proposes a certificate plus complete manifest in a new absolute directory outside this project:

```sh
python -B verify.py compute
python -B verify.py regenerate --output-dir /absolute/outside/path/new-ode-proposal
```

The parent must exist and the output must not exist. The proposal references the unchanged project inputs; it is not a copied source tree. It never overwrites this archive. Changing the patch requires new review; new hashes alone are not a correctness argument.

Main baseline: f1957c65b6cb318bdd3701af14297892583b4297. The inspected develop head is 382d160667b3138ac8ca4b4444f5dc82869783a8. This package vendors the complete 138-file gillespy2 library for each tree plus the two required upstream test modules and GPL notices, all checked against the frozen research source hashes. No generated builds or caches are included.

The original research run also passed 192 core unit tests; this smaller package does not replay that full core suite. The full cross-solver pause/resume module was not run because it imports C++ solvers requiring unavailable SCons. ODE-specific resume, duplicate trajectories and retry paths are covered. Timeout handling and other solvers are outside the proposed correction. No full-suite-green claim is made.

The fresh patch review independently reconstructed the official source and tested all five documented integrators in two analytic models, with 32 observations per source tree. That separate review found no mandatory correctness revision. This package replays the narrower finite cases listed in VERIFICATION.md; its certificate itself does not imply review acceptance. No maintainer or external human approval is claimed.
