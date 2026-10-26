# COBRApy loopless analysis: empty-selection repair

At frozen development commit
[`a9c5ced`](https://github.com/opencobra/cobrapy/tree/a9c5ced29670bed5035fb19e87e7b0b588adb62d),
`add_loopless` can raise `IndexError` on an ordinary feasible acyclic metabolic
pathway. The cyclic finder correctly returns no internal cycles, then an empty
NumPy array is used as a floating-point index. A boundary-only model triggers
the same problem in `find_cyclic_reactions`, also affecting loopless FVA.

Two guards, six production-source lines, handle the vacuous cases. The proposed
upstream patch also adds regression tests and a release note. A basic uptake →
A → B → drain model needs no empty biochemical reaction to reproduce the defect.

## Reproduce without changing the source checkout

Use a clean Python 3.11 or 3.12 environment and clone the frozen upstream:

```sh
python -m pip install -r projects/cobra-loopless-auditor/requirements.txt
git clone https://github.com/opencobra/cobrapy.git /tmp/cobrapy-source
git -C /tmp/cobrapy-source checkout a9c5ced29670bed5035fb19e87e7b0b588adb62d
python projects/cobra-loopless-auditor/reproduce.py --checkout /tmp/cobrapy-source
```

The runner checks source bytes, applies the patch in a temporary copy, and runs
an independent diamond-network contract suite against both sources. Expected:
baseline **24 failed / 23 passed**, patched **47 passed**. It verifies the exact
exception type, not just a failing exit code. `--logs DIR` preserves separate
logs; `--output PATH` writes a separate JSON receipt. The supplied checkout and
archived evidence are left unchanged.

The independently constructed network has analytic optimum 9 and nontrivial
FVA bounds. Controls cover a forced loop, partial/empty selections, generators,
both loopless formulations, both FVA algorithms, user constraints and context
restoration. The affected upstream suites also pass **37 tests**, with 9
benchmark/parallel cases excluded, under GLPK.

[Patch](upstream.patch) · [archived receipt](evidence/certificate.json) · [credit and licensing](THIRD_PARTY.md)

This is a software correctness repair. It does not establish a new biological
mechanism, the prevalence of affected real models, or correctness under every
solver. The separate zero-stoichiometric-column/GLPK-status investigation is
parked and is not part of this fix. The inspected PR #1472 branch does not
contain these guards; broader priority and maintainer acceptance remain pending.
