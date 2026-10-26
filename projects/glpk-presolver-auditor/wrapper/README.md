# Optlang GLPK MIP residual validation

This is a separate downstream wrapper candidate. It **does not modify or upgrade libGLPK** and does not remove the native branch guard's documented near-tolerance residual. The native C audit, its limited guard and the scaled residual remain independently reproducible from the parent project.

The pinned Optlang wrapper can expose native `GLP_OPT` as `optimal` even when an original ranged row is violated. The two binary examples are `1 <= 11*x <= 10` and `250000 <= 1000000*x <= 999995`; exact enumeration makes both infeasible. The latter remains a residual of the native guard, with row violation 5. The wrapper candidate catches both at the default feasibility tolerance.

The [production patch](optlang-mip-validation.patch) checks original MIP primal equalities and bounds through GLPK's documented `glp_check_kkt(GLP_MIP, GLP_KKT_PE/PB)`. A condition fails only when **both** its maximum absolute and maximum relative errors exceed `configuration.tolerances.feasibility`. This preserves legitimate large-scale roundoff. It is a floating diagnostic, not an exact certificate or an independent integrality/global-optimality check.

Reported residual maxima must also be finite and nonnegative; nonfinite diagnostics are rejected even if the other error is zero. This check still relies on the native diagnostic: extreme overflow inside GLPK can conceal an error by reporting `(0,0)`. An independent review demonstrated that limitation on an explicitly loaded bad assignment with extreme finite coefficients, not on a solver-generated result. No full exact recomputation is claimed here.

The existing enabled-presolve fallback after undefined/infeasible results is preserved, including when configuration is False. If the final native optimal assignment fails validation, the wrapper makes one direct MIP retry with presolve disabled and restores configuration in `finally`. A second invalid optimal assignment becomes `undefined`, never an inferred infeasibility. A legitimate status returned by the retry retains its meaning. Value access after nonoptimal statuses is unchanged and can expose stale native values.

The new recovery is bounded to one per `_optimize` invocation. Counting the old fallback gives at most three MIP calls per invocation. The unchanged shared auto retry can invoke it twice, giving a worst-case bound of six calls; existing per-call timeout settings are retained. The minimal contradiction actually takes four calls with auto, two with True and three with False.

## Portable reproduction

Use Python 3.11 or 3.12 and Git. Install the frozen minimal dependencies into a separate environment, then supply a local Optlang clone containing the pinned commit. No network, dependency installation or system-solver replacement occurs inside the verifier.

```sh
python -m venv .venv-wrapper
.venv-wrapper/bin/python -m pip install -r projects/glpk-presolver-auditor/wrapper/requirements.txt
git clone https://github.com/opencobra/optlang.git artifacts/optlang-source
.venv-wrapper/bin/python projects/glpk-presolver-auditor/reproduce_wrapper.py \
  --checkout artifacts/optlang-source \
  --output artifacts/glpk-wrapper/receipt.json --logs artifacts/glpk-wrapper/logs
```

The verifier uses `git archive` of commit `6427c6c10220ebff56b1e72591bc4ee43549d81b`, ignoring checkout modifications. It verifies the original GLPK file hash, builds temporary baseline/candidate copies, applies only the wrapper patch, and checks that each imported interface comes from its intended temporary source. Every child receives exactly seven selected environment keys. Test traceback mode is short; saved logs replace temporary, checkout and project paths. Existing receipts/log directories cannot be overwritten, and outputs must stay outside the archived project.

The fixed targeted outcomes are:

| Suite | Result |
|---|---|
| Baseline public behavior | 27 failed, 27 passed, 8 candidate-helper checks deselected |
| Candidate, including helper checks | 62 passed |

The tests cover signed and scaled contradictions, exact feasible coupled integer objectives and offsets, configuration restoration after recovery/exception, second-invalid-result handling, and a large-row roundoff control for auto/False/True. The latter has equality residual about `1.22e-4` absolute but `6.1e-17` relative and is deliberately accepted.

Existing GLPK and exact-interface compatibility suites were separately observed at 187 passed, 10 skipped and two existing clone-with-LP warnings; their receipt is separate from the targeted verifier. They are not a required repeat in the lightweight wrapper CI. Two independently written native C grids remain a separate validation layer. Neither passing tests nor AI source review constitutes external human peer review or a biological discovery.

The wrapper patch and tests use Apache-2.0, matching Optlang; see [LICENSE.optlang](../LICENSE.optlang). The original verifier and explanations use the lab's MIT license. The official GNU reporting route remains the bug mailing list; no GitHub mirror is described as GNU upstream.
