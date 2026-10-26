# Native GLPK integer presolver audit

Official GNU GLPK 5.0 can report an exactly infeasible integer model as optimal. The minimal example is a binary variable `x` with `1 <= 11*x <= 10`: neither permitted value is feasible, but MIP presolve returns `OPT, x=1`, violating the original upper row bound by one.

Pure C reproduction from the official source distribution isolates the problem from Optlang, SWIG and COBRApy. In `npp_ineq_singlet`, rounding the implied lower bound fixes the integer variable at its old upper bound, after which the code drops the contradictory implied upper bound. An exact construction and the local [branch guard](native-guard.patch) are documented in [SCOPE.md](SCOPE.md).

The patch is a **limited candidate**. Its inherited `1e-5` column tolerance still admits the exactly infeasible scaled case `250000 <= 1000000*x <= 999995`. It returns `x=1` and violates the original row by **5**. This residual is deliberately retained in the verifier; the patch is not presented as a complete repair of all ranged integer models.

Two independently implemented integer-enumeration experiments give:

| Admitted grid | Original wrong optimal statuses | Candidate |
|---|---:|---:|
| 72,000 one-variable models | 2,600 | 0 |
| 24,000 models with fixed offsets, coupling and objectives | 12,030 | 0 |

All 6,134 feasible independent cases recover the enumerated optimum and original rows correctly. These are finite checks on exact rational inputs, with floating solver operations. A separate AI source review confirmed the mechanism and identified the residual; external maintainer review and historical priority remain pending. No biological mechanism is claimed.

## Reproduce

Python 3.12+, Git, a C compiler and Make are required. Supply the [official GNU archive](https://ftp.gnu.org/gnu/glpk/glpk-5.0.tar.gz); the script verifies its SHA-256 before extracting it.

```sh
python projects/glpk-presolver-auditor/reproduce_native.py \
  --archive /path/to/glpk-5.0.tar.gz \
  --output artifacts/glpk-native/receipt.json --logs artifacts/glpk-native/logs
```

The input archive and archived project stay unchanged. Builds run in temporary directories. Default output is stdout; a receipt path must be outside this project and must not already exist. The verifier checks both admitted grids, feasible recovery, and the known residual/control cases. It does not download, install or replace a system solver.

The GNU project identifies a bug-report mailing list; no canonical GNU GitHub PR target was verified. Mirrors are not described as GNU upstream. [Sources](SOURCE.md) · [Licensing](THIRD_PARTY.md).

## Separate Optlang wrapper candidate

The [wrapper audit](wrapper/README.md) packages a GLPK-only primal residual check and bounded disabled-presolve recovery in Optlang. It leaves the native library unchanged, preserves Optlang's existing enabled-presolve fallback, and retains the native guard's known residual above. Its portable verifier pins the Optlang commit and dependencies, imports temporary baseline/candidate trees, and checks the expected 27 failing baseline behaviors against 62 passing candidate tests. The Linux/macOS wrapper CI is separate from the native C workflow.
