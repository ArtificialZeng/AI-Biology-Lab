# Exact scope

The native defect is the failure of the continuous implication `ll < uu => ceil(ll) <= uu` used after integer lower-bound rounding. A feasible fixed integer value always obeys its implied upper bound, so the guard safely rejects the admitted strict contradictions. It changes no successful-path transformation data and does not affect continuous-LP preprocessing.

The first grid enumerates nonzero coefficients -12..12, integer column domains -2<=lb<ub<=2, and integer ranged row bounds -12<=lo<hi<=12. The independent 24,000-case grid adds a second integer variable, fixed continuous offsets, coupled rows, nonzero objectives, negative coefficients and denominator-32 exact data.

The inherited 1e-5 column-value tolerance remains outside the candidate's repair. Scaled near-boundary residuals are documented and tested, not omitted. There is no claim of all MIPs, all real coefficients, all tolerance policies or all platforms being correct. Presolve OFF after simplex, or splitting the minimal ranged row, provides correct infeasibility for the two explicit examples; this is an observed control, not a universal workaround theorem.

This does not settle the biological modeling contract for zero stoichiometric columns. A false feasible solver status and an invalid model assumption are separate questions. No external human acceptance or first-discovery claim is made.

The separate Optlang wrapper candidate checks only floating primal equality/bound diagnostics before exposing an optimal MIP status. Its absolute-or-relative tolerance policy can accept scale-small residuals, so it is not an exact feasibility certificate. It adds one disabled-presolve recovery per solver invocation and reports undefined after a second invalid optimal result. It preserves the preexisting enabled-presolve recovery path, shared auto retry and nonoptimal value-access behavior. It does not patch libGLPK, certify all integer models, or replace the native guard/residual scope above. See [wrapper/README.md](wrapper/README.md).
