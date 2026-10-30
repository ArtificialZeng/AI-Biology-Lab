# Frozen primary sources — 26 October 2026

- [GNU GLPK project](https://www.gnu.org/software/glpk/) and [official 5.0 release](https://ftp.gnu.org/gnu/glpk/glpk-5.0.tar.gz), SHA-256 `4a1013eebb50f728fc601bdd833b0b2870333c3b3e5a816eeba921d95bec6f15`.
- The implicated function is `src/npp/npp3.c:npp_ineq_singlet`, with implied integer-bound handling in the same file. Source return 4 is propagated as preprocessing infeasibility. The candidate keeps the surrounding integer tolerance policy unchanged.
- [Optlang pinned wrapper](https://github.com/opencobra/optlang/blob/6427c6c10220ebff56b1e72591bc4ee43549d81b/src/optlang/glpk_interface.py), identical to the installed 1.9.1 wrapper during isolation, is a separate downstream exposure. Native C reproduction does not import it.

The bounded bug-list/source search found no inspected matching repair; that is not a priority certificate. GNU's project page lists its bug mailing list. A canonical GitHub PR repository was not verified; GitHub mirrors and distribution packaging repositories are not represented as GNU upstream.

For the wrapper candidate, Optlang source is pinned to commit `6427c6c10220ebff56b1e72591bc4ee43549d81b`, original `glpk_interface.py` SHA-256 `7932237d98ee07c344cc6cf12d067aa8b0496132a46d83fce5c45893e500b8e8`. The official GNU 5.0 manual's `glp_check_kkt` subsection explicitly permits MIP primal equality/bound checks; dual checks are not called. Actual relative-error normalization is grounded in `src/draft/glpapi10.c`. The wrapper's configured-feasibility threshold is a stated downstream policy, not a claim that native MIP preprocessing internally applies simplex `tol_bnd` identically. Frozen swiglpk/pytest/SymPy dependencies are in `wrapper/requirements.txt`.
