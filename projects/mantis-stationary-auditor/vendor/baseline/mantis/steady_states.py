"""Exhaustive steady-state enumeration for mass-action networks.

The hybrid solver in :mod:`mantis.analysis` (ODE integration + multi-start
``least_squares``) *can silently miss* steady states — in particular unstable
fixed points that forward integration flows away from.  At mass-action
equilibrium the positive steady states are the positive real roots of a
polynomial system, so they can instead be found *exhaustively* with completeness
guarantees:

* **Default pure-Python backend** (``backend="sympy"``).  The reaction-rate
  polynomials are reduced against the conservation laws, a lexicographic Gröbner
  basis triangularises the system, every complex root is recovered by
  back-substitution (``numpy.roots``) and refined by a Newton polish, and the
  real non-negative roots are returned.  For a zero-dimensional ideal this finds
  *all* solutions (their count is bounded by the Bézout number), so it certifies
  the full steady-state set rather than sampling it.
* **Optional homotopy backend** (``backend="phcpy"``).  Delegates to PHCpack via
  :mod:`phcpy` (Verschelde) — numerical polynomial homotopy continuation that
  scales to larger / stiffer systems.  Requires ``pip install mantis-delta[homotopy]``.

The classic demonstration is the Schlögl model: this enumerator returns all
three positive steady states (including the unstable middle branch that
integration never reaches), whereas the multi-start solver typically returns
only the two stable ones.
"""
from __future__ import annotations

import numpy as np

from .parsing import Reaction
from .stoichiometry import build_stoichiometry_matrix, conservation_laws_sympy, matrix_rank


def _scaled_polynomial_system(reactions, species, rate_values, initial_conditions,
                              chemostatted_values, scale):
    """Build the rationalised, non-dimensionalised steady-state polynomial system.

    Returns ``(reduced_eqs, free_syms, x_syms, param, scale)`` where ``param`` maps
    every scaled species variable to its expression in the free variables (after
    eliminating the conservation laws), or the identity when there are none.
    """
    import sympy

    from .symbolic import (
        make_species_symbols, make_rate_symbols, build_odes, substitute_rates,
    )

    chem = chemostatted_values or {}
    sp_syms = make_species_symbols(species)
    rate_syms, key_to_sym = make_rate_symbols(reactions)
    odes_sym = build_odes(reactions, species, sp_syms, rate_syms, key_to_sym,
                          chemostatted_values=chem or None)
    odes = {sp: substitute_rates(expr, key_to_sym, rate_values) for sp, expr in odes_sym.items()}

    N = build_stoichiometry_matrix(reactions, species)
    rank = matrix_rank(N)
    laws = conservation_laws_sympy(N, species)

    # Characteristic concentration scale → variables of order 1 (well conditioned).
    ic_subs = {sympy.Symbol(s): float(initial_conditions.get(s, 0.0)) for s in species}
    if scale is None:
        totals = [abs(float(law.subs(ic_subs))) for law in laws]
        cands = [t for t in totals if t > 1e-300]
        cands += [float(v) for v in initial_conditions.values() if v and v > 0]
        cands.append(1.0)
        scale = max(cands)
    C = sympy.Rational(scale).limit_denominator(10 ** 9)

    syms = [sympy.Symbol(s) for s in species]
    xs = [sympy.Symbol(f"_x{i}") for i in range(len(species))]
    sub = {syms[i]: C * xs[i] for i in range(len(species))}

    eqs = []
    for s in species:
        e = sympy.expand(odes[s].subs(sub))
        e = sympy.nsimplify(e, rational=True)
        poly = sympy.Poly(e, *xs) if e != 0 else None
        if poly is not None and poly.coeffs():
            mc = max(abs(c) for c in poly.coeffs())
            e = sympy.expand(e / mc)
        eqs.append(e)

    # Pick rank(N) independent stoichiometric rows (QR column pivoting on Nᵀ).
    if rank > 0:
        from scipy.linalg import qr
        _, _, perm = qr(N.T, pivoting=True)
        ind = list(perm[:rank])
    else:
        ind = []

    if laws:
        cons = []
        for law in laws:
            tot = float(law.subs(ic_subs))
            cons.append(sympy.nsimplify(sympy.expand(law.subs(sub) - tot), rational=True))
        sol_set = sympy.linsolve(cons, xs)
        if not sol_set:
            raise ValueError("Conservation totals are infeasible for this network.")
        param_tuple = list(sol_set)[0]
        free_syms = sorted(param_tuple.free_symbols, key=lambda s: s.name)
        consub = {xs[i]: param_tuple[i] for i in range(len(xs))}
        param = {xs[i]: param_tuple[i] for i in range(len(xs))}
        reduced = []
        for i in ind:
            r = sympy.expand(eqs[i].subs(consub))
            if r != 0:
                reduced.append(sympy.nsimplify(r, rational=True))
    else:
        free_syms = list(xs)
        param = {x: x for x in xs}
        reduced = [eqs[i] for i in ind if eqs[i] != 0]

    return reduced, free_syms, xs, param, float(C)


def _groebner_all_roots(reduced, free_syms):
    """All complex solutions of a 0-dimensional system via lex Gröbner triangulation."""
    import sympy

    if not free_syms:
        return [{}]
    if not reduced:
        raise ValueError(
            "Steady-state set is positive-dimensional (a continuum of equilibria); "
            "supply conservation totals that pin it down, or use backend='phcpy'."
        )

    basis = sympy.groebner(reduced, *free_syms, order="lex")
    polys = list(basis.polys)

    solutions = [{}]
    for var in reversed(free_syms):
        next_solutions = []
        for partial in solutions:
            univariate = None
            for p in polys:
                expr = p.as_expr().subs(partial)
                if (expr.free_symbols & set(free_syms)) == {var}:
                    univariate = expr
                    break
            if univariate is None:
                raise ValueError(
                    "Could not triangularise the system (positive-dimensional or "
                    "degenerate); try backend='phcpy'."
                )
            coeffs = [complex(c) for c in sympy.Poly(univariate, var).all_coeffs()]
            if len(coeffs) <= 1:
                continue
            for root in np.roots(coeffs):
                child = dict(partial)
                child[var] = complex(root)
                next_solutions.append(child)
        solutions = next_solutions
    return solutions


def _newton_polish(reduced, free_syms, z0, max_iter=60):
    """Refine a complex root of ``reduced`` with Newton's method; return (z, residual)."""
    import sympy

    F = sympy.lambdify(free_syms, reduced, "numpy")
    J = sympy.Matrix(reduced).jacobian(free_syms)
    Jf = sympy.lambdify(free_syms, J, "numpy")
    n = len(free_syms)
    z = np.array(z0, dtype=complex)
    for _ in range(max_iter):
        fv = np.array(F(*z), dtype=complex).flatten()
        res = float(np.linalg.norm(fv))
        if res < 1e-14:
            break
        Jv = np.array(Jf(*z), dtype=complex).reshape(len(reduced), n)
        try:
            dz = np.linalg.lstsq(Jv, fv, rcond=None)[0]
        except np.linalg.LinAlgError:
            break
        z = z - dz
    fv = np.array(F(*z), dtype=complex).flatten()
    return z, float(np.linalg.norm(fv))


def _sympy_all_roots(reactions, species, rate_values, initial_conditions,
                     chemostatted_values, scale, residual_tol, imag_tol):
    """Return all positive real steady-state concentration vectors (pure-Python)."""
    reduced, free_syms, xs, param, C = _scaled_polynomial_system(
        reactions, species, rate_values, initial_conditions, chemostatted_values, scale
    )
    raw = _groebner_all_roots(reduced, free_syms)

    results: list[np.ndarray] = []
    for sol in raw:
        if free_syms:
            z0 = [sol[s] for s in free_syms]
            z, res = _newton_polish(reduced, free_syms, z0)
            if res > residual_tol:
                continue
            fsub = {free_syms[i]: z[i] for i in range(len(free_syms))}
        else:
            fsub = {}
        full = np.array([complex(param[xs[i]].subs(fsub)) for i in range(len(species))])
        # Keep real, non-negative roots only.
        if np.any(np.abs(full.imag) > imag_tol * (np.abs(full.real) + 1.0)):
            continue
        real = full.real
        if np.any(real < -1e-6 * (np.max(np.abs(real)) + 1e-30)):
            continue
        results.append(np.maximum(real, 0.0) * C)

    # Deduplicate (relative L2 < 0.1%).
    unique: list[np.ndarray] = []
    for y in results:
        if not any(np.linalg.norm(y - u) <= 1e-3 * (np.linalg.norm(u) + 1e-30) for u in unique):
            unique.append(y)
    return unique


def _phcpy_all_roots(reactions, species, rate_values, initial_conditions,
                     chemostatted_values, scale, residual_tol, imag_tol):  # pragma: no cover
    """Optional homotopy-continuation backend via PHCpack/phcpy (experimental)."""
    try:
        from phcpy.solver import solve as phc_solve
        from phcpy.solutions import strsol2dict
    except ImportError as exc:
        raise ImportError(
            "backend='phcpy' requires the optional homotopy backend. "
            "Install with:  pip install mantis-delta[homotopy]"
        ) from exc

    reduced, free_syms, xs, param, C = _scaled_polynomial_system(
        reactions, species, rate_values, initial_conditions, chemostatted_values, scale
    )
    import sympy

    if not free_syms:
        fsub = {}
        full = np.array([complex(param[xs[i]].subs(fsub)) for i in range(len(species))])
        return [np.maximum(full.real, 0.0) * C]

    var_names = [s.name for s in free_syms]
    eq_strings = [str(sympy.expand(eq)).replace("**", "^") + ";" for eq in reduced]
    sols = phc_solve(eq_strings, verbose=False)

    results = []
    for s in sols:
        d = strsol2dict(s)
        z = np.array([complex(d[v]) for v in var_names])
        fsub = {free_syms[i]: z[i] for i in range(len(free_syms))}
        full = np.array([complex(param[xs[i]].subs(fsub)) for i in range(len(species))])
        if np.any(np.abs(full.imag) > imag_tol * (np.abs(full.real) + 1.0)):
            continue
        if np.any(full.real < -1e-6 * (np.max(np.abs(full.real)) + 1e-30)):
            continue
        results.append(np.maximum(full.real, 0.0) * C)

    unique = []
    for y in results:
        if not any(np.linalg.norm(y - u) <= 1e-3 * (np.linalg.norm(u) + 1e-30) for u in unique):
            unique.append(y)
    return unique


def all_steady_states(
    reactions: list[Reaction],
    species: list[str],
    rate_values: dict[str, float],
    initial_conditions: dict[str, float],
    chemostatted_values: dict[str, float] | None = None,
    backend: str = "auto",
    scale: float | None = None,
    residual_tol: float = 1e-8,
    imag_tol: float = 1e-6,
):
    """Exhaustively enumerate the positive steady states of a mass-action network.

    Returns a list of :class:`~mantis.analysis.SteadyState` objects (with
    eigenvalues and stability classification), sorted by residual.  See the module
    docstring for the algorithm and the meaning of ``backend``.
    """
    from scipy.linalg import eigvals
    from .analysis import SteadyState, _full_jacobian_fn, classify_steady_state, build_ode_function
    from .stoichiometry import fold_chemostatted_into_rates

    chem_vals = chemostatted_values or {}
    chem_keys = set(chem_vals.keys())
    eff_rates = (
        fold_chemostatted_into_rates(reactions, rate_values, chem_vals)
        if chem_vals else rate_values
    )

    if backend == "phcpy":
        roots = _phcpy_all_roots(reactions, species, rate_values, initial_conditions,
                                 chem_vals, scale, residual_tol, imag_tol)
    elif backend in ("auto", "sympy"):
        try:
            roots = _sympy_all_roots(reactions, species, rate_values, initial_conditions,
                                     chem_vals, scale, residual_tol, imag_tol)
        except Exception:
            if backend == "auto":
                try:
                    roots = _phcpy_all_roots(reactions, species, rate_values,
                                             initial_conditions, chem_vals, scale,
                                             residual_tol, imag_tol)
                except Exception:
                    raise
            else:
                raise
    else:
        raise ValueError(f"Unknown backend {backend!r}; use 'auto', 'sympy', or 'phcpy'.")

    f_ode = build_ode_function(reactions, species, eff_rates, chem_keys)
    full_jac = _full_jacobian_fn(reactions, species, eff_rates, chem_keys)

    states: list[SteadyState] = []
    for y in roots:
        eigs = eigvals(full_jac(y))
        is_stable, is_osc = classify_steady_state(eigs)
        states.append(SteadyState(
            concentrations=dict(zip(species, y.tolist())),
            eigenvalues=eigs,
            is_stable=is_stable,
            is_oscillatory=is_osc,
            residual=float(np.linalg.norm(f_ode(0, y))),
        ))
    states.sort(key=lambda s: s.residual)
    return states
