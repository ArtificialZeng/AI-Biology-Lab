"""Parameter regions for multistationarity — Conradi, Feliu, Mincheva & Wiuf
(*PLoS Comput. Biol.* 2017).

Injectivity (:mod:`mantis.injectivity`) answers a yes/no question — *can* this
network ever have more than one positive steady state in a compatibility class?
Conradi *et al.* (2017) go further and ask *where in parameter space* it does:
they give the determinant criterion that carves the rate-/total-concentration
space into a monostationary region and a multistationary region.

The object of study is the **critical function**

    φ(c, k) = (−1)^s · det J(c, k),

where ``J`` is the Jacobian of the canonical steady-state map (``s = rank(N)``
independent mass-action ODE rows stacked on the ``n − s`` conservation laws) and
``s`` is the dimension of the stoichiometric subspace.  Reading ``φ`` as a
polynomial in the steady-state concentrations ``c`` with coefficients that are
polynomials in the rate constants ``k``:

* If **every coefficient is positive** for all positive ``k`` then ``φ > 0`` on
  the whole positive orthant, the steady-state map is injective, and the network
  is **monostationary** for *all* rates and totals (Craciun–Feinberg).
* If **some coefficient can be negative**, ``φ`` attains negative values and the
  network is **multistationary** for parameters in that region.  Each such
  coefficient inequality ``a_α(k) < 0`` is a face of the multistationarity
  region — the explicit "where" that the 2017 paper contributes.

For a network whose positive steady states admit a monomial parametrisation
(weakly reversible / toric steady states, the paper's setting) a negative
coefficient is both necessary *and* sufficient for multistationarity; in general
it is the necessary condition (a sign-change witness).  This module returns the
critical function, the per-coefficient sign certificate, and — when numeric
rates are supplied — whether those specific rates land in the multistationary
region.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .parsing import Reaction
from .stoichiometry import build_stoichiometry_matrix, conservation_laws_sympy, matrix_rank


@dataclass
class MultistationarityResult:
    """Outcome of the Conradi–Feliu–Mincheva–Wiuf critical-function analysis.

    Attributes
    ----------
    monostationary : bool
        True ⇒ certified at most one positive steady state per compatibility
        class for *all* positive rate constants (the critical function is
        sign-definite positive).
    multistationary_possible : bool
        True ⇒ some coefficient of the critical function can be negative, so
        multistationarity occurs for parameters in ``region_conditions``.
    critical_function : sympy.Expr
        ``φ = (−1)^s · det J`` in concentrations and rate constants.
    coefficient_terms : list[tuple]
        ``(monomial, coefficient_expr, sign)`` for each concentration-monomial of
        ``φ``; ``sign`` ∈ {``"+"``, ``"−"``, ``"±"``} is the coefficient's sign
        over positive rate constants.
    region_conditions : list[sympy.Expr]
        Coefficient expressions whose negativity delimits the multistationarity
        region (``expr < 0`` ⇒ multistationarity).  Empty when monostationary.
    multistationary_at_rates : bool | None
        If numeric rates were supplied: whether they actually yield ≥ 2 positive
        steady states.  When the steady-state set reduces to a univariate
        polynomial this is *exact* (a positive-root count); otherwise it falls
        back to the necessary sign condition (a critical coefficient is negative).
        ``None`` if no rates were supplied.
    region_boundary : sympy.Expr | None
        The **exact** multistationarity-region boundary when the steady states
        reduce to a univariate polynomial ``p(x)``: the discriminant
        ``disc_x(p)``.  Crossing ``region_boundary = 0`` changes the number of
        (positive) steady states, so it is the genuine fold locus in parameter
        space — unlike ``region_conditions``, which is only the necessary sign
        condition.  ``None`` when the reduction is unavailable.
    steady_state_polynomial : sympy.Expr | None
        The reduced univariate steady-state polynomial ``p(x)`` (in one species,
        with rate / total parameters), when available.
    summary_lines : list[str]
    """
    monostationary: bool
    multistationary_possible: bool
    critical_function: object = None
    coefficient_terms: list = field(default_factory=list)
    region_conditions: list = field(default_factory=list)
    multistationary_at_rates: "bool | None" = None
    summary_lines: list[str] = field(default_factory=list)
    region_boundary: object = None
    steady_state_polynomial: object = None

    def __str__(self) -> str:
        return "\n".join(self.summary_lines)


def _coefficient_sign(expr, positive_syms) -> str:
    """Sign of a polynomial ``expr`` over the positive orthant of ``positive_syms``.

    Returns ``"+"`` / ``"−"`` if every monomial coefficient shares a sign (so the
    expression is sign-definite for positive arguments), else ``"±"``.
    """
    import sympy

    expr = sympy.expand(expr)
    if expr == 0:
        return "0"
    if not expr.free_symbols:
        v = float(expr)
        return "+" if v > 0 else ("−" if v < 0 else "0")
    syms = [s for s in positive_syms if s in expr.free_symbols]
    if not syms:
        v = float(expr)
        return "+" if v > 0 else ("−" if v < 0 else "0")
    poly = sympy.Poly(expr, *syms)
    signs = {int(np.sign(float(c))) for c in poly.coeffs() if c != 0}
    if signs == {1}:
        return "+"
    if signs == {-1}:
        return "−"
    return "±"


def _steady_state_polynomial(reactions, species, sp_syms, rate_syms, key_to_sym,
                             dyn_syms, N, s):
    """Reduce the steady-state system to a single univariate polynomial, if possible.

    Eliminates the ``n − s`` conservation-law species (with symbolic totals ``T_j``)
    from the ``s`` independent mass-action steady-state equations.  Returns
    ``(poly, var, totals)`` when this collapses to one univariate polynomial — the
    common case being ``s = 1`` (one-dimensional stoichiometric subspace, e.g.
    Schlögl) — else ``None``.  Restricting the analysis to this polynomial puts the
    critical-function test *on the steady-state variety*, where a sign change is
    genuinely equivalent to multistationarity (Conradi *et al.* 2017).
    """
    import sympy

    from .symbolic import build_odes

    if s <= 0:
        return None

    # Mass-action ODEs with symbolic rates and symbolic chemostatted concentrations.
    odes = build_odes(reactions, species, sp_syms, rate_syms, key_to_sym,
                      chemostatted_values=None)

    from scipy.linalg import qr
    _, _, perm = qr(N.T, pivoting=True)
    ind = list(perm[:s])
    ode_eqs = [sympy.expand(odes[species[i]]) for i in ind]

    laws = conservation_laws_sympy(N, species)
    totals = [sympy.Symbol(f"T_{j + 1}") for j in range(len(laws))]
    free = list(dyn_syms)
    if laws:
        cons = [law - T for law, T in zip(laws, totals)]
        sol = sympy.linsolve(cons, dyn_syms)
        if not sol:
            return None
        tup = list(sol)[0]
        subs = {dyn_syms[i]: tup[i] for i in range(len(dyn_syms))}
        free = sorted(
            set().union(*[t.free_symbols for t in tup]) & set(dyn_syms),
            key=lambda x: x.name,
        )
        ode_eqs = [sympy.expand(e.subs(subs)) for e in ode_eqs]

    if len(free) == 1:
        var = free[0]
        num = sympy.together(ode_eqs[0]).as_numer_denom()[0]
        p = sympy.Poly(sympy.expand(num), var)
        return p, var, totals

    # Higher-dimensional: try a lex Gröbner elimination to one univariate generator.
    if 2 <= len(free) <= 3:
        try:
            nums = [sympy.together(e).as_numer_denom()[0] for e in ode_eqs]
            gb = sympy.groebner(nums, *free, order="lex")
            var = free[-1]
            for g in gb.polys:
                expr = g.as_expr()
                if expr.free_symbols & set(free) == {var}:
                    return sympy.Poly(expr, var), var, totals
        except Exception:
            return None
    return None


def _count_positive_roots(poly, var, subs):
    """Number of positive real roots of ``poly`` after substituting ``subs``.

    Returns ``None`` if the substitution leaves free parameters (e.g. unspecified
    conservation totals), so the count cannot be evaluated.
    """
    import sympy

    expr = poly.as_expr().subs(subs)
    if expr.free_symbols - {var}:
        return None
    coeffs = [complex(c) for c in sympy.Poly(expr, var).all_coeffs()]
    if len(coeffs) <= 1:
        return 0
    roots = np.roots(coeffs)
    return int(sum(1 for r in roots
                   if abs(r.imag) <= 1e-7 * (abs(r) + 1.0) and r.real > 1e-9))


def multistationarity_region(
    reactions: list[Reaction],
    species: list[str],
    rate_values: dict[str, float] | None = None,
    chemostatted_values: dict[str, float] | None = None,
) -> MultistationarityResult:
    """Compute the multistationarity parameter region via the CFMW critical function."""
    import sympy

    from .symbolic import make_species_symbols, make_rate_symbols, mass_action_flux

    chem = chemostatted_values or {}
    chem_keys = set(chem)

    # Symbols for every species in the reactions; chemostatted ones enter φ as
    # positive parameters (never folded to a number, so they appear in the region).
    all_names: list[str] = []
    seen: set[str] = set()
    for rxn in reactions:
        for name, _ in list(rxn.reactants) + list(rxn.products):
            if name not in seen:
                seen.add(name)
                all_names.append(name)
    sp_syms = make_species_symbols(sorted(all_names))
    rate_syms, key_to_sym = make_rate_symbols(reactions)

    dyn_syms = [sp_syms[s] for s in species]
    n = len(species)
    N = build_stoichiometry_matrix(reactions, species)
    s = matrix_rank(N)

    if n == 0:
        return MultistationarityResult(
            True, False, sympy.Integer(0), [], [], None,
            ["Multistationarity (Conradi–Feliu–Mincheva–Wiuf): no dynamic species."],
        )

    fluxes = [mass_action_flux(rxn, sp_syms, key_to_sym[rxn.rate_key], None)
              for rxn in reactions]

    # rank(N) independent ODE rows (QR column pivot on Nᵀ) + the conservation laws.
    if s > 0:
        from scipy.linalg import qr
        _, _, perm = qr(N.T, pivoting=True)
        ind = list(perm[:s])
    else:
        ind = []

    rows = []
    for i in ind:
        expr = sympy.Integer(0)
        for j in range(len(reactions)):
            coeff = int(round(N[i, j]))
            if coeff:
                expr += coeff * fluxes[j]
        rows.append(sympy.expand(expr))
    rows += list(conservation_laws_sympy(N, species))

    if len(rows) != n:
        return MultistationarityResult(
            False, True, None, [], [], None,
            ["Multistationarity: inconclusive (degenerate stoichiometry)."],
        )

    J = sympy.Matrix([[sympy.diff(r, c) for c in dyn_syms] for r in rows])
    det = sympy.expand(J.det())
    phi = sympy.expand((-1) ** s * det)

    # Positive parameters: rate constants and chemostatted concentrations.
    chem_param_syms = [sp_syms[s_] for s_ in sorted(all_names) if s_ in chem_keys]
    positive_syms = list(rate_syms.values()) + chem_param_syms

    if phi == 0:
        return MultistationarityResult(
            False, True, phi, [], [], None,
            [
                "Multistationarity (Conradi–Feliu–Mincheva–Wiuf): INCONCLUSIVE",
                "  → The critical function vanishes identically (degenerate Jacobian).",
            ],
        )

    # φ as a polynomial in the concentrations; group coefficients (in k, chemostats).
    conc_poly = sympy.Poly(phi, *dyn_syms) if any(c in phi.free_symbols for c in dyn_syms) \
        else None
    coeff_terms = []
    region = []
    if conc_poly is not None:
        for monom, coeff in conc_poly.terms():
            mon_expr = sympy.prod([c ** e for c, e in zip(dyn_syms, monom)])
            sign = _coefficient_sign(coeff, positive_syms)
            coeff_terms.append((mon_expr, sympy.expand(coeff), sign))
            if sign in ("−", "±"):
                region.append(sympy.expand(coeff))
    else:
        sign = _coefficient_sign(phi, positive_syms)
        coeff_terms.append((sympy.Integer(1), phi, sign))
        if sign in ("−", "±"):
            region.append(phi)

    all_positive = all(t[2] in ("+", "0") for t in coeff_terms)
    monostationary = all_positive
    multistationary_possible = not all_positive

    # ── Exact region: restrict to the steady-state variety (univariate reduction) ─
    # Reading φ's coefficients over the whole positive orthant is only *necessary*
    # for multistationarity (φ < 0 may occur off the steady-state set).  When the
    # steady states reduce to a univariate polynomial p(x), multistationarity ⇔ p
    # has ≥ 2 positive roots, whose boundary is the discriminant disc_x(p) = 0 — the
    # exact fold locus.
    region_boundary = None
    ss_poly = None
    var = None
    try:
        poly_info = _steady_state_polynomial(
            reactions, species, sp_syms, rate_syms, key_to_sym, dyn_syms, N, s
        )
    except Exception:
        poly_info = None
    if poly_info is not None:
        p, var, _totals = poly_info
        if p.degree() >= 2:
            ss_poly = p.as_expr()
            try:
                region_boundary = sympy.factor(sympy.discriminant(p, var))
            except Exception:
                region_boundary = None

    # ── Locate the supplied numeric rates relative to the region ──────────────
    at_rates = None
    rate_subs = None
    if rate_values is not None:
        rate_subs = {key_to_sym[k]: sympy.Float(v)
                     for k, v in rate_values.items() if k in key_to_sym}
        rate_subs.update({sp_syms[s_]: sympy.Float(chem[s_])
                          for s_ in chem_keys if s_ in sp_syms})

    # Exact numeric verdict via positive-root count, when the polynomial is fully
    # determined by the supplied rates (no free conservation totals remain).
    exact_at_rates = False
    if rate_subs is not None and ss_poly is not None:
        npos = _count_positive_roots(p, var, rate_subs)
        if npos is not None:
            at_rates = npos >= 2
            exact_at_rates = True

    # Fallback: the necessary sign condition (a critical coefficient goes negative).
    if not exact_at_rates and rate_subs is not None and multistationary_possible:
        at_rates = False
        for _mon, coeff, _sign in coeff_terms:
            val = coeff.subs(rate_subs)
            if val.free_symbols:
                continue
            if float(val) < 0:
                at_rates = True
                break

    # ── Human-readable summary ───────────────────────────────────────────────
    if monostationary:
        lines = [
            "Multistationarity (Conradi–Feliu–Mincheva–Wiuf): MONOSTATIONARY",
            f"  → The critical function φ = (−1)^{s}·det J is sign-definite positive,",
            "    so it never vanishes for positive concentrations and rate constants.",
            "  → At most ONE positive steady state per compatibility class, for ALL",
            "    rate constants and total concentrations — multistationarity is",
            "    structurally impossible.",
        ]
    else:
        lines = [
            "Multistationarity (Conradi–Feliu–Mincheva–Wiuf): REGION EXISTS",
            f"  → The critical function φ = (−1)^{s}·det J has coefficient(s) that can",
            "    be negative, so multistationarity occurs in a region of parameter space.",
        ]
        if region_boundary is not None:
            # Exact boundary on the steady-state variety (univariate reduction).
            lines.append("  → Exact region boundary (steady states reduce to a univariate")
            lines.append(f"    polynomial in {var}): the discriminant")
            lines.append(f"        disc({var}) = {region_boundary}")
            lines.append("    vanishes at the folds; the sign of disc fixes the number of")
            lines.append("    (positive) steady states on either side.")
        else:
            lines.append("  → Necessary sign condition (φ attains a negative value) requires:")
            for coeff in region[:8]:
                lines.append(f"        {sympy.simplify(coeff)} < 0")
            if len(region) > 8:
                lines.append(f"        … and {len(region) - 8} more coefficient condition(s)")
        if at_rates is not None:
            if exact_at_rates:
                verdict = ("multistationary (≥ 2 positive steady states)"
                           if at_rates else
                           "monostationary (one positive steady state)")
                lines.append(f"  → At the supplied rates the network is {verdict} "
                             "— exact positive-root count.")
            else:
                verdict = ("satisfy the necessary sign condition (multistationarity "
                           "occurs for some total concentrations)" if at_rates else
                           "do NOT make any critical coefficient negative on their own")
                lines.append(f"  → The supplied rate constants {verdict}.")

    return MultistationarityResult(
        monostationary=monostationary,
        multistationary_possible=multistationary_possible,
        critical_function=phi,
        coefficient_terms=coeff_terms,
        region_conditions=region,
        multistationary_at_rates=at_rates,
        summary_lines=lines,
        region_boundary=region_boundary,
        steady_state_polynomial=ss_poly,
    )
