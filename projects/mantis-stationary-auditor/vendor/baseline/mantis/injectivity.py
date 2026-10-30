"""Injectivity test — Craciun & Feinberg (*SIAM J. Appl. Math.* 2005, 2006).

A mass-action network is **injective** when its species-formation-rate map is
injective on each stoichiometric compatibility class; an injective network has *at
most one* positive steady state per class, so it **cannot be multistationary**.
Crucially this rules out multistationarity for *higher-deficiency* networks too,
where the Deficiency One Theorem is silent.

Craciun & Feinberg reduce injectivity to a sign condition on a determinant.  Form
the square "steady-state map"

    F(c) = ( rank(N) independent components of N·v(c) ,  the conservation totals ),

whose Jacobian ``J(c)`` is an n×n matrix whose entries are polynomials in the
(positive) concentrations ``c`` and the (positive) rate constants ``k``.  If
``det J(c)`` is **sign-definite** — every monomial coefficient shares one sign, so
the determinant never vanishes for any positive ``c`` and ``k`` — then ``F`` is
injective and the network admits at most one positive steady state per class.

The test is *sufficient*: a sign-definite determinant certifies injectivity, while a
sign-indefinite one is inconclusive (the network may still be injective for specific
rates, or may genuinely be multistationary, as for the Schlögl model).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .parsing import Reaction
from .stoichiometry import build_stoichiometry_matrix, conservation_laws_sympy, matrix_rank


@dataclass
class InjectivityResult:
    """Outcome of the Craciun–Feinberg injectivity test.

    Attributes
    ----------
    injective : bool
        True ⇒ the network is certified injective (at most one positive steady state
        per compatibility class). False ⇒ inconclusive (not sign-definite).
    is_sign_definite : bool
        Whether the steady-state Jacobian determinant is sign-definite.
    determinant : sympy.Expr
        The symbolic determinant (in concentrations and rate constants).
    reason : str
    summary_lines : list[str]
    """
    injective: bool
    is_sign_definite: bool
    determinant: "object" = None
    reason: str = ""
    summary_lines: list[str] = field(default_factory=list)

    def __str__(self) -> str:
        return "\n".join(self.summary_lines)


def test_injectivity(
    reactions: list[Reaction],
    species: list[str],
    chemostatted_values: dict[str, float] | None = None,
) -> InjectivityResult:
    """Run the determinant sign-definiteness injectivity test (symbolic, rate-free)."""
    import sympy

    from .symbolic import make_species_symbols, make_rate_symbols, mass_action_flux

    chem = set((chemostatted_values or {}).keys())

    # Symbols for every species appearing in the reactions (dynamic + chemostatted),
    # so chemostatted species enter the determinant as positive parameters.
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
    if n == 0:
        return InjectivityResult(True, True, sympy.Integer(0), "no dynamic species",
                                 ["Injectivity: trivially injective (no dynamic species)."])

    # Symbolic mass-action fluxes (all species symbolic — nothing folded to a number).
    fluxes = [mass_action_flux(rxn, sp_syms, key_to_sym[rxn.rate_key], None)
              for rxn in reactions]

    N = build_stoichiometry_matrix(reactions, species)
    rank = matrix_rank(N)

    # rank(N) independent rows of N·v, chosen by QR column pivoting on Nᵀ.
    if rank > 0:
        from scipy.linalg import qr
        _, _, perm = qr(N.T, pivoting=True)
        ind = list(perm[:rank])
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

    # Conservation-law rows (their gradients pin down the remaining n - rank directions).
    laws = conservation_laws_sympy(N, species)
    for law in laws:
        rows.append(law)

    if len(rows) != n:
        return InjectivityResult(
            False, False, None,
            "could not form a square steady-state map",
            ["Injectivity: inconclusive (degenerate stoichiometry)."],
        )

    J = sympy.Matrix([[sympy.diff(r, c) for c in dyn_syms] for r in rows])
    det = sympy.expand(J.det())

    if det == 0:
        return InjectivityResult(
            False, False, det,
            "determinant is identically zero",
            [
                "Injectivity (Craciun–Feinberg): NOT certified",
                "  → The steady-state Jacobian determinant vanishes identically "
                "(degenerate); the test is inconclusive.",
            ],
        )

    poly = sympy.Poly(det, *(dyn_syms + list(rate_syms.values())
                             + [sp_syms[s] for s in sorted(all_names) if s in chem]))
    signs = {int(np.sign(float(c))) for c in poly.coeffs() if c != 0}
    sign_definite = len(signs) == 1

    if sign_definite:
        lines = [
            "Injectivity (Craciun–Feinberg): CERTIFIED injective",
            "  → The steady-state Jacobian determinant is sign-definite, so it never "
            "vanishes for positive concentrations and rate constants.",
            "  → The network admits AT MOST ONE positive steady state per stoichiometric",
            "    compatibility class — multistationarity is structurally impossible",
            "    (for all rate constants).",
        ]
    else:
        lines = [
            "Injectivity (Craciun–Feinberg): NOT certified (inconclusive)",
            "  → The determinant is sign-indefinite; this sufficient test cannot rule out",
            "    multiple steady states (the network may be multistationary for some rates).",
        ]

    return InjectivityResult(
        injective=sign_definite,
        is_sign_definite=sign_definite,
        determinant=det,
        reason="sign-definite determinant" if sign_definite else "sign-indefinite determinant",
        summary_lines=lines,
    )
