"""SymPy ODE generation, Jacobian, and substitution helpers."""
from __future__ import annotations

import sympy

from .parsing import Reaction
from .stoichiometry import build_stoichiometry_matrix


def make_species_symbols(species: list[str]) -> dict[str, sympy.Symbol]:
    """Create a SymPy Symbol for each species name."""
    return {s: sympy.Symbol(s) for s in species}


def make_rate_symbols(
    reactions: list[Reaction],
) -> tuple[dict[str, sympy.Symbol], dict[str, sympy.Symbol]]:
    """
    Create rate symbols k_1, k_2, … for each reaction.
    Returns:
        rate_syms  : dict of symbol_name → Symbol  (e.g. 'k_1' → k_1)
        key_to_sym : dict of rate_key → Symbol
    """
    rate_syms: dict[str, sympy.Symbol] = {}
    key_to_sym: dict[str, sympy.Symbol] = {}
    for i, rxn in enumerate(reactions, 1):
        name = f"k_{i}"
        sym = sympy.Symbol(name)
        rate_syms[name] = sym
        key_to_sym[rxn.rate_key] = sym
    return rate_syms, key_to_sym


def mass_action_flux(
    rxn: Reaction,
    species_syms: dict[str, sympy.Symbol],
    rate_sym: sympy.Symbol,
    chemostatted_values: dict[str, float] | None = None,
) -> sympy.Expr:
    """Return rate_sym * Π(species^coeff) for each reactant.

    Chemostatted species are substituted as Float constants rather than symbols.
    """
    chem = chemostatted_values or {}
    expr = rate_sym
    for name, coeff in rxn.reactants:
        if name in chem:
            expr = expr * sympy.Float(chem[name]) ** coeff
        else:
            expr = expr * species_syms[name] ** coeff
    return expr


def build_odes(
    reactions: list[Reaction],
    species: list[str],
    species_syms: dict[str, sympy.Symbol],
    rate_syms: dict[str, sympy.Symbol],
    key_to_sym: dict[str, sympy.Symbol],
    chemostatted_values: dict[str, float] | None = None,
) -> dict[str, sympy.Expr]:
    """
    Return dict: species_name → d[species]/dt as SymPy expression.
    Uses stoichiometry matrix: d[Xi]/dt = Σⱼ N[i,j] * flux_j(y).

    Chemostatted species appear as Float constants in the flux expressions.
    """
    N = build_stoichiometry_matrix(reactions, species)
    fluxes = [
        mass_action_flux(rxn, species_syms, key_to_sym[rxn.rate_key], chemostatted_values)
        for rxn in reactions
    ]
    odes: dict[str, sympy.Expr] = {}
    for i, sp in enumerate(species):
        expr = sympy.Integer(0)
        for j, flux in enumerate(fluxes):
            coeff = int(round(N[i, j]))
            if coeff != 0:
                expr += coeff * flux
        odes[sp] = sympy.expand(expr)
    return odes


def build_jacobian(
    odes: dict[str, sympy.Expr],
    species: list[str],
    species_syms: dict[str, sympy.Symbol],
) -> sympy.Matrix:
    """Return the symbolic Jacobian matrix (n_species × n_species)."""
    return sympy.Matrix([
        [sympy.diff(odes[si], species_syms[sj]) for sj in species]
        for si in species
    ])


def substitute_rates(
    expr,
    key_to_sym: dict[str, sympy.Symbol],
    rate_values: dict[str, float],
) -> sympy.Expr:
    """Substitute numerical rate constant values into a symbolic expression."""
    subs = {}
    for key, sym in key_to_sym.items():
        if key in rate_values:
            subs[sym] = sympy.Float(rate_values[key])
    if isinstance(expr, sympy.Matrix):
        return expr.subs(subs)
    return expr.subs(subs)


def substitute_steady_state(
    J: sympy.Matrix,
    species_syms: dict[str, sympy.Symbol],
    ss: dict[str, float],
) -> sympy.Matrix:
    """Substitute steady-state concentrations into Jacobian."""
    subs = {species_syms[s]: sympy.Float(v) for s, v in ss.items()}
    return J.subs(subs)
