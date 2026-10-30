"""Stoichiometry matrix, species list, matrix rank, and conservation-law computation."""
from math import gcd
from functools import reduce

import numpy as np
import sympy

from .parsing import Complex, Reaction


def build_species_list(reactions: list[Reaction], chemostatted: set[str] | None = None) -> list[str]:
    """Collect all species from all reactions, return sorted list.

    Parameters
    ----------
    reactions : list[Reaction]
    chemostatted : set[str], optional
        Species names to exclude from the returned list (they are held fixed).
    """
    chem = chemostatted or set()
    species: set[str] = set()
    for rxn in reactions:
        for name, _ in rxn.reactants:
            if name not in chem:
                species.add(name)
        for name, _ in rxn.products:
            if name not in chem:
                species.add(name)
    return sorted(species)


def build_stoichiometry_matrix(
    reactions: list[Reaction],
    species: list[str],
) -> np.ndarray:
    """
    Return N of shape (n_species, n_reactions).
    N[i,j] = sum(product coeff) - sum(reactant coeff) for species i in reaction j.

    Species that appear in reactions but are not in `species` (e.g. chemostatted
    species excluded from the dynamic list) are silently ignored — they have no
    stoichiometry row.
    """
    sp_idx = {s: i for i, s in enumerate(species)}
    N = np.zeros((len(species), len(reactions)), dtype=float)
    for j, rxn in enumerate(reactions):
        for name, coeff in rxn.reactants:
            if name in sp_idx:
                N[sp_idx[name], j] -= coeff
        for name, coeff in rxn.products:
            if name in sp_idx:
                N[sp_idx[name], j] += coeff
    return N


def fold_chemostatted_into_rates(
    reactions: list[Reaction],
    rate_values: dict[str, float],
    chemostatted_values: dict[str, float],
) -> dict[str, float]:
    """
    Pre-multiply fixed chemostatted concentrations into rate constants.
    For A chemostatted at c_A: k_eff = k * c_A^stoich so the ODE sees only dynamic species.
    """
    eff: dict[str, float] = {}
    for rxn in reactions:
        base = rate_values.get(rxn.rate_key, 0.0)
        for name, coeff in rxn.reactants:
            if name in chemostatted_values:
                base = base * (chemostatted_values[name] ** coeff)
        eff[rxn.rate_key] = base
    return eff


def matrix_rank(N: np.ndarray) -> int:
    return int(np.linalg.matrix_rank(N))


def _lcm(a: int, b: int) -> int:
    return a * b // gcd(a, b)


def _scale_to_integers(vec: list[sympy.Rational]) -> list[int]:
    """Scale a rational vector so all entries are integers with GCD=1."""
    if all(v == 0 for v in vec):
        return [0] * len(vec)
    denoms = [int(sympy.denom(v)) for v in vec if v != 0]
    lcm = reduce(_lcm, denoms)
    scaled = [int(v * lcm) for v in vec]
    common = reduce(gcd, [abs(x) for x in scaled if x != 0])
    return [x // common for x in scaled]


def _make_nonneg_basis(vecs: list[list[int]]) -> list[list[int]]:
    """
    Transform a list of integer null-space vectors to have all non-negative entries.
    Physical conservation laws for mass-balanced networks always admit a non-negative
    integer basis.  Uses two passes: (1) add/subtract multiples to eliminate negatives,
    (2) reduce by subtracting non-negative vectors to find the minimal form.
    """
    V = [list(v) for v in vecs]
    k = len(V)
    n = len(V[0]) if k else 0

    # Step 1: flip sign of any vector that has more negatives than positives
    for i in range(k):
        n_neg = sum(1 for x in V[i] if x < 0)
        n_pos = sum(1 for x in V[i] if x > 0)
        if n_neg > n_pos:
            V[i] = [-x for x in V[i]]

    # Step 2: add ±multiples of other vectors to eliminate remaining negative entries
    for _ in range(20 * k):
        changed = False
        for i in range(k):
            if not any(x < 0 for x in V[i]):
                continue
            n_neg_i = sum(1 for x in V[i] if x < 0)
            for j in range(k):
                if i == j:
                    continue
                for sign in (1, -1):
                    for c in (1, 2, 3):
                        candidate = [V[i][d] + sign * c * V[j][d] for d in range(n)]
                        n_neg_c = sum(1 for x in candidate if x < 0)
                        if n_neg_c < n_neg_i:
                            V[i] = candidate
                            n_neg_i = n_neg_c
                            changed = True
                            if n_neg_i == 0:
                                break
                    if n_neg_i == 0:
                        break
                if n_neg_i == 0:
                    break
        if not changed:
            break

    # Step 3: reduce — subtract non-negative multiples of other vectors to minimize
    # the number of non-zero entries (find the "most elementary" form).
    for _ in range(10 * k):
        changed = False
        for i in range(k):
            for j in range(k):
                if i == j:
                    continue
                # Can we subtract V[j] from V[i] and stay non-negative?
                if all(V[i][d] >= V[j][d] >= 0 for d in range(n)):
                    reduced = [V[i][d] - V[j][d] for d in range(n)]
                    if any(x > 0 for x in reduced):  # don't reduce to zero
                        if sum(x > 0 for x in reduced) < sum(x > 0 for x in V[i]):
                            V[i] = reduced
                            changed = True
        if not changed:
            break

    return V


def conservation_laws_sympy(
    N: np.ndarray,
    species: list[str],
) -> list[sympy.Expr]:
    """
    Return left null space of N as SymPy expressions (one per conservation law).
    Uses exact rational arithmetic via SymPy; post-processes the basis to have
    all non-negative integer coefficients (physical moiety conservation form).
    """
    # Convert to integer SymPy matrix — stoich matrices always have integer entries
    M = sympy.Matrix([[int(round(x)) for x in row] for row in N])
    # Left null space: vectors c such that c^T N = 0 ↔ N^T c = 0
    null_vecs = M.T.nullspace()
    raw = [_scale_to_integers(list(vec)) for vec in null_vecs]
    positive = _make_nonneg_basis(raw)
    syms = [sympy.Symbol(s) for s in species]
    laws = []
    for coeffs in positive:
        expr = sum(c * s for c, s in zip(coeffs, syms) if c != 0)
        laws.append(expr)
    return laws
