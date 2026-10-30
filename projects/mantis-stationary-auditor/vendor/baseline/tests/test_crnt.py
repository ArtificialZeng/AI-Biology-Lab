import numpy as np
import pytest
from mantis.parsing import parse_reactions
from mantis.stoichiometry import build_species_list, build_stoichiometry_matrix, matrix_rank
from mantis.crnt import crnt_analysis


def _analyze(strings):
    rxns = parse_reactions(strings)
    species = build_species_list(rxns)
    N = build_stoichiometry_matrix(rxns, species)
    return crnt_analysis(rxns, N, species)


# ── 1. Simple reversible A ↔ B ────────────────────────────────────────────────

def test_simple_reversible_deficiency():
    r = _analyze(["A <-> B"])
    assert r.n_species == 2
    assert r.n_complexes == 2
    assert r.n_linkage_classes == 1
    assert r.rank_N == 1
    assert r.deficiency == 0     # δ = 2 - 1 - 1 = 0


def test_simple_reversible_weak_rev():
    r = _analyze(["A <-> B"])
    assert r.is_weakly_reversible is True


def test_simple_reversible_dzt():
    r = _analyze(["A <-> B"])
    assert r.deficiency_zero_applies is True


def test_simple_reversible_conservation():
    from mantis.stoichiometry import conservation_laws_sympy
    rxns = parse_reactions(["A <-> B"])
    species = build_species_list(rxns)
    N = build_stoichiometry_matrix(rxns, species)
    laws = conservation_laws_sympy(N, species)
    assert len(laws) == 1
    import sympy
    A, B = sympy.Symbol("A"), sympy.Symbol("B")
    # Law should be A + B = const
    assert laws[0].equals(A + B) or laws[0].equals(B + A)


# ── 2. Michaelis-Menten (explicit, irreversible final step) ───────────────────

def test_michaelis_menten_deficiency():
    r = _analyze(["E + S <-> ES", "ES -> E + P"])
    # Complexes: {E+S}, {ES}, {E+P}  → n=3
    # Reactions: E+S→ES, ES→E+S, ES→E+P  → irreversible final step
    # All complexes are connected in one linkage class
    assert r.n_complexes == 3
    assert r.n_linkage_classes == 1
    assert r.deficiency == 0     # δ = 3 - 1 - 2 = 0


def test_michaelis_menten_not_weakly_reversible():
    r = _analyze(["E + S <-> ES", "ES -> E + P"])
    # ES -> E+P has no reverse path back
    assert r.is_weakly_reversible is False


def test_michaelis_menten_dzt_not_applicable():
    r = _analyze(["E + S <-> ES", "ES -> E + P"])
    assert r.deficiency_zero_applies is False


# ── 3. Brusselator ────────────────────────────────────────────────────────────
# When A, B, D, E are treated as full dynamic species (our default), CRNT gives
# n=7, l=3, rank=4 → δ = 7-3-4 = 0.
# The literature's δ=2 result assumes A, B, D, E are chemostatted (fixed-
# concentration parameters), which reduces the effective stoichiometry rank to 2.
# The key property for the contrast with CHA is: NOT weakly reversible → DZT
# does not apply regardless of δ, and numerical oscillations can exist.

_BRUSSELATOR = [
    "A -> X",
    "2X + Y -> 3X",
    "B + X -> Y + D",
    "X -> E",
]


def test_brusselator_deficiency():
    r = _analyze(_BRUSSELATOR)
    # Full-species formulation: n=7, l=3, rank=4 → δ=0
    assert r.n_complexes == 7
    assert r.n_linkage_classes == 3
    assert r.rank_N == 4
    assert r.deficiency == 0


def test_brusselator_not_weakly_reversible():
    r = _analyze(_BRUSSELATOR)
    assert r.is_weakly_reversible is False


def test_brusselator_dzt_not_applicable():
    r = _analyze(_BRUSSELATOR)
    # DZT doesn't apply: not weakly reversible
    assert r.deficiency_zero_applies is False


# ── 4. CHA cascade ────────────────────────────────────────────────────────────

CHA_REACTIONS = [
    "miR21 + H1 <-> miR21_H1",
    "miR21_H1 + H2 <-> H1H2 + miR21",
    "H1H2 + CP <-> H1H2_CP",
    "H1 + H2 <-> H1H2",
]


def test_cha_complexes():
    r = _analyze(CHA_REACTIONS)
    assert r.n_species == 7
    assert r.n_complexes == 8


def test_cha_linkage_classes():
    r = _analyze(CHA_REACTIONS)
    assert r.n_linkage_classes == 4


def test_cha_rank_and_deficiency():
    r = _analyze(CHA_REACTIONS)
    assert r.rank_N == 3
    assert r.deficiency == 1    # δ = 8 - 4 - 3 = 1


def test_cha_weakly_reversible():
    r = _analyze(CHA_REACTIONS)
    assert r.is_weakly_reversible is True


def test_cha_dzt_not_applicable():
    r = _analyze(CHA_REACTIONS)
    assert r.deficiency_zero_applies is False


def test_cha_d1t_applicable():
    r = _analyze(CHA_REACTIONS)
    assert r.deficiency_one_applies is True


def test_cha_conservation_laws():
    from mantis.stoichiometry import conservation_laws_sympy
    import sympy
    rxns = parse_reactions(CHA_REACTIONS)
    species = build_species_list(rxns)
    N = build_stoichiometry_matrix(rxns, species)
    laws = conservation_laws_sympy(N, species)
    assert len(laws) == 4
    # Each law should be a SymPy expression; collect species that appear
    all_syms = set()
    for law in laws:
        all_syms |= {str(s) for s in law.free_symbols}
    # All 7 CHA species must appear somewhere across the 4 laws
    expected = {"CP", "H1", "H1H2", "H1H2_CP", "H2", "miR21", "miR21_H1"}
    assert all_syms == expected
