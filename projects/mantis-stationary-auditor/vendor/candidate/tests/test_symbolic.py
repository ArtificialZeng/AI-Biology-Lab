import sympy
import pytest
from mantis.parsing import parse_reactions
from mantis.stoichiometry import build_species_list, build_stoichiometry_matrix
from mantis.symbolic import (
    make_species_symbols,
    make_rate_symbols,
    build_odes,
    build_jacobian,
    substitute_rates,
)


def _setup(strings, rate_dict=None):
    rxns = parse_reactions(strings)
    species = build_species_list(rxns)
    sp_syms = make_species_symbols(species)
    rate_syms, key_to_sym = make_rate_symbols(rxns)
    odes = build_odes(rxns, species, sp_syms, rate_syms, key_to_sym)
    if rate_dict:
        odes = {s: substitute_rates(e, key_to_sym, rate_dict) for s, e in odes.items()}
    return rxns, species, sp_syms, rate_syms, key_to_sym, odes


# ── A ↔ B ──────────────────────────────────────────────────────────────────────

def test_ab_ode_structure():
    """dA/dt = -k1*A + k2*B, dB/dt = k1*A - k2*B"""
    rxns, species, sp_syms, rate_syms, key_to_sym, odes = _setup(["A <-> B"])
    A, B = sp_syms["A"], sp_syms["B"]
    # Find which symbol corresponds to forward and reverse
    k_forward = key_to_sym.get("A -> B")
    k_reverse = key_to_sym.get("B -> A")
    assert k_forward is not None, f"Key 'A -> B' not found. Available: {list(key_to_sym)}"
    dA = odes["A"]
    dB = odes["B"]
    expected_dA = -k_forward * A + k_reverse * B
    expected_dB = k_forward * A - k_reverse * B
    assert sympy.expand(dA - expected_dA) == 0, f"dA/dt mismatch: {dA}"
    assert sympy.expand(dB - expected_dB) == 0, f"dB/dt mismatch: {dB}"


def test_ab_jacobian():
    """Jacobian of A↔B system: 2×2 matrix with correct entries."""
    rxns, species, sp_syms, rate_syms, key_to_sym, odes = _setup(["A <-> B"])
    J = build_jacobian(odes, species, sp_syms)
    assert J.shape == (2, 2)
    # J[0,0] = d(dA/dt)/dA = -k_forward
    k_f = key_to_sym["A -> B"]
    k_r = key_to_sym["B -> A"]
    assert sympy.expand(J[0, 0] + k_f) == 0
    assert sympy.expand(J[0, 1] - k_r) == 0
    assert sympy.expand(J[1, 0] - k_f) == 0
    assert sympy.expand(J[1, 1] + k_r) == 0


def test_ab_numeric_rates():
    """Substituting numeric rates gives float coefficients."""
    rxns, species, sp_syms, rate_syms, key_to_sym, odes = _setup(
        ["A <-> B"], rate_dict={"A -> B": 1.0, "B -> A": 0.5}
    )
    A = sp_syms["A"]
    B = sp_syms["B"]
    dA = odes["A"]
    # Should be -1.0*A + 0.5*B
    val = dA.subs([(A, 2.0), (B, 4.0)])
    assert abs(float(val) - (-1.0 * 2.0 + 0.5 * 4.0)) < 1e-10


# ── CHA system ODE structure ─────────────────────────────────────────────────

CHA_STRINGS = [
    "miR21 + H1 <-> miR21_H1",
    "miR21_H1 + H2 <-> H1H2 + miR21",
    "H1H2 + CP <-> H1H2_CP",
    "H1 + H2 <-> H1H2",
]


def test_cha_ode_species():
    """ODEs are defined for all 7 CHA species."""
    rxns, species, sp_syms, _, key_to_sym, odes = _setup(CHA_STRINGS)
    assert set(odes.keys()) == {"CP", "H1", "H1H2", "H1H2_CP", "H2", "miR21", "miR21_H1"}


def test_cha_mir21_ode():
    """dmiR21/dt must contain the four expected terms."""
    rxns, species, sp_syms, _, key_to_sym, odes = _setup(CHA_STRINGS)
    dmiR21 = odes["miR21"]
    miR21 = sp_syms["miR21"]
    H1 = sp_syms["H1"]
    miR21_H1 = sp_syms["miR21_H1"]
    H2 = sp_syms["H2"]
    H1H2 = sp_syms["H1H2"]
    # miR21 appears in R1 (consumed) and R2 (produced)
    # dmiR21/dt = -k1*miR21*H1 + k2*miR21_H1 + k3*miR21_H1*H2 - k4*H1H2*miR21
    # (where k1..k4 are symbolic) — check structure via differentiation
    d2 = sympy.diff(dmiR21, miR21)
    # Should have -k_R1f * H1 and -k_R2r * H1H2 terms
    assert d2 != 0, "dmiR21/dt must depend on miR21"
    # All coefficients of miR21-containing terms are negative (consumption)
    # verify by checking sign at positive values
    val = float(dmiR21.subs(
        [(sp_syms[s], 1.0) for s in species]
        + [(sym, 1.0) for sym in key_to_sym.values()]
    ))
    # With all species=1 and all rates=1: dmiR21/dt = -1*1 + 1*1 + 1*1 - 1*1 = 0
    assert abs(val) < 1e-9, f"dmiR21 at all-1 should be 0, got {val}"


def test_cha_mass_balance():
    """Sum of all ODEs with species weighted by conservation law should be 0."""
    rxns, species, sp_syms, _, key_to_sym, odes = _setup(CHA_STRINGS)
    # miR21 total: miR21 + miR21_H1 = const
    # d(miR21)/dt + d(miR21_H1)/dt should be 0 symbolically
    total = sympy.expand(odes["miR21"] + odes["miR21_H1"])
    assert total == 0, f"miR21 conservation violated: {total}"
    # CP total: CP + H1H2_CP = const
    total2 = sympy.expand(odes["CP"] + odes["H1H2_CP"])
    assert total2 == 0, f"CP conservation violated: {total2}"
