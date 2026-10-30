"""
Validation case: Goldbeter-Koshland covalent modification switch.

Reference: Goldbeter A, Koshland DE Jr. (1981) An amplified sensitivity arising
from covalent modification in biological systems. PNAS 78(11): 6840–6844.

Mechanism (full mass-action, no QSS approximation):
  W + E1 ⇌ WE1  →  Wp + E1     (kinase arm)
  Wp + E2 ⇌ WpE2 →  W + E2     (phosphatase arm)

Key CRNT result: δ=1, NOT weakly reversible, D1T applies
  → at most one steady state per stoichiometry class
  → bistability is structurally impossible for ANY positive rate constants

This is consistent with the original GK result: the switch produces sigmoidal
(ultrasensitive) but MONOSTABLE dose-response — no hysteresis, no bistability.
"""
import numpy as np
import pytest
import sympy
from mantis import CRNetwork


GK_STRINGS = [
    "W + E1 <-> WE1",
    "WE1 -> Wp + E1",
    "Wp + E2 <-> WpE2",
    "WpE2 -> W + E2",
]

# Simple (non-stiff) parameters for fast tests — ultrasensitivity not required here
GK_RATES = {
    "E1 + W -> WE1":   1.0,    # k_on  (kinase)
    "WE1 -> E1 + W":   0.5,    # k_off (kinase)
    "WE1 -> E1 + Wp":  1.0,    # kcat  (phosphorylation)
    "E2 + Wp -> WpE2": 1.0,    # k_on  (phosphatase)
    "WpE2 -> E2 + Wp": 0.5,    # k_off (phosphatase)
    "WpE2 -> E2 + W":  1.0,    # kcat  (dephosphorylation)
}

# Equal enzyme amounts → symmetric kinase/phosphatase activity
GK_IC = {
    "W":    1.0,
    "E1":   0.1,
    "WE1":  0.0,
    "Wp":   0.0,
    "E2":   0.1,
    "WpE2": 0.0,
}


@pytest.fixture(scope="module")
def rn():
    return CRNetwork.from_string(GK_STRINGS, rates=GK_RATES)


# ── CRNT structural analysis ───────────────────────────────────────────────────

def test_n_species(rn):
    assert rn.n_species == 6

def test_n_complexes(rn):
    assert rn.n_complexes == 6

def test_n_linkage_classes(rn):
    # LC1: {W+E1, WE1, Wp+E1}  LC2: {Wp+E2, WpE2, W+E2}
    assert rn.n_linkage_classes == 2

def test_stoichiometry_rank(rn):
    from mantis.stoichiometry import matrix_rank
    assert matrix_rank(rn.stoichiometry_matrix) == 3

def test_deficiency(rn):
    # δ = 6 - 2 - 3 = 1
    assert rn.deficiency == 1

def test_not_weakly_reversible(rn):
    # WE1 → Wp+E1 and WpE2 → W+E2 are irreversible — no path back in either LC
    assert rn.is_weakly_reversible is False

def test_dzt_not_applicable(rn):
    assert rn._crnt_result.deficiency_zero_applies is False

def test_d1t_applicable(rn):
    # D1T structural conditions: δ=1, per-LC deficiency ≤ 1 (both = 0 here)
    assert rn._crnt_result.deficiency_one_applies is True


# ── Conservation laws ──────────────────────────────────────────────────────────

def test_n_conservation_laws(rn):
    # dim(left null space of N) = n_species − rank = 6 − 3 = 3
    assert len(rn.conservation_laws) == 3

def test_conservation_law_e1(rn):
    """E1 + WE1 = const (kinase is neither created nor destroyed)."""
    E1   = sympy.Symbol("E1")
    WE1  = sympy.Symbol("WE1")
    e1_law = [e for e in rn.conservation_laws if e.has(E1) and e.has(WE1)]
    assert len(e1_law) == 1

def test_conservation_law_e2(rn):
    """E2 + WpE2 = const (phosphatase is neither created nor destroyed)."""
    E2   = sympy.Symbol("E2")
    WpE2 = sympy.Symbol("WpE2")
    e2_law = [e for e in rn.conservation_laws if e.has(E2) and e.has(WpE2)]
    assert len(e2_law) == 1

def test_conservation_law_substrate(rn):
    """W + WE1 + Wp + WpE2 = const (total substrate conserved)."""
    W    = sympy.Symbol("W")
    WE1  = sympy.Symbol("WE1")
    Wp   = sympy.Symbol("Wp")
    WpE2 = sympy.Symbol("WpE2")
    substrate_law = [
        e for e in rn.conservation_laws
        if e.has(W) and e.has(WE1) and e.has(Wp) and e.has(WpE2)
    ]
    assert len(substrate_law) == 1


# ── Numerical steady state ─────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def ss(rn):
    ss_list = rn.steady_states(GK_IC, n_attempts=20, seed=0)
    assert len(ss_list) >= 1, "Solver found no steady state"
    return ss_list[0]


def test_ss_nonneg(ss):
    for sp, c in ss.concentrations.items():
        assert c >= -1e-12, f"{sp} = {c}"

def test_ss_stable(ss):
    assert ss.is_stable, f"Eigenvalues: {ss.eigenvalues}"

def test_ss_residual(ss):
    assert ss.residual < 1e-6, f"Residual: {ss.residual}"

def test_ss_e1_conservation(ss):
    """E1 + WE1 = E1_total at steady state."""
    c = ss.concentrations
    total = c["E1"] + c["WE1"]
    assert abs(total - GK_IC["E1"]) / GK_IC["E1"] < 0.01

def test_ss_e2_conservation(ss):
    """E2 + WpE2 = E2_total at steady state."""
    c = ss.concentrations
    total = c["E2"] + c["WpE2"]
    assert abs(total - GK_IC["E2"]) / GK_IC["E2"] < 0.01

def test_ss_substrate_conservation(ss):
    """W + WE1 + Wp + WpE2 = W_total at steady state."""
    c = ss.concentrations
    total = c["W"] + c["WE1"] + c["Wp"] + c["WpE2"]
    Wt = GK_IC["W"]
    assert abs(total - Wt) / Wt < 0.01

def test_ss_symmetry(ss):
    """At equal kinase/phosphatase activity, total phosphorylated = total
    unphosphorylated (= Wt/2 exactly, by symmetry of the network and ICs)."""
    c = ss.concentrations
    phosphorylated   = c["Wp"]  + c["WpE2"]
    unphosphorylated = c["W"]   + c["WE1"]
    Wt = GK_IC["W"]
    # Both should equal Wt/2 = 0.5 within 2%
    assert abs(phosphorylated / Wt - 0.5) < 0.02, (
        f"Phosphorylated fraction = {phosphorylated/Wt:.4f}, expected 0.5"
    )
    assert abs(unphosphorylated / Wt - 0.5) < 0.02


# ── Monostability scan ─────────────────────────────────────────────────────────

def test_single_ss_at_multiple_kinase_rates(rn):
    """D1T guarantees at most one SS per stoichiometry class.
    Scanning kcat1 over a 10x range should yield exactly one SS each time."""
    for kcat1 in [0.2, 0.5, 1.0, 2.0, 5.0]:
        rates = dict(GK_RATES)
        rates["WE1 -> E1 + Wp"] = kcat1
        rn_scan = CRNetwork.from_string(GK_STRINGS, rates=rates)
        ss_list = rn_scan.steady_states(GK_IC, n_attempts=15, seed=42)
        assert len(ss_list) == 1, (
            f"Found {len(ss_list)} steady states at kcat1={kcat1}; expected 1 (D1T)"
        )
