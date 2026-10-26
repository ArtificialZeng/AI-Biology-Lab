import numpy as np
import pytest
from mantis import CRNetwork


# ── 1. Simple reversible A ↔ B ────────────────────────────────────────────────

def test_ab_steady_state():
    """
    A <-> B with k_f=1, k_r=0.5.
    Analytical SS: A* = k_r/(k_f+k_r) * total, B* = k_f/(k_f+k_r) * total.
    Starting from [A]=2, [B]=0 → total=2 → A*=2/3, B*=4/3.
    """
    rn = CRNetwork.from_string(
        ["A <-> B"],
        rates={"A -> B": 1.0, "B -> A": 0.5},
    )
    ic = {"A": 2.0, "B": 0.0}
    ss_list = rn.steady_states(ic, n_attempts=20, seed=0)
    assert len(ss_list) >= 1, "Should find at least one steady state"
    ss = ss_list[0]
    A_star = ss.concentrations["A"]
    B_star = ss.concentrations["B"]
    # Analytical: total=2, A*=2*(0.5/1.5)=2/3, B*=2*(1.0/1.5)=4/3
    assert abs(A_star - 2.0 / 3.0) < 1e-5, f"A* = {A_star}, expected {2/3}"
    assert abs(B_star - 4.0 / 3.0) < 1e-5, f"B* = {B_star}, expected {4/3}"
    assert ss.is_stable


def test_ab_conservation_at_ss():
    """A + B = const at steady state."""
    rn = CRNetwork.from_string(
        ["A <-> B"],
        rates={"A -> B": 1.0, "B -> A": 0.5},
    )
    ic = {"A": 3.0, "B": 1.0}
    ss_list = rn.steady_states(ic, n_attempts=20, seed=0)
    ss = ss_list[0]
    total = ss.concentrations["A"] + ss.concentrations["B"]
    assert abs(total - 4.0) < 1e-5, f"A+B = {total}, expected 4"


def test_ab_eigenvalues_stable():
    """Stable SS eigenvalues for A<->B both have Re < 0."""
    rn = CRNetwork.from_string(
        ["A <-> B"],
        rates={"A -> B": 1.0, "B -> A": 0.5},
    )
    ic = {"A": 1.0, "B": 1.0}
    ss_list = rn.steady_states(ic, n_attempts=20, seed=0)
    ss = ss_list[0]
    assert ss.is_stable


# ── 2. Michaelis-Menten ────────────────────────────────────────────────────────

def test_mm_conservation():
    """E + ES = E_total at steady state."""
    rn = CRNetwork.from_string(
        ["E + S <-> ES", "ES -> E + P"],
        rates={
            "E + S -> ES": 1e6,
            "ES -> E + S": 1e3,
            "ES -> E + P": 100.0,
        },
    )
    E0 = 1e-6
    S0 = 1e-4
    ic = {"E": E0, "S": S0, "ES": 0.0, "P": 0.0}
    ss_list = rn.steady_states(ic, n_attempts=30, seed=0)
    assert len(ss_list) >= 1
    ss = ss_list[0]
    E_total = ss.concentrations["E"] + ss.concentrations["ES"]
    assert abs(E_total - E0) / E0 < 0.01, f"E_total = {E_total}, expected {E0}"


# ── 3. CHA system ──────────────────────────────────────────────────────────────

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../hairpin_design_scripts'))

CHA_STRINGS = [
    "miR21 + H1 <-> miR21_H1",
    "miR21_H1 + H2 <-> H1H2 + miR21",
    "H1H2 + CP <-> H1H2_CP",
    "H1 + H2 <-> H1H2",
]

CHA_RATES_10NM = {
    "miR21 + H1 -> miR21_H1": 3.0e5,
    "miR21_H1 -> miR21 + H1": 2.687e-3,
    "miR21_H1 + H2 -> H1H2 + miR21": 1.0e6,
    "H1H2 + miR21 -> miR21_H1 + H2": 3.383e-3,
    "H1H2 + CP -> H1H2_CP": 2.0e5,
    "H1H2_CP -> H1H2 + CP": 5.79,
    "H1 + H2 -> H1H2": 2.379,
    "H1H2 -> H1 + H2": 7.208e-17,
}

CHA_IC_10NM = {
    "H1": 100e-9,
    "H2": 100e-9,
    "CP": 100e-9,
    "miR21": 10e-9,
    "miR21_H1": 0.0,
    "H1H2": 0.0,
    "H1H2_CP": 0.0,
}


def test_cha_ss_nonneg():
    """CHA steady-state concentrations are all non-negative."""
    rn = CRNetwork.from_string(CHA_STRINGS, rates=CHA_RATES_10NM)
    ss_list = rn.steady_states(CHA_IC_10NM, n_attempts=30, seed=0)
    assert len(ss_list) >= 1
    ss = ss_list[0]
    for s, c in ss.concentrations.items():
        assert c >= -1e-15, f"{s} concentration is negative: {c}"


def test_cha_ss_conservation():
    """CHA conservation laws are satisfied at steady state."""
    rn = CRNetwork.from_string(CHA_STRINGS, rates=CHA_RATES_10NM)
    ss_list = rn.steady_states(CHA_IC_10NM, n_attempts=30, seed=0)
    ss = ss_list[0]
    c = ss.concentrations
    # miR21 total: miR21 + miR21_H1 = 10 nM
    total_mir21 = c["miR21"] + c["miR21_H1"]
    assert abs(total_mir21 - 10e-9) / 10e-9 < 0.01, f"miR21 total = {total_mir21}"
    # CP total: CP + H1H2_CP = 100 nM
    total_cp = c["CP"] + c["H1H2_CP"]
    assert abs(total_cp - 100e-9) / 100e-9 < 0.01, f"CP total = {total_cp}"


def test_cha_ss_signal():
    """With 10 nM miR21, significant H1H2_CP should form at steady state."""
    rn = CRNetwork.from_string(CHA_STRINGS, rates=CHA_RATES_10NM)
    ss_list = rn.steady_states(CHA_IC_10NM, n_attempts=30, seed=0)
    ss = ss_list[0]
    signal = ss.concentrations["H1H2_CP"]
    # At 10 nM miR21, kinetic_simulation.py gives ~0.34 nM at t=7200s
    assert signal > 0.1e-9, f"H1H2_CP = {signal:.2e}, expected > 0.1 nM"


def test_cha_ss_stable():
    """CHA steady state should be stable (all significant eigenvalues Re < 0)."""
    rn = CRNetwork.from_string(CHA_STRINGS, rates=CHA_RATES_10NM)
    ss_list = rn.steady_states(CHA_IC_10NM, n_attempts=30, seed=0)
    ss = ss_list[0]
    assert ss.is_stable, f"SS not stable. Eigenvalues: {ss.eigenvalues}"
