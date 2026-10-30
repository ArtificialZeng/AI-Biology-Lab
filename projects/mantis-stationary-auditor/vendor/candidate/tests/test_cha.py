"""Integration tests for the CHA miR-21 system via CRNetwork public API."""
import pytest
from mantis import CRNetwork

CHA_STRINGS = [
    "miR21 + H1 <-> miR21_H1",
    "miR21_H1 + H2 <-> H1H2 + miR21",
    "H1H2 + CP <-> H1H2_CP",
    "H1 + H2 <-> H1H2",
]

CHA_RATES = {
    "miR21 + H1 -> miR21_H1": 3.0e5,
    "miR21_H1 -> miR21 + H1": 2.687e-3,
    "miR21_H1 + H2 -> H1H2 + miR21": 1.0e6,
    "H1H2 + miR21 -> miR21_H1 + H2": 3.383e-3,
    "H1H2 + CP -> H1H2_CP": 2.0e5,
    "H1H2_CP -> H1H2 + CP": 5.79,
    "H1 + H2 -> H1H2": 2.379,
    "H1H2 -> H1 + H2": 7.208e-17,
}

CHA_IC = {
    "H1": 100e-9,
    "H2": 100e-9,
    "CP": 100e-9,
    "miR21": 10e-9,
    "miR21_H1": 0.0,
    "H1H2": 0.0,
    "H1H2_CP": 0.0,
}


@pytest.fixture(scope="module")
def rn():
    return CRNetwork.from_string(CHA_STRINGS, rates=CHA_RATES)


# ── Structural properties ──────────────────────────────────────────────────────

def test_n_species(rn):
    assert rn.n_species == 7


def test_n_complexes(rn):
    assert rn.n_complexes == 8


def test_n_reactions(rn):
    assert rn.n_reactions == 8


def test_n_linkage_classes(rn):
    assert rn.n_linkage_classes == 4


def test_stoichiometry_rank(rn):
    from mantis.stoichiometry import matrix_rank
    assert matrix_rank(rn.stoichiometry_matrix) == 3


def test_deficiency(rn):
    assert rn.deficiency == 1


def test_weakly_reversible(rn):
    assert rn.is_weakly_reversible is True


# ── Conservation laws ──────────────────────────────────────────────────────────

def test_n_conservation_laws(rn):
    assert len(rn.conservation_laws) == 4


def test_conservation_law_mir21(rn):
    import sympy
    cl_exprs = rn.conservation_laws
    mir21 = sympy.Symbol("miR21")
    mir21_h1 = sympy.Symbol("miR21_H1")
    # One law must contain both miR21 and miR21_H1
    mir21_law = [e for e in cl_exprs if e.has(mir21) and e.has(mir21_h1)]
    assert len(mir21_law) == 1, f"No miR21 conservation law found; laws = {cl_exprs}"


def test_conservation_law_cp(rn):
    import sympy
    cl_exprs = rn.conservation_laws
    cp = sympy.Symbol("CP")
    h1h2_cp = sympy.Symbol("H1H2_CP")
    cp_law = [e for e in cl_exprs if e.has(cp) and e.has(h1h2_cp)]
    assert len(cp_law) == 1, f"No CP conservation law found; laws = {cl_exprs}"


# ── CRNT summary ───────────────────────────────────────────────────────────────

def test_crnt_summary_contains_delta(rn):
    summary = rn.crnt_summary()
    assert "δ = 1" in summary or "δ=1" in summary


def test_crnt_summary_weakly_reversible(rn):
    summary = rn.crnt_summary()
    assert "Yes" in summary or "yes" in summary


def test_crnt_summary_d1t(rn):
    summary = rn.crnt_summary()
    assert "Deficiency One" in summary or "D1T" in summary


# ── Steady state ───────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def ss(rn):
    ss_list = rn.steady_states(CHA_IC, n_attempts=30, seed=0)
    assert len(ss_list) >= 1
    return ss_list[0]


def test_ss_nonneg(ss):
    for sp, c in ss.concentrations.items():
        assert c >= -1e-15, f"{sp} = {c}"


def test_ss_stable(ss):
    assert ss.is_stable, f"Eigenvalues: {ss.eigenvalues}"


def test_ss_mir21_conservation(ss):
    c = ss.concentrations
    total = c["miR21"] + c["miR21_H1"]
    assert abs(total - 10e-9) / 10e-9 < 0.01


def test_ss_cp_conservation(ss):
    c = ss.concentrations
    total = c["CP"] + c["H1H2_CP"]
    assert abs(total - 100e-9) / 100e-9 < 0.01


def test_ss_signal(ss):
    # kinetic_simulation.py gives ~0.34 nM at t=7200s with these ICs
    assert ss.concentrations["H1H2_CP"] > 0.1e-9
