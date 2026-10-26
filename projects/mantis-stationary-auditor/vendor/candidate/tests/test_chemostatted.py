"""Tests for chemostatted (fixed-concentration) species support."""
import pytest
from mantis import CRNetwork


BRUSSELATOR_REACTIONS = [
    "A -> X",
    "2X + Y -> 3X",
    "B + X -> Y + D",
    "X -> E",
]

BRUSSELATOR_RATES = {
    "A -> X": 1.0,
    "2X + Y -> 3X": 1.0,
    "B + X -> Y + D": 1.0,
    "X -> E": 1.0,
}

# A=1, B=3: Hopf condition B > 1+A^2 = 2 is satisfied → unstable oscillatory fixed point
BRUSSELATOR_CHEM = CRNetwork.from_string(
    BRUSSELATOR_REACTIONS,
    rates=BRUSSELATOR_RATES,
    chemostatted={"A": 1.0, "B": 3.0, "D": 0.0, "E": 0.0},
)


# ── Species and property tests ─────────────────────────────────────────────────

def test_chemostatted_species_excluded():
    """Only dynamic species X and Y appear in rn.species."""
    assert BRUSSELATOR_CHEM.species == ["X", "Y"]


def test_chemostatted_property_returns_dict():
    """rn.chemostatted is a dict with all four fixed species."""
    chem = BRUSSELATOR_CHEM.chemostatted
    assert isinstance(chem, dict)
    assert set(chem.keys()) == {"A", "B", "D", "E"}
    assert chem["A"] == 1.0
    assert chem["B"] == 3.0


def test_chemostatted_property_is_copy():
    """Mutating the returned dict does not affect the network."""
    chem = BRUSSELATOR_CHEM.chemostatted
    chem["A"] = 999.0
    assert BRUSSELATOR_CHEM.chemostatted["A"] == 1.0


# ── ODE tests ─────────────────────────────────────────────────────────────────

def test_chemostatted_odes_have_only_dynamic():
    """ODEs are only produced for dynamic species X and Y."""
    odes = BRUSSELATOR_CHEM.odes()
    assert set(odes.keys()) == {"X", "Y"}


def test_chemostatted_odes_no_dynamic_chemostatted_symbols():
    """Chemostatted species A, B, D, E do not appear as free symbols in the ODEs."""
    import sympy
    odes = BRUSSELATOR_CHEM.odes()
    for sp, expr in odes.items():
        free = {str(s) for s in expr.free_symbols}
        for chem_sp in ("A", "B", "D", "E"):
            assert chem_sp not in free, (
                f"Chemostatted species {chem_sp} appears in d[{sp}]/dt: {expr}"
            )


def test_chemostatted_ode_values():
    """Verify correct dX/dt and dY/dt expressions for the Brusselator.

    Classical Brusselator with chemostatted A=a, B=b:
      dX/dt = a - (b+1)*X + X^2*Y
      dY/dt = b*X - X^2*Y
    With a=1, b=3 and rates k=1:
      dX/dt = 1.0 - 4*X + X^2*Y
      dY/dt = 3*X - X^2*Y
    """
    import sympy
    odes = BRUSSELATOR_CHEM.odes()
    X = sympy.Symbol("X")
    Y = sympy.Symbol("Y")
    # Test numerical values at X=1, Y=2
    val_X = float(odes["X"].subs([(X, 1.0), (Y, 2.0)]))
    val_Y = float(odes["Y"].subs([(X, 1.0), (Y, 2.0)]))
    # dX/dt = 1 - 4*1 + 1*2 = -1
    assert abs(val_X - (-1.0)) < 1e-10, f"dX/dt at (1,2) = {val_X}, expected -1"
    # dY/dt = 3*1 - 1*2 = 1
    assert abs(val_Y - 1.0) < 1e-10, f"dY/dt at (1,2) = {val_Y}, expected 1"


# ── Steady-state tests ─────────────────────────────────────────────────────────

def test_chemostatted_steady_state_unstable():
    """A=1, B=3: B > 1+A^2=2 → unstable fixed point at (X*=A=1, Y*=B/A=3)."""
    ss_list = BRUSSELATOR_CHEM.steady_states({"X": 1.5, "Y": 2.0}, seed=0)
    assert len(ss_list) >= 1, "Should find at least one steady state"
    ss = ss_list[0]
    assert not ss.is_stable, (
        f"Steady state should be unstable for B=3 > 1+A^2=2. "
        f"Eigenvalues: {ss.eigenvalues}"
    )
    assert ss.is_oscillatory, (
        f"Unstable fixed point should be oscillatory (complex eigenvalues with Re>0). "
        f"Eigenvalues: {ss.eigenvalues}"
    )


def test_chemostatted_steady_state_stable():
    """A=1, B=1.5: B < 1+A^2=2 → stable spiral at (X*=1, Y*=1.5)."""
    rn_stable = CRNetwork.from_string(
        BRUSSELATOR_REACTIONS,
        rates=BRUSSELATOR_RATES,
        chemostatted={"A": 1.0, "B": 1.5, "D": 0.0, "E": 0.0},
    )
    ss_list = rn_stable.steady_states({"X": 1.5, "Y": 2.0}, seed=0)
    assert len(ss_list) >= 1, "Should find at least one steady state"
    ss = ss_list[0]
    assert ss.is_stable, (
        f"Steady state should be stable for B=1.5 < 1+A^2=2. "
        f"Eigenvalues: {ss.eigenvalues}"
    )
    assert ss.is_oscillatory, (
        f"Stable spiral should be oscillatory (complex eigenvalues with Re<0). "
        f"Eigenvalues: {ss.eigenvalues}"
    )


def test_chemostatted_fixed_point_location():
    """The chemostatted Brusselator fixed point is X*=A=1, Y*=B/A=3."""
    ss_list = BRUSSELATOR_CHEM.steady_states({"X": 1.0, "Y": 3.0}, seed=0)
    assert len(ss_list) >= 1
    ss = ss_list[0]
    # Fixed point at (1, 3) for A=1, B=3
    assert abs(ss.concentrations["X"] - 1.0) < 0.01, (
        f"X* = {ss.concentrations['X']}, expected ~1.0"
    )
    assert abs(ss.concentrations["Y"] - 3.0) < 0.01, (
        f"Y* = {ss.concentrations['Y']}, expected ~3.0"
    )


# ── CRNT structure tests ───────────────────────────────────────────────────────

def test_chemostatted_n_species():
    """Only 2 dynamic species in the reduced network."""
    assert BRUSSELATOR_CHEM.n_species == 2


def test_chemostatted_crnt_summary_mentions_chemostatted():
    """CRNT summary lists chemostatted species."""
    summary = BRUSSELATOR_CHEM.crnt_summary()
    assert "Chemostatted" in summary or "chemostatted" in summary


# ── Backward-compatibility tests ──────────────────────────────────────────────

def test_no_chemostatted_is_unchanged():
    """Network without chemostatted behaves exactly as before."""
    rn = CRNetwork.from_string(
        BRUSSELATOR_REACTIONS,
        rates=BRUSSELATOR_RATES,
    )
    # All 6 species present
    assert set(rn.species) == {"A", "B", "D", "E", "X", "Y"}
    # Empty chemostatted property
    assert rn.chemostatted == {}


def test_no_chemostatted_full_odes():
    """Without chemostatted, all 6 species get ODEs."""
    rn = CRNetwork.from_string(
        BRUSSELATOR_REACTIONS,
        rates=BRUSSELATOR_RATES,
    )
    odes = rn.odes()
    assert set(odes.keys()) == {"A", "B", "D", "E", "X", "Y"}
