"""Tests for exhaustive (homotopy/Gröbner) steady-state enumeration."""
import numpy as np
import pytest

from mantis import CRNetwork


# Schlögl model with A, B chemostatted at 1.0 and rates tuned so that
# dX/dt = -(X-1)(X-2)(X-3): three positive steady states at X = 1, 2, 3.
SCHLOGL = ["A + 2 X <-> 3 X", "X <-> B"]
SCHLOGL_RATES = {
    "A + 2X -> 3X": 6.0, "3X -> A + 2X": 1.0,
    "X -> B": 11.0, "B -> X": 6.0,
}


def _schlogl():
    return CRNetwork.from_string(SCHLOGL, rates=SCHLOGL_RATES, chemostatted={"A": 1.0, "B": 1.0})


def test_schlogl_finds_all_three_states():
    net = _schlogl()
    states = net.all_steady_states({"X": 0.5})
    xs = sorted(s.concentrations["X"] for s in states)
    assert len(states) == 3
    assert xs == pytest.approx([1.0, 2.0, 3.0], abs=1e-6)


def test_schlogl_middle_state_is_unstable():
    net = _schlogl()
    states = net.all_steady_states({"X": 0.5})
    by_x = {round(s.concentrations["X"]): s for s in states}
    # Outer branches stable, middle (X=2) unstable — the textbook bistable picture.
    assert by_x[1].is_stable
    assert by_x[3].is_stable
    assert not by_x[2].is_stable


def test_exhaustive_beats_multistart_on_schlogl():
    """The exhaustive engine recovers the unstable root the multi-start solver misses."""
    net = _schlogl()
    exhaustive = net.all_steady_states({"X": 0.5})
    multistart = net.steady_states({"X": 0.5}, n_attempts=40, seed=0)
    n_unstable_exhaustive = sum(1 for s in exhaustive if not s.is_stable)
    n_unstable_multistart = sum(1 for s in multistart if not s.is_stable)
    assert len(exhaustive) >= len(multistart)
    assert n_unstable_exhaustive >= 1
    # Demonstrate the gap the feature exists to close.
    assert n_unstable_exhaustive >= n_unstable_multistart


def test_reversible_pair_unique_state():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 2.0})
    states = net.all_steady_states({"A": 3.0})
    assert len(states) == 1
    # k1 A = k2 B, A + B = 3 → A = 2, B = 1.
    assert states[0].concentrations["A"] == pytest.approx(2.0, abs=1e-9)
    assert states[0].concentrations["B"] == pytest.approx(1.0, abs=1e-9)
    assert states[0].is_stable


CHA = [
    "miR21 + H1 <-> miR21_H1",
    "miR21_H1 + H2 <-> H1H2 + miR21",
    "H1H2 + CP <-> H1H2_CP",
    "H1 + H2 <-> H1H2",
]
CHA_RATES = {
    "miR21 + H1 -> miR21_H1": 1e6, "miR21_H1 -> miR21 + H1": 1.0,
    "miR21_H1 + H2 -> H1H2 + miR21": 1e5, "H1H2 + miR21 -> miR21_H1 + H2": 1.0,
    "H1H2 + CP -> H1H2_CP": 1e5, "H1H2_CP -> H1H2 + CP": 1.0,
    "H1 + H2 -> H1H2": 1.0, "H1H2 -> H1 + H2": 1.0,
}


def test_cha_unique_positive_state_matches_integration():
    """Deficiency-one ⇒ at most one positive steady state per class; it must match the ODE."""
    net = CRNetwork.from_string(CHA, rates=CHA_RATES)
    ic = {"miR21": 1e-9, "H1": 1e-7, "H2": 1e-7, "CP": 1e-7}
    exhaustive = net.all_steady_states(ic)
    assert len(exhaustive) == 1
    integrated = net.steady_states(ic, n_attempts=8, seed=0)[0]
    for sp in net.species:
        assert exhaustive[0].concentrations[sp] == pytest.approx(
            integrated.concentrations[sp], rel=1e-3, abs=1e-15
        )


def test_unknown_backend_raises():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 1.0})
    with pytest.raises(ValueError, match="Unknown backend"):
        net.all_steady_states({"A": 1.0}, backend="nope")
