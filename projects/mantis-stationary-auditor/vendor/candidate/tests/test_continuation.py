"""Tests for pseudo-arclength continuation (fold / Hopf detection)."""
import numpy as np
import pytest

from mantis import CRNetwork


# Schlögl with A, B chemostatted at 1.0.  Steady states solve
#   X³ − 6X² + k·X − 6 = 0   (k = rate of X → B),
# whose folds (also 3X² − 12X + k = 0) sit at the analytic points
#   (X, k) = (1.34727…, 10.7218…) and (2.53209…, 11.1505…).
SCHLOGL = ["A + 2 X <-> 3 X", "X <-> B"]
SCHLOGL_RATES = {
    "A + 2X -> 3X": 6.0, "3X -> A + 2X": 1.0,
    "X -> B": 11.0, "B -> X": 6.0,
}


def _schlogl():
    return CRNetwork.from_string(SCHLOGL, rates=SCHLOGL_RATES,
                                 chemostatted={"A": 1.0, "B": 1.0})


def test_continuation_traces_the_s_curve():
    net = _schlogl()
    res = net.continuation("X -> B", (1.0, 20.0), {"X": 0.5})
    # The branch must fold back: λ is non-monotonic in arclength.
    lam = res.parameter_values
    assert np.any(np.diff(lam) < 0)
    # It spans the lower (X≈0.3) and upper (X≈6) branches.
    xs = res.species_branch("X")
    assert xs.min() < 0.6
    assert xs.max() > 5.0


def test_continuation_finds_both_folds_accurately():
    net = _schlogl()
    res = net.continuation("X -> B", (1.0, 20.0), {"X": 0.5})
    folds = res.folds()
    assert len(folds) == 2
    by_x = sorted(folds, key=lambda f: f.state["X"])
    # Analytic fold locations.
    assert by_x[0].state["X"] == pytest.approx(1.34727, abs=1e-2)
    assert by_x[0].parameter == pytest.approx(10.7218, abs=2e-2)
    assert by_x[1].state["X"] == pytest.approx(2.53209, abs=1e-2)
    assert by_x[1].parameter == pytest.approx(11.1505, abs=2e-2)
    # At a fold a real eigenvalue passes through zero.
    for f in folds:
        assert abs(f.eigenvalue.real) < 1e-4
        assert abs(f.eigenvalue.imag) < 1e-6


def test_fold_window_brackets_the_bistable_region():
    net = _schlogl()
    res = net.continuation("X -> B", (1.0, 20.0), {"X": 0.5})
    ks = sorted(f.parameter for f in res.folds())
    # The nominal rate k = 11 (three steady states) lies inside the fold window.
    assert ks[0] < 11.0 < ks[1]


def test_monostable_branch_has_no_folds():
    # A simple reversible conversion has a single steady state for all parameters.
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 2.0})
    res = net.continuation("A -> B", (0.1, 10.0), {"A": 3.0})
    assert res.folds() == []
    assert res.bifurcations == []
    # Every traced point is stable.
    assert np.all(res.stable)


def test_hopf_detected_on_the_brusselator():
    # Brusselator with A, B chemostatted: a Hopf occurs when the effective
    # parameter k·B crosses 1 + A² = 2 (A = 1).  Continuing the rate of
    # "B + X -> Y + D" (B = 3) the threshold is k = 2/3.
    bru = CRNetwork.from_string(
        ["A -> X", "2 X + Y -> 3 X", "B + X -> Y + D", "X -> E"],
        rates={"A -> X": 1.0, "2X + Y -> 3X": 1.0,
               "B + X -> Y + D": 1.0, "X -> E": 1.0},
        chemostatted={"A": 1.0, "B": 3.0, "D": 0.0, "E": 0.0},
    )
    res = bru.continuation("B + X -> Y + D", (0.5, 4.0), {"X": 1.0, "Y": 3.0})
    hopfs = res.hopfs()
    assert len(hopfs) >= 1
    h = hopfs[0]
    assert h.parameter == pytest.approx(2.0 / 3.0, abs=2e-2)
    # A Hopf crossing: a complex pair on (near) the imaginary axis.
    assert abs(h.eigenvalue.real) < 1e-2
    assert abs(h.eigenvalue.imag) > 0.1


def test_explicit_initial_state_selects_a_branch():
    net = _schlogl()
    # Seed on the upper branch at k = 11 (X = 3) and continue.
    res = net.continuation("X -> B", (1.0, 20.0), {"X": 0.5},
                           initial_state={"X": 3.0})
    assert len(res.branch) > 3
    assert res.species_branch("X").max() > 2.5
