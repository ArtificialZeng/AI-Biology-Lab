"""Tests for Horn–Jackson global-stability certification."""
import numpy as np
import sympy
import pytest

from mantis import CRNetwork


def test_reversible_pair_is_complex_balanced():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 2.0})
    # δ=0 and weakly reversible ⇒ complex-balanced for all rates (no ICs needed).
    assert net.is_complex_balanced() is True


def test_reversible_pair_global_stability_certificate():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 2.0})
    cert = net.certify_global_stability({"A": 3.0}, n_samples=100)
    assert cert.is_complex_balanced
    assert cert.globally_stable
    assert cert.max_dVdt <= 1e-9
    # Equilibrium k1 A = k2 B, A + B = 3 → A = 2, B = 1.
    assert cert.equilibrium["A"] == pytest.approx(2.0, abs=1e-6)
    assert cert.equilibrium["B"] == pytest.approx(1.0, abs=1e-6)
    assert cert.lyapunov_function is not None
    assert cert.lyapunov_derivative is not None


def test_lyapunov_derivative_is_nonpositive_and_zero_at_equilibrium():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 2.0})
    cert = net.certify_global_stability({"A": 3.0}, n_samples=10)
    A, B = sympy.Symbol("A"), sympy.Symbol("B")
    dV = cert.lyapunov_derivative
    # Zero at the equilibrium c*.
    at_star = float(dV.subs({A: cert.equilibrium["A"], B: cert.equilibrium["B"]}))
    assert at_star == pytest.approx(0.0, abs=1e-9)
    # Strictly negative away from c* but inside the same class (A+B=3).
    away = float(dV.subs({A: 2.5, B: 0.5}))
    assert away < 0.0


def test_triangle_cycle_global_stability():
    tri = CRNetwork.from_string(
        ["A <-> B", "B <-> C", "C <-> A"],
        rates={"A -> B": 1.0, "B -> A": 1.0, "B -> C": 2.0,
               "C -> B": 1.0, "C -> A": 1.0, "A -> C": 1.0},
    )
    assert tri.deficiency == 0 and tri.is_weakly_reversible
    cert = tri.certify_global_stability({"A": 1.0, "B": 1.0, "C": 1.0}, n_samples=100)
    assert cert.globally_stable
    assert cert.max_dVdt <= 1e-9


def test_non_weakly_reversible_not_certified():
    mm = CRNetwork.from_string(
        ["E + S <-> ES", "ES -> E + P"],
        rates={"E + S -> ES": 1.0, "ES -> E + S": 1.0, "ES -> E + P": 1.0},
    )
    assert mm.is_complex_balanced({"E": 1.0, "S": 1.0}) is False
    cert = mm.certify_global_stability({"E": 1.0, "S": 1.0})
    assert cert.is_complex_balanced is False
    assert cert.globally_stable is False
    assert "NOT certified" in str(cert)


def test_complex_balance_requires_weak_reversibility():
    # An irreversible open reaction is not weakly reversible → cannot be complex-balanced.
    net = CRNetwork.from_string(["A -> B"], rates={"A -> B": 1.0})
    assert net.is_complex_balanced() is False


def test_birch_point_is_the_positive_equilibrium():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 2.0})
    birch = net.complex_balanced_equilibrium({"A": 3.0})
    assert birch["A"] == pytest.approx(2.0, abs=1e-6)
    assert birch["B"] == pytest.approx(1.0, abs=1e-6)


def test_birch_point_requires_complex_balanced():
    mm = CRNetwork.from_string(
        ["E + S <-> ES", "ES -> E + P"],
        rates={"E + S -> ES": 1.0, "ES -> E + S": 1.0, "ES -> E + P": 1.0},
    )
    with pytest.raises(ValueError, match="not complex-balanced"):
        mm.complex_balanced_equilibrium({"E": 1.0, "S": 1.0})
