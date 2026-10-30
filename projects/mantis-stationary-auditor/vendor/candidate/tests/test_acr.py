"""Tests for Absolute Concentration Robustness detection (Shinar–Feinberg 2010)."""
import pytest

from mantis import CRNetwork


# Canonical minimal ACR network: X + Y -> 2Y, Y -> X.
# Deficiency 1; non-terminal complexes {X+Y} and {Y} differ only in X ⇒ ACR in X,
# with robust value X* = k2/k1.
ACR_NET = ["X + Y -> 2 Y", "Y -> X"]


def test_minimal_acr_detected():
    net = CRNetwork.from_string(ACR_NET, rates={"X + Y -> 2Y": 3.0, "Y -> X": 6.0})
    assert net.deficiency == 1
    res = net.detect_acr({"X": 5.0, "Y": 5.0})
    assert res.has_acr
    assert res.species == ["X"]


def test_acr_value_is_k_ratio():
    net = CRNetwork.from_string(ACR_NET, rates={"X + Y -> 2Y": 3.0, "Y -> X": 6.0})
    res = net.detect_acr({"X": 5.0, "Y": 5.0})
    assert res.acr_values["X"] == pytest.approx(6.0 / 3.0, abs=1e-6)


def test_acr_value_independent_of_totals():
    """The robust species value must not change when the initial totals change."""
    net = CRNetwork.from_string(ACR_NET, rates={"X + Y -> 2Y": 3.0, "Y -> X": 6.0})
    v1 = net.detect_acr({"X": 5.0, "Y": 5.0}).acr_values["X"]
    v2 = net.detect_acr({"X": 50.0, "Y": 200.0}).acr_values["X"]
    assert v1 == pytest.approx(v2, rel=1e-6)


def test_deficiency_zero_not_applicable():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 1.0})
    res = net.detect_acr()
    assert not res.applies
    assert not res.has_acr


def test_cha_deficiency_one_but_no_acr():
    cha = [
        "miR21 + H1 <-> miR21_H1",
        "miR21_H1 + H2 <-> H1H2 + miR21",
        "H1H2 + CP <-> H1H2_CP",
        "H1 + H2 <-> H1H2",
    ]
    net = CRNetwork.from_string(cha)
    res = net.detect_acr()
    assert res.applies          # δ = 1
    assert not res.has_acr      # no single-species-difference among non-terminal complexes
