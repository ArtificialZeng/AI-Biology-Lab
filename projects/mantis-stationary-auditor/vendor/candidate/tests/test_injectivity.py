"""Tests for the Craciun–Feinberg injectivity test."""
from mantis import CRNetwork


def test_reversible_pair_is_injective():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 2.0})
    res = net.is_injective()
    assert res.injective
    assert res.is_sign_definite
    assert "CERTIFIED" in str(res)


def test_michaelis_menten_is_injective():
    # MM has a unique steady state per class → injective, even though it is not
    # weakly reversible (deficiency theorems alone would not settle this).
    mm = CRNetwork.from_string(
        ["E + S <-> ES", "ES -> E + P"],
        rates={"E + S -> ES": 1.0, "ES -> E + S": 1.0, "ES -> E + P": 1.0},
    )
    assert mm.is_injective().injective


def test_schlogl_is_not_certified_injective():
    # The Schlögl model is genuinely multistationary → the test must NOT certify it.
    sch = CRNetwork.from_string(
        ["A + 2 X <-> 3 X", "X <-> B"],
        rates={"A + 2X -> 3X": 6.0, "3X -> A + 2X": 1.0, "X -> B": 11.0, "B -> X": 6.0},
        chemostatted={"A": 1.0, "B": 1.0},
    )
    res = sch.is_injective()
    assert not res.injective
    assert not res.is_sign_definite
    assert "inconclusive" in str(res).lower()


def test_injectivity_is_rate_independent():
    # The verdict comes from a symbolic, rate-free determinant: same answer for any rates.
    net1 = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 1.0})
    net2 = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1e6, "B -> A": 1e-3})
    assert net1.is_injective().injective == net2.is_injective().injective
