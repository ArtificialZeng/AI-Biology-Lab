"""Tests for Finite State Projection (Munsky–Khammash 2006)."""
import numpy as np
from scipy.stats import binom
import pytest

from mantis import CRNetwork

NA = 6.02214076e23
UNIT_VOL = 1.0 / NA


def test_fsp_stationary_matches_binomial():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 2.0, "B -> A": 3.0})
    N = 30
    fsp = net.fsp({"A": N, "B": 0}, volume_L=UNIT_VOL, initial_as="count")
    assert fsp.states.shape == (N + 1, 2)
    assert fsp.truncation_error == pytest.approx(0.0, abs=1e-9)
    vals, probs = fsp.marginal("A")
    assert np.max(np.abs(probs - binom.pmf(vals, N, 0.6))) < 1e-12


def test_fsp_stationary_matches_ack_product_form():
    """FSP (linear-algebra) and the ACK closed form must agree exactly."""
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 4.0})
    N = 25
    fsp = net.fsp({"A": N, "B": 0}, volume_L=UNIT_VOL, initial_as="count")
    ack = net.stationary_distribution({"A": N, "B": 0}, volume_L=UNIT_VOL, initial_as="count")
    _, fsp_marg = fsp.marginal("A")
    _, ack_marg = ack.marginal("A")
    assert np.max(np.abs(fsp_marg - ack_marg)) < 1e-10


def test_fsp_transient_conserves_probability_for_closed_network():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 2.0, "B -> A": 3.0})
    N = 30
    res = net.fsp({"A": N, "B": 0}, volume_L=UNIT_VOL, t=0.1, initial_as="count")
    assert res.probabilities.sum() == pytest.approx(1.0, abs=1e-9)
    assert res.truncation_error < 1e-9
    assert res.time == 0.1


def test_fsp_transient_mean_matches_ssa():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 2.0, "B -> A": 3.0})
    N = 30
    res = net.fsp({"A": N, "B": 0}, volume_L=UNIT_VOL, t=0.1, initial_as="count")
    finals = [
        net.stochastic_simulate({"A": N, "B": 0}, (0.0, 0.1), volume_L=UNIT_VOL,
                                initial_as="count", seed=s).counts["A"][-1]
        for s in range(400)
    ]
    assert res.expected_counts()["A"] == pytest.approx(np.mean(finals), rel=0.05)


def test_fsp_probability_lookup():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 2.0, "B -> A": 3.0})
    N = 20
    fsp = net.fsp({"A": N, "B": 0}, volume_L=UNIT_VOL, initial_as="count")
    # Sum over the line A + B = N must be 1.
    total = sum(fsp.probability({"A": a, "B": N - a}) for a in range(N + 1))
    assert total == pytest.approx(1.0, abs=1e-9)
    # A state off the compatibility class has zero probability.
    assert fsp.probability({"A": N, "B": N}) == 0.0


def test_fsp_three_species_chain_stationary_normalised():
    net = CRNetwork.from_string(
        ["A <-> B", "B <-> C"],
        rates={"A -> B": 1.0, "B -> A": 2.0, "B -> C": 1.5, "C -> B": 1.0},
    )
    fsp = net.fsp({"A": 12, "B": 0, "C": 0}, volume_L=UNIT_VOL, initial_as="count")
    assert fsp.probabilities.sum() == pytest.approx(1.0, abs=1e-9)
    assert fsp.truncation_error == pytest.approx(0.0, abs=1e-9)
