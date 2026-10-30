"""Tests for the Anderson–Craciun–Kurtz product-form stationary distribution."""
import numpy as np
from scipy.stats import binom, poisson
import pytest

from mantis import CRNetwork

NA = 6.02214076e23
UNIT_VOL = 1.0 / NA  # so that molecule counts equal the (scaled) concentrations


def test_closed_network_is_conditioned_binomial():
    # A ⇌ B with k1·A = k2·B and A + B = N → Binomial(N, p), p = k2/(k1+k2).
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 2.0, "B -> A": 3.0})
    N = 40
    sd = net.stationary_distribution({"A": N, "B": 0}, volume_L=UNIT_VOL, initial_as="count")
    assert sd.is_conditioned and sd.is_exact
    mA, mB = sd.means
    p = mA / (mA + mB)
    assert p == pytest.approx(3.0 / 5.0)
    vals, probs = sd.marginal("A")
    assert np.max(np.abs(probs - binom.pmf(vals, N, p))) < 1e-12


def test_conditioned_probabilities_sum_to_one():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 4.0})
    N = 30
    sd = net.stationary_distribution({"A": N, "B": 0}, volume_L=UNIT_VOL, initial_as="count")
    total = sum(sd.probability({"A": a, "B": N - a}) for a in range(N + 1))
    assert total == pytest.approx(1.0, abs=1e-10)


def test_conditioned_expected_counts_match_binomial_mean():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 2.0, "B -> A": 3.0})
    N = 40
    sd = net.stationary_distribution({"A": N, "B": 0}, volume_L=UNIT_VOL, initial_as="count")
    p = sd.means[0] / sd.means.sum()
    assert sd.expected_counts()["A"] == pytest.approx(N * p, rel=1e-9)


def test_open_birth_death_is_poisson():
    # S -> A, A -> S with S chemostatted: A ~ Poisson(λ), λ = k1·[S]/k2 (counts).
    net = CRNetwork.from_string(
        ["S -> A", "A -> S"], rates={"S -> A": 5.0, "A -> S": 1.0}, chemostatted={"S": 2.0}
    )
    sd = net.stationary_distribution({"A": 0}, volume_L=UNIT_VOL, initial_as="count")
    assert not sd.is_conditioned
    assert sd.poisson_means()["A"] == pytest.approx(10.0, rel=1e-9)
    vals, probs = sd.marginal("A")
    assert np.max(np.abs(probs - poisson.pmf(vals, 10.0))) < 1e-12


def test_matches_ssa_ensemble_mean():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 2.0, "B -> A": 3.0})
    N = 40
    sd = net.stationary_distribution({"A": N, "B": 0}, volume_L=UNIT_VOL, initial_as="count")
    finals = [
        net.stochastic_simulate({"A": N, "B": 0}, (0.0, 5.0), volume_L=UNIT_VOL,
                                initial_as="count", seed=s).counts["A"][-1]
        for s in range(300)
    ]
    assert np.mean(finals) == pytest.approx(sd.expected_counts()["A"], rel=0.05)


def test_not_complex_balanced_raises():
    mm = CRNetwork.from_string(
        ["E + S <-> ES", "ES -> E + P"],
        rates={"E + S -> ES": 1.0, "ES -> E + S": 1.0, "ES -> E + P": 1.0},
    )
    with pytest.raises(ValueError, match="not complex-balanced"):
        mm.stationary_distribution({"E": 10, "S": 10}, volume_L=UNIT_VOL, initial_as="count")
