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


def test_conditioned_probability_is_zero_outside_support():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 1.0})
    sd = net.stationary_distribution({"A": 2, "B": 0}, volume_L=UNIT_VOL, initial_as="count")
    assert sd.probability({"A": 0, "B": 0}) == 0.0
    assert sd.probability({"A": 3, "B": 0}) == 0.0


@pytest.mark.parametrize("initial_a, expected_states", [(2, {(2, 0), (0, 2)}), (1, {(1, 1)})])
def test_stationary_support_preserves_discrete_communicating_class(initial_a, expected_states):
    net = CRNetwork.from_string(["2A <-> 2B"], rates={"2A -> 2B": 1.0, "2B -> 2A": 1.0})
    sd = net.stationary_distribution(
        {"A": initial_a, "B": 2 - initial_a}, volume_L=UNIT_VOL, initial_as="count",
    )
    assert sd.is_exact and sd.is_conditioned
    assert {tuple(row) for row in sd._states} == expected_states
    assert sum(sd.probability({"A": a, "B": b}) for a, b in expected_states) == pytest.approx(1.0)
    for a in range(3):
        if (a, 2 - a) not in expected_states:
            assert sd.probability({"A": a, "B": 2 - a}) == 0.0
    values, probabilities = sd.marginal("A")
    mean = values @ probabilities
    variance = (values - mean) ** 2 @ probabilities
    assert variance == pytest.approx(0.0 if initial_a == 1 else 1.0)
    draws = sd.sample(size=40, seed=17)
    assert set(zip(draws["A"], draws["B"])) <= expected_states


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


def test_catalyst_below_required_count_blocks_reachability():
    net = CRNetwork.from_string(
        ["A + 2C <-> B + 2C"],
        rates={"A + 2C -> B + 2C": 1.0, "B + 2C -> A + 2C": 1.0},
    )
    sd = net.stationary_distribution(
        {"A": 2, "B": 0, "C": 1}, volume_L=UNIT_VOL, initial_as="count",
    )
    assert {tuple(row) for row in sd._states} == {(2, 0, 1)}
    assert sd.probability({"A": 2, "B": 0, "C": 1}) == pytest.approx(1.0)
    assert sd.probability({"A": 0, "B": 2, "C": 1}) == 0.0


def test_chemostatted_catalyst_uses_fixed_concentration_not_count_threshold():
    net = CRNetwork.from_string(
        ["A + 2C <-> B + 2C"],
        rates={"A + 2C -> B + 2C": 1.0, "B + 2C -> A + 2C": 1.0},
        chemostatted={"C": 1.0},
    )
    sd = net.stationary_distribution(
        {"A": 2, "B": 0}, volume_L=UNIT_VOL, initial_as="count",
    )
    assert {tuple(row) for row in sd._states} == {(0, 2), (1, 1), (2, 0)}
    assert sd.probability({"A": 1, "B": 1}) == pytest.approx(0.5)
