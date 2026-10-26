"""τ-leap stochastic simulator tests."""
import numpy as np
import pytest

import mantis


def _bind_unbind_network():
    return mantis.CRNetwork.from_string(
        ["A + B <-> C"],
        rates={"A + B -> C": 1e6, "C -> A + B": 1e-2},
    )


class TestTauLeapAPI:
    def test_returns_stochastic_result(self):
        rn = _bind_unbind_network()
        res = rn.tau_leap_simulate(
            initial_conditions={"A": 500, "B": 500, "C": 0},
            t_span=(0.0, 1.0),
            volume_L=1e-15,
            initial_as="count",
            seed=0,
        )
        assert isinstance(res, mantis.StochasticResult)
        assert res.times[0] == 0.0
        assert res.times[-1] == pytest.approx(1.0)
        assert "A" in res.counts and "C" in res.counts

    def test_counts_are_integer_valued(self):
        rn = _bind_unbind_network()
        res = rn.tau_leap_simulate(
            initial_conditions={"A": 500, "B": 500, "C": 0},
            t_span=(0.0, 0.5),
            volume_L=1e-15,
            initial_as="count",
            seed=1,
        )
        for arr in res.counts.values():
            assert arr.dtype.kind == "i"

    def test_mass_balance_preserved(self):
        """A + C must be constant (mass balance for the reaction A + B <-> C)."""
        rn = _bind_unbind_network()
        res = rn.tau_leap_simulate(
            initial_conditions={"A": 1000, "B": 1000, "C": 0},
            t_span=(0.0, 2.0),
            volume_L=1e-15,
            initial_as="count",
            seed=2,
        )
        sum_arr = res.counts["A"] + res.counts["C"]
        assert np.all(sum_arr == 1000)

    def test_seed_reproducibility(self):
        rn = _bind_unbind_network()
        kw = dict(initial_conditions={"A": 300, "B": 300, "C": 0},
                  t_span=(0.0, 1.0), volume_L=1e-15,
                  initial_as="count", seed=42)
        r1 = rn.tau_leap_simulate(**kw)
        r2 = rn.tau_leap_simulate(**kw)
        assert np.array_equal(r1.times, r2.times)
        assert np.array_equal(r1.counts["A"], r2.counts["A"])

    def test_fixed_tau(self):
        """Fixed τ honors the user's choice."""
        rn = _bind_unbind_network()
        res = rn.tau_leap_simulate(
            initial_conditions={"A": 500, "B": 500, "C": 0},
            t_span=(0.0, 1.0),
            volume_L=1e-15,
            initial_as="count",
            tau=0.01,
            seed=0,
        )
        # With fixed tau=0.01 over 1 sec, ~100 leaps — well below any pathological case
        assert res.times[-1] == pytest.approx(1.0)

    def test_no_negative_counts(self):
        """Leap rejection must prevent species counts from going below zero."""
        rn = _bind_unbind_network()
        res = rn.tau_leap_simulate(
            initial_conditions={"A": 50, "B": 50, "C": 0},  # small pop → stiff
            t_span=(0.0, 1.0),
            volume_L=1e-15,
            initial_as="count",
            tau=0.01,  # deliberately too coarse
            seed=0,
        )
        for arr in res.counts.values():
            assert (arr >= 0).all(), "tau-leap drove a species below zero"

    def test_n_record(self):
        rn = _bind_unbind_network()
        res = rn.tau_leap_simulate(
            initial_conditions={"A": 500, "B": 500, "C": 0},
            t_span=(0.0, 1.0),
            volume_L=1e-15,
            initial_as="count",
            n_record=50,
            seed=0,
        )
        assert len(res.times) == 50


class TestTauLeapVsSSA:
    """τ-leap ensemble mean should approach SSA ensemble mean for large populations."""

    def test_mean_matches_ssa_at_large_population(self):
        rn = _bind_unbind_network()
        ic = {"A": 2000, "B": 2000, "C": 0}
        t_span = (0.0, 0.5)
        volume = 1e-15

        ssa_finals, leap_finals = [], []
        for s in range(20):
            r_ssa = rn.stochastic_simulate(
                initial_conditions=ic, t_span=t_span, volume_L=volume,
                initial_as="count", seed=s,
            )
            r_leap = rn.tau_leap_simulate(
                initial_conditions=ic, t_span=t_span, volume_L=volume,
                initial_as="count", seed=s,
            )
            ssa_finals.append(r_ssa.counts["C"][-1])
            leap_finals.append(r_leap.counts["C"][-1])

        ssa_mean = np.mean(ssa_finals)
        leap_mean = np.mean(leap_finals)
        # τ-leap should match SSA within ~10% at this population size
        rel_err = abs(leap_mean - ssa_mean) / max(ssa_mean, 1.0)
        assert rel_err < 0.15, f"τ-leap mean {leap_mean:.0f} vs SSA mean {ssa_mean:.0f}"
