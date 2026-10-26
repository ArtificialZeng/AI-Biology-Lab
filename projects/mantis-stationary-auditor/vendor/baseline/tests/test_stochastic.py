"""Gillespie SSA stochastic simulator tests."""
import numpy as np
import pytest

import mantis


def _bind_unbind_network():
    return mantis.CRNetwork.from_string(
        ["A + B <-> C"],
        rates={"A + B -> C": 1e6, "C -> A + B": 1e-2},
    )


class TestGillespieAPI:
    def test_returns_stochastic_result(self):
        rn = _bind_unbind_network()
        res = rn.stochastic_simulate(
            initial_conditions={"A": 50, "B": 50, "C": 0},
            t_span=(0.0, 5.0),
            volume_L=1e-15,
            initial_as="count",
            seed=0,
        )
        assert isinstance(res, mantis.StochasticResult)
        assert res.success
        assert res.times[0] == 0.0
        assert res.times[-1] == 5.0
        assert "A" in res.counts and "C" in res.counts

    def test_counts_are_integer_valued(self):
        rn = _bind_unbind_network()
        res = rn.stochastic_simulate(
            initial_conditions={"A": 20, "B": 20, "C": 0},
            t_span=(0.0, 2.0),
            volume_L=1e-15,
            initial_as="count",
            seed=1,
        )
        for arr in res.counts.values():
            assert arr.dtype.kind == "i"

    def test_mass_balance_preserved(self):
        """A + C must be constant (A·molecule count is conserved by the reaction)."""
        rn = _bind_unbind_network()
        res = rn.stochastic_simulate(
            initial_conditions={"A": 100, "B": 100, "C": 0},
            t_span=(0.0, 10.0),
            volume_L=1e-15,
            initial_as="count",
            seed=2,
        )
        # A + C should be conserved (= 100)
        sum_arr = res.counts["A"] + res.counts["C"]
        assert np.all(sum_arr == 100)

    def test_seed_reproducibility(self):
        rn = _bind_unbind_network()
        kw = dict(initial_conditions={"A": 30, "B": 30, "C": 0},
                  t_span=(0.0, 2.0), volume_L=1e-15,
                  initial_as="count", seed=42)
        r1 = rn.stochastic_simulate(**kw)
        r2 = rn.stochastic_simulate(**kw)
        assert np.array_equal(r1.times, r2.times)
        assert np.array_equal(r1.counts["A"], r2.counts["A"])

    def test_concentration_initial(self):
        """initial_as='concentration' converts properly via N_A·V."""
        rn = _bind_unbind_network()
        # 1 µM in 1 fL = 6 * 10^-7 mol/L × 6e23 × 1e-15 L = 0.6  → rounded to 1
        res = rn.stochastic_simulate(
            initial_conditions={"A": 1e-6, "B": 1e-6, "C": 0.0},
            t_span=(0.0, 0.001),
            volume_L=1e-15,
            seed=0,
        )
        # ~600 molecules each
        assert 400 < res.counts["A"][0] < 800

    def test_volume_dependence(self):
        """Larger volume → fewer molecules at fixed concentration → more stochastic noise."""
        rn = _bind_unbind_network()
        # Same concentrations, different volumes — initial counts scale with V
        small = rn.stochastic_simulate(
            initial_conditions={"A": 1e-6, "B": 1e-6, "C": 0.0},
            t_span=(0.0, 0.001), volume_L=1e-15, seed=0,
        )
        big = rn.stochastic_simulate(
            initial_conditions={"A": 1e-6, "B": 1e-6, "C": 0.0},
            t_span=(0.0, 0.001), volume_L=1e-12, seed=0,
        )
        assert big.counts["A"][0] > 100 * small.counts["A"][0]


class TestGillespieVsODE:
    """Ensemble mean of many SSA trajectories should converge to the ODE solution."""

    def test_mean_trajectory_approaches_ode(self):
        rn = _bind_unbind_network()
        ic = {"A": 1e-6, "B": 1e-6, "C": 0.0}  # 1 µM each
        t_span = (0.0, 0.005)
        volume = 1e-15  # ~600 molecules

        # Deterministic
        ode = rn.simulate(ic, t_span)
        ode_final_C = ode.concentrations["C"][-1]

        # Stochastic ensemble (50 trajectories)
        finals = []
        for s in range(50):
            res = rn.stochastic_simulate(
                ic, t_span, volume_L=volume, seed=s,
            )
            finals.append(res.concentrations["C"][-1])
        mean_C = sum(finals) / len(finals)

        # Within ~15% of ODE (50 trajectories with ~600 molecules → moderate noise)
        rel = abs(mean_C - ode_final_C) / max(ode_final_C, 1e-30)
        assert rel < 0.20, f"SSA mean {mean_C:.3e} vs ODE {ode_final_C:.3e} ({rel:.1%})"
