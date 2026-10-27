import numpy as np
import pytest
import msprime

class TestTerminalCoalescenceRate:
    @pytest.mark.parametrize("validation", [False, True])
    @pytest.mark.parametrize("steps", [[100], [0, 100], [0, 50, 100]])
    @pytest.mark.parametrize("sizes", [(1000, 10000), (10000, 1000), (0.5, 0.25)])
    def test_size_change(self, validation, steps, sizes):
        old_size, new_size = sizes
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=old_size)
        demography.add_population_parameters_change(
            time=100, population="A", initial_size=new_size
        )
        debugger = demography.debug()
        rates, survival = debugger.coalescence_rate_trajectory(
            steps, {"A": 2}, double_step_validation=validation
        )
        extended_rates, extended_survival = debugger.coalescence_rate_trajectory(
            steps + [101], {"A": 2}, double_step_validation=validation
        )
        np.testing.assert_allclose(rates, extended_rates[:-1])
        np.testing.assert_allclose(survival, extended_survival[:-1])
        assert rates[-1] == pytest.approx(
            1 / (2 * max(1, new_size)), rel=1e-12, abs=0
        )
        assert survival[-1] == pytest.approx(np.exp(-100 / (2 * max(1, old_size))))

    @pytest.mark.parametrize("validation", [False, True])
    @pytest.mark.parametrize("proportion", [0.25, 0.5, 1])
    def test_lineage_movement(self, validation, proportion):
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=1000)
        demography.add_population(name="B", initial_size=10000)
        demography.add_mass_migration(
            time=100, source="A", dest="B", proportion=proportion
        )
        debugger = demography.debug()
        rates, survival = debugger.coalescence_rate_trajectory(
            [0, 100], {"A": 2}, double_step_validation=validation
        )
        extended_rates, extended_survival = debugger.coalescence_rate_trajectory(
            [0, 100, 101], {"A": 2}, double_step_validation=validation
        )
        expected_rate = (1 - proportion) ** 2 / 2000 + proportion**2 / 20000
        assert rates[-1] == pytest.approx(expected_rate)
        assert survival[-1] == pytest.approx(np.exp(-0.05))
        np.testing.assert_allclose(rates, extended_rates[:-1])
        np.testing.assert_allclose(survival, extended_survival[:-1])

    @pytest.mark.parametrize("validation", [False, True])
    def test_multiple_lineage_movements(self, validation):
        demography = msprime.Demography()
        for name, size in [("A", 1000), ("B", 10000), ("C", 100000)]:
            demography.add_population(name=name, initial_size=size)
        demography.add_mass_migration(time=100, source="A", dest="B", proportion=0.5)
        demography.add_mass_migration(time=100, source="B", dest="C", proportion=0.5)
        rates, survival = demography.debug().coalescence_rate_trajectory(
            [0, 100], {"A": 2}, double_step_validation=validation
        )
        # At the event time each surviving lineage has probabilities
        # (1/2, 1/4, 1/4) of being in A, B, C, respectively.
        expected_rate = 0.5**2 / 2000 + 0.25**2 / 20000 + 0.25**2 / 200000
        assert rates[-1] == pytest.approx(expected_rate)
        assert survival[-1] == pytest.approx(np.exp(-0.05))

    @pytest.mark.parametrize("steps", [[0], [100], [0, 100]])
    def test_constant_size_singleton_and_zero(self, steps):
        debugger = msprime.Demography.isolated_model([1000]).debug()
        rates, survival = debugger.coalescence_rate_trajectory(steps, {0: 2})
        np.testing.assert_allclose(rates, 1 / 2000)
        np.testing.assert_allclose(survival, np.exp(-np.array(steps) / 2000))

    def test_zero_survival(self):
        debugger = msprime.Demography.isolated_model([1]).debug()
        rates, survival = debugger.coalescence_rate_trajectory([0, 2000], {0: 2})
        assert survival[-1] == 0
        assert np.isnan(rates[-1])


