"""Unit tests for the stochastic kernel layer (mantis._kernels)."""
import numpy as np

from mantis import CRNetwork
from mantis._kernels import build_kernel_arrays, normalize_seed, _propensities
from mantis.parsing import parse_reactions
from mantis.analysis import AVOGADRO


def test_normalize_seed_deterministic_and_masked():
    assert normalize_seed(42) == 42
    assert 0 <= normalize_seed(2**40 + 7) <= 0xFFFFFFFF
    # None yields *some* valid 32-bit seed.
    s = normalize_seed(None)
    assert 0 <= s <= 0xFFFFFFFF


def test_build_kernel_arrays_change_and_ceff():
    # 2A -> B with k = 4.0, volume so that N_A*V = 10.
    rxns = parse_reactions(["2 A -> B"])
    species = ["A", "B"]
    na_v = 10.0
    c_eff, react_sp, react_co, change = build_kernel_arrays(
        rxns, species, {"2A -> B": 4.0}, set(), na_v
    )
    # change vector: A:-2, B:+1
    assert list(change[0]) == [-2, 1]
    # order 2 → c = k/(N_A V)^(1) = 0.4; divided by 2! = 0.2
    assert c_eff[0] == np.float64(4.0 / na_v / 2.0)
    assert react_sp[0, 0] == 0 and react_co[0, 0] == 2


def test_propensity_matches_hand_calculation():
    # A + B -> C, k chosen so c_eff = 1.0 for a clean check.
    rxns = parse_reactions(["A + B -> C"])
    species = ["A", "B", "C"]
    na_v = 1.0
    c_eff, react_sp, react_co, _ = build_kernel_arrays(
        rxns, species, {"A + B -> C": 1.0}, set(), na_v
    )
    n = np.array([5, 3, 0], dtype=np.int64)
    out = np.empty(1, dtype=np.float64)
    _propensities(n, c_eff, react_sp, react_co, out)
    # bimolecular A·B propensity = c_eff * 5 * 3
    assert out[0] == 15.0


def test_chemostatted_species_folded_out_of_kernel():
    # A chemostatted; reaction A + X -> Y becomes pseudo-first-order in X.
    net = CRNetwork.from_string(
        ["A + X -> Y"], rates={"A + X -> Y": 2.0}, chemostatted={"A": 3.0},
    )
    r = net.stochastic_simulate(
        {"X": 100}, (0.0, 10.0), volume_L=1.0 / AVOGADRO,
        initial_as="count", seed=0,
    )
    # X is consumed, Y produced; A never appears as a dynamic species.
    assert "A" not in net.species
    assert r.counts["Y"][-1] > 0
    assert r.counts["X"][-1] < 100


def test_ssa_reproducible_and_conserves_mass():
    net = CRNetwork.from_string(
        ["A <-> B"], rates={"A -> B": 1.0, "B -> A": 1.0},
    )
    kw = dict(t_span=(0.0, 5.0), volume_L=1.0 / AVOGADRO, initial_as="count", seed=7)
    r1 = net.stochastic_simulate({"A": 200}, **kw)
    r2 = net.stochastic_simulate({"A": 200}, **kw)
    assert np.array_equal(r1.counts["A"], r2.counts["A"])
    total = r1.counts["A"] + r1.counts["B"]
    assert np.all(total == 200)
