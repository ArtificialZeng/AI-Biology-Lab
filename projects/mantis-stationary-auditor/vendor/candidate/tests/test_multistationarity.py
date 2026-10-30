"""Tests for the Conradi–Feliu–Mincheva–Wiuf multistationarity-region analysis."""
from mantis import CRNetwork


def test_schlogl_region_exists():
    # The Schlögl model is genuinely multistationary → a region must exist.
    sch = CRNetwork.from_string(
        ["A + 2 X <-> 3 X", "X <-> B"],
        rates={"A + 2X -> 3X": 6.0, "3X -> A + 2X": 1.0, "X -> B": 11.0, "B -> X": 6.0},
        chemostatted={"A": 1.0, "B": 1.0},
    )
    res = sch.multistationarity_region()
    assert not res.monostationary
    assert res.multistationary_possible
    assert res.region_conditions          # at least one coefficient can be negative
    assert "REGION EXISTS" in str(res)


def test_reversible_pair_is_monostationary():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 2.0})
    res = net.multistationarity_region()
    assert res.monostationary
    assert not res.multistationary_possible
    assert res.region_conditions == []
    assert "MONOSTATIONARY" in str(res)


def test_michaelis_menten_is_monostationary():
    # MM is injective (unique steady state per class) → monostationary, and the
    # determinant test settles it even though it is not weakly reversible.
    mm = CRNetwork.from_string(
        ["E + S <-> ES", "ES -> E + P"],
        rates={"E + S -> ES": 1.0, "ES -> E + S": 1.0, "ES -> E + P": 1.0},
    )
    assert mm.multistationarity_region().monostationary


def test_agrees_with_injectivity():
    # The critical-function verdict must be consistent with the injectivity test:
    # an injective network is necessarily monostationary.
    for reactions, rates, chem in [
        (["A <-> B"], {"A -> B": 1.0, "B -> A": 2.0}, None),
        (["E + S <-> ES", "ES -> E + P"],
         {"E + S -> ES": 1.0, "ES -> E + S": 1.0, "ES -> E + P": 1.0}, None),
    ]:
        net = CRNetwork.from_string(reactions, rates=rates, chemostatted=chem)
        if net.is_injective().injective:
            assert net.multistationarity_region().monostationary


def _schlogl(k):
    return CRNetwork.from_string(
        ["A + 2 X <-> 3 X", "X <-> B"],
        rates={"A + 2X -> 3X": 6.0, "3X -> A + 2X": 1.0, "X -> B": k, "B -> X": 6.0},
        chemostatted={"A": 1.0, "B": 1.0},
    )


def test_numeric_rates_classified_against_region():
    res = _schlogl(11.0).multistationarity_region()
    # k = 11 sits inside the bistable window ⇒ exact verdict is True.
    assert res.multistationary_at_rates is True


def test_exact_discriminant_region_boundary():
    import sympy
    res = _schlogl(11.0).multistationarity_region()
    # Steady states reduce to a univariate cubic ⇒ an exact boundary is available.
    assert res.region_boundary is not None
    assert res.steady_state_polynomial is not None
    # Substituting the Schlögl constants must give the cubic's discriminant in k:
    #   Δ(k) = −4k³ + 36k² + 648k − 6156   (zero at the two folds k = 10.72, 11.15).
    k = sympy.Symbol("k_3")
    subs = {sympy.Symbol("k_1"): 6, sympy.Symbol("k_2"): 1, sympy.Symbol("k_4"): 6,
            sympy.Symbol("A"): 1, sympy.Symbol("B"): 1}
    disc = sympy.expand(res.region_boundary.subs(subs))
    expected = sympy.expand(-4 * k**3 + 36 * k**2 + 648 * k - 6156)
    assert sympy.simplify(disc - expected) == 0


def test_exact_at_rates_tracks_the_bistable_window():
    # The exact positive-root count flips across the analytic fold window
    # (10.72, 11.15) — the structural sign condition alone could not resolve this.
    assert _schlogl(10.0).multistationarity_region().multistationary_at_rates is False
    assert _schlogl(11.0).multistationarity_region().multistationary_at_rates is True
    assert _schlogl(12.0).multistationarity_region().multistationary_at_rates is False


def test_is_multistationary_exact_by_enumeration():
    # Fix 2: definitive count of positive steady states in the given class.
    assert _schlogl(11.0).is_multistationary({"X": 0.5}) is True
    assert _schlogl(10.0).is_multistationary({"X": 0.5}) is False
    assert _schlogl(12.0).is_multistationary({"X": 0.5}) is False
    # A monostationary network is never multistationary, for any totals.
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 2.0})
    assert net.is_multistationary({"A": 3.0}) is False
