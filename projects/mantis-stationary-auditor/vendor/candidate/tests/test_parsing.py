import pytest
from mantis.parsing import (
    parse_complex, parse_reaction_string, parse_reactions,
    canonical_rate_key, normalize_rate_key, Reaction
)


def test_parse_complex_basic():
    c = parse_complex("A + 2B")
    assert c == frozenset({("A", 1), ("B", 2)})


def test_parse_complex_no_coeff():
    c = parse_complex("A")
    assert c == frozenset({("A", 1)})


def test_parse_complex_reordered():
    assert parse_complex("2B + A") == parse_complex("A + 2B")


def test_parse_complex_underscore_species():
    c = parse_complex("miR21_H1 + H2")
    assert c == frozenset({("miR21_H1", 1), ("H2", 1)})


def test_parse_complex_invalid_name():
    # Pure digits with no species name following is invalid
    with pytest.raises(ValueError, match="Invalid species name"):
        parse_complex("123 + A")


def test_parse_reaction_irreversible():
    rxns = parse_reaction_string("A + 2B -> C")
    assert len(rxns) == 1
    r = rxns[0]
    assert r.reactants == frozenset({("A", 1), ("B", 2)})
    assert r.products == frozenset({("C", 1)})


def test_parse_reaction_reversible():
    rxns = parse_reaction_string("A + 2B <-> C")
    assert len(rxns) == 2
    keys = {r.rate_key for r in rxns}
    assert "A + 2B -> C" in keys
    assert "C -> A + 2B" in keys


def test_rate_key_normalization():
    # "2B + A -> C" and "A + 2B -> C" should produce the same canonical key
    r1 = parse_reaction_string("2B + A -> C")[0]
    r2 = parse_reaction_string("A + 2B -> C")[0]
    assert r1.rate_key == r2.rate_key
    assert r1.rate_key == "A + 2B -> C"


def test_normalize_rate_key():
    assert normalize_rate_key("H2 + miR21_H1 -> H1H2 + miR21") == \
           normalize_rate_key("miR21_H1 + H2 -> miR21 + H1H2")


def test_parse_reactions_list():
    rxns = parse_reactions(["A <-> B", "B -> C"])
    assert len(rxns) == 3


def test_missing_arrow():
    with pytest.raises(ValueError):
        parse_reaction_string("A + B")
