"""SBML import/export round-trip tests (optional python-libsbml dependency)."""
import pytest

pytest.importorskip("libsbml", reason="python-libsbml not installed")

from mantis import CRNetwork


CHA_REACTIONS = [
    "miR21 + H1 <-> miR21_H1",
    "miR21_H1 + H2 <-> H1H2 + miR21",
    "H1H2 + CP <-> H1H2_CP",
    "H1 + H2 <-> H1H2",
]


def _structure(net):
    return (net.n_species, net.n_complexes, net.n_linkage_classes,
            net.deficiency, net.is_weakly_reversible)


def test_export_returns_sbml_string():
    net = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 2.0, "B -> A": 0.5})
    xml = net.to_sbml()
    assert isinstance(xml, str)
    assert "<sbml" in xml
    assert "listOfReactions" in xml


def test_roundtrip_preserves_structure():
    net = CRNetwork.from_string(CHA_REACTIONS)
    xml = net.to_sbml()
    net2 = CRNetwork.from_sbml(xml)
    assert _structure(net) == _structure(net2)


def test_roundtrip_preserves_rates():
    rates = {
        "miR21 + H1 -> miR21_H1": 3.0e6,
        "miR21_H1 -> miR21 + H1": 0.1,
        "H1 + H2 -> H1H2": 1.2e3,
    }
    net = CRNetwork.from_string(CHA_REACTIONS, rates=rates)
    net2 = CRNetwork.from_sbml(net.to_sbml())
    # Rate keys are stored in canonical (alphabetically-sorted) form; compare the
    # full canonical dicts rather than the raw user keys.
    assert net._rates.keys() == net2._rates.keys()
    for key in net._rates:
        assert net2._rates[key] == pytest.approx(net._rates[key], rel=1e-9)


def test_roundtrip_stoichiometric_coefficients():
    net = CRNetwork.from_string(["2 X + Y -> 3 X"], rates={"2X + Y -> 3X": 1.5})
    net2 = CRNetwork.from_sbml(net.to_sbml())
    # The stoichiometry matrix is invariant to the round-trip.
    import numpy as np
    assert np.array_equal(net.stoichiometry_matrix, net2.stoichiometry_matrix)
    assert net2._rates["2X + Y -> 3X"] == pytest.approx(1.5)


def test_roundtrip_chemostatted_species():
    net = CRNetwork.from_string(
        ["A -> X", "X -> A"],
        rates={"A -> X": 1.0, "X -> A": 2.0},
        chemostatted={"A": 5.0},
    )
    xml = net.to_sbml()
    assert "boundaryCondition=\"true\"" in xml
    net2 = CRNetwork.from_sbml(xml)
    assert net2.chemostatted == pytest.approx({"A": 5.0})
    assert "A" not in net2.species


def test_file_write_and_read(tmp_path):
    net = CRNetwork.from_string(CHA_REACTIONS, rates={"H1 + H2 -> H1H2": 7.0})
    path = tmp_path / "model.xml"
    net.to_sbml(str(path))
    assert path.exists()
    net2 = CRNetwork.from_sbml(str(path))
    assert _structure(net) == _structure(net2)
    assert net2._rates["H1 + H2 -> H1H2"] == pytest.approx(7.0)


def test_missing_dependency_message(monkeypatch):
    """If libsbml is absent, the error names the optional extra."""
    import builtins
    import mantis.sbml as sbml_mod

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "libsbml":
            raise ImportError("no libsbml")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="mantis-delta\\[sbml\\]"):
        sbml_mod._require_libsbml()
