"""Independent public-contract checks; contains no empty biochemical reactions."""
import numpy as np
import pytest

from cobra import Configuration, Metabolite, Model, Reaction
from cobra.flux_analysis.find_cyclic_reactions import find_cyclic_reactions
from cobra.flux_analysis.loopless import add_loopless
from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.util import create_stoichiometric_matrix

Configuration().processes = 1


def reaction(name, metabolites, bounds=(0, 10)):
    result = Reaction(name)
    result.add_metabolites(metabolites)
    result.bounds = bounds
    return result


def diamond(cycle=False):
    """DAG with a mixed-sign null vector, unequal arm caps, custom constraint."""
    m = Model("independent_diamond")
    m.solver = "glpk"
    a, b, c, d = [Metabolite(x) for x in "ABCD"]
    m.add_reactions([
        reaction("supply", {a: 1}),
        reaction("AB", {a: -1, b: 1}, (0, 3)),
        reaction("BD", {b: -1, d: 1}, (0, 3)),
        reaction("AC", {a: -1, c: 1}, (0, 7)),
        reaction("CD", {c: -1, d: 1}, (0, 7)),
        reaction("drain", {d: -1}),
    ])
    m.objective = m.reactions.drain
    m.add_cons_vars(m.problem.Constraint(m.reactions.AB.flux_expression,
                                       ub=2, name="external_capacity"))
    if cycle:
        e, f, g = [Metabolite(x) for x in "EFG"]
        m.add_reactions([
            reaction("EF", {e: -1, f: 1}, (1, 10)),
            reaction("FG", {f: -1, g: 1}, (1, 10)),
            reaction("GE", {g: -1, e: 1}, (1, 10)),
        ])
    return m


def boundary_model():
    m = Model("independent_boundary")
    m.solver = "glpk"
    a = Metabolite("feed")
    m.add_reactions([reaction("supply", {a: 1}, (0, 9)),
                     reaction("drain", {a: -1}, (0, 9))])
    m.objective = m.reactions.drain
    return m


def fingerprint(m):
    """Check structure, bounds, objective, and user constraints, not solve status."""
    return (
        tuple((v.name, v.lb, v.ub, v.type) for v in m.variables),
        tuple((c.name, str(c.expression), c.lb, c.ub) for c in m.constraints),
        str(m.objective.expression), m.objective.direction,
        tuple((r.id, r.bounds) for r in m.reactions),
    )


@pytest.mark.parametrize("method", [None, "fastSNP", "original"])
def test_diamond_default_and_original(method):
    m = diamond()
    before = fingerprint(m)
    assert m.slim_optimize() == pytest.approx(9)
    kwargs = {} if method is None else {"method": method}
    assert add_loopless(m, **kwargs) is None
    assert m.slim_optimize() == pytest.approx(9)
    if method != "original":
        assert fingerprint(m) == before
    else:
        # The original matrix nullspace has a mixed-sign cycle; retaining the
        # pre-existing formulation must not cut the feasible DAG optimum.
        assert str(m.objective.expression) == before[2]
        assert "external_capacity" in m.constraints


@pytest.mark.parametrize("method", [None, "fastSNP", "original"])
def test_empty_model_no_auxiliary_state(method):
    m = Model("empty")
    m.solver = "glpk"
    before = fingerprint(m)
    kwargs = {} if method is None else {"method": method}
    assert add_loopless(m, **kwargs) is None
    assert fingerprint(m) == before


@pytest.mark.parametrize("method", [None, "fastSNP", "original"])
def test_default_full_network_cyclic_control(method):
    m = diamond(cycle=True)
    assert m.slim_optimize() == pytest.approx(9)
    kwargs = {} if method is None else {"method": method}
    add_loopless(m, **kwargs)
    m.slim_optimize()
    assert m.solver.status == "infeasible"


@pytest.mark.parametrize("method", ["fastSNP", "original"])
@pytest.mark.parametrize("selection", [[], ["supply", "drain"], ["AB"],
                                         ["EF", "FG"], ["EF", "FG", "GE"]])
def test_subset_contract_and_cyclic_control(method, selection):
    m = diamond(cycle=True)
    before = fingerprint(m)
    assert m.slim_optimize() == pytest.approx(9)
    assert add_loopless(m, method=method, reactions=selection) is None
    m.slim_optimize()
    if selection == ["EF", "FG", "GE"]:
        assert m.solver.status == "infeasible"
    else:
        assert m.solver.status == "optimal"
        assert m.slim_optimize() == pytest.approx(9)
        assert fingerprint(m) == before


@pytest.mark.parametrize("method", ["fastSNP", "original"])
@pytest.mark.parametrize("selection", [[], ["EF", "FG", "GE"]])
def test_selection_generator_consumed_once(method, selection):
    m = diamond(cycle=True)
    consumed = []
    def ids():
        for name in selection:
            consumed.append(name)
            yield name
    gen = ids()
    add_loopless(m, method=method, reactions=gen)
    assert consumed == selection
    assert list(gen) == []
    m.slim_optimize()
    assert m.solver.status == ("infeasible" if selection else "optimal")


@pytest.mark.parametrize("method", ["fastSNP", "original"])
@pytest.mark.parametrize("empty", [True, False])
def test_context_restoration(method, empty):
    m = diamond(cycle=True)
    before = fingerprint(m)
    with m:
        add_loopless(m, method=method, reactions=[] if empty else ["EF", "FG", "GE"])
        m.slim_optimize()
        assert m.solver.status == ("optimal" if empty else "infeasible")
    assert fingerprint(m) == before
    assert m.slim_optimize() == pytest.approx(9)


@pytest.mark.parametrize("method", ["basic", "optimized"])
@pytest.mark.parametrize("kind", ["empty", "boundary", "dag", "cycle"])
def test_finder_closes_boundaries_preserves_model(method, kind):
    m = Model() if kind == "empty" else (boundary_model() if kind == "boundary"
        else diamond(cycle=kind == "cycle"))
    m.solver = "glpk"
    before = fingerprint(m)
    np.random.seed(832)
    found, directions = find_cyclic_reactions(m, method=method)
    assert fingerprint(m) == before
    expected = {"EF", "FG", "GE"} if kind == "cycle" else set()
    assert set(found) == expected
    assert len(directions) == len(found)
    assert all(d == (False, True) for d in directions)
    # Both models have through-flow but it is not an internal circulation.
    assert not {"supply", "drain"}.intersection(found)


@pytest.mark.parametrize("method", ["fastSNP", "original"])
def test_boundary_only_add_loopless(method):
    m = boundary_model()
    before = fingerprint(m)
    add_loopless(m, method=method)
    assert fingerprint(m) == before
    assert m.slim_optimize() == pytest.approx(9)


@pytest.mark.parametrize("kind", ["dag", "boundary"])
@pytest.mark.parametrize("algorithm", ["fastSNP", "cycleFreeFlux"])
@pytest.mark.parametrize("all_fluxes", [False, True])
def test_fva_analytic_bounds_vectors_and_no_side_effects(kind, algorithm, all_fluxes):
    m = diamond() if kind == "dag" else boundary_model()
    before = fingerprint(m)
    value = flux_variability_analysis(m, loopless=algorithm, all_fluxes=all_fluxes,
                                     fraction_of_optimum=0.5, processes=1)
    bounds = value[0] if all_fluxes else value
    expected = {"supply": (4.5, 9), "drain": (4.5, 9)}
    if kind == "dag":
        expected.update(AB=(0, 2), BD=(0, 2), AC=(2.5, 7), CD=(2.5, 7))
    for rid, limits in expected.items():
        assert bounds.loc[rid, "minimum"] == pytest.approx(limits[0])
        assert bounds.loc[rid, "maximum"] == pytest.approx(limits[1])
    if all_fluxes:
        s = create_stoichiometric_matrix(m)
        ids = [r.id for r in m.reactions]
        for direction, frame in value[1].items():
            for target, row in frame.iterrows():
                vector = row.loc[ids].to_numpy()
                assert np.max(np.abs(s @ vector)) < 1e-7
                assert row["drain"] >= 4.5 - 1e-7
                assert row[target] == pytest.approx(bounds.loc[target, direction])
                for rid in ids:
                    lb, ub = m.reactions.get_by_id(rid).bounds
                    assert lb - 1e-7 <= row[rid] <= ub + 1e-7
                if kind == "dag":
                    assert row["AB"] <= 2 + 1e-7
    assert fingerprint(m) == before
    assert m.slim_optimize() == pytest.approx(9)


@pytest.mark.parametrize("routine", ["add", "finder"])
def test_validation_precedes_empty_return(routine):
    m = Model()
    call = add_loopless if routine == "add" else find_cyclic_reactions
    with pytest.raises(ValueError):
        call(m, method="invalid")
    with pytest.raises(ValueError):
        call(m, zero_cutoff=m.tolerance / 2)
    if routine == "finder":
        with pytest.raises(ValueError):
            call(m, required_stop_checks_num=0)
