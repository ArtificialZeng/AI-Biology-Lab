"""Independent exact pathway and state-preservation review tests."""
from fractions import Fraction
import warnings

import pytest
from cobra import Metabolite, Model, Reaction
from cobra.flux_analysis import fastcc


def signed_path(reverse, minimum):
    model = Model('independent_signed_path')
    model.solver = 'glpk'
    a, b, c = [Metabolite(x) for x in ('a', 'b', 'c')]
    specs = [('u', {a: 1}), ('p', {a: -1, b: 1}),
             ('d', {b: -1}), ('dead', {a: -1, c: 1}),
             ('closed', {a: -1, b: 1})]
    reactions = []
    for name, stoich in specs:
        reaction = Reaction(name)
        reaction.add_metabolites(stoich)
        reaction.bounds = (minimum, 10) if name in {'u', 'p', 'd'} else (0, 10)
        reactions.append(reaction)
    model.add_reactions(reactions)
    model.reactions.closed.bounds = (0, 0)
    if reverse:
        for reaction in reactions[:3]:
            reaction *= -1
    model.objective = model.reactions.d
    # The C balance proves dead=0; the closed bound proves closed=0.
    # A and B then force u=p=d (after undoing the orientation).
    # This rational witness proves that all remaining reactions can be nonzero.
    sign = -1 if reverse else 1
    witness = [Fraction(7 * sign)] * 3 + [Fraction(0)] * 2
    for metabolite in model.metabolites:
        assert sum(Fraction(r.metabolites.get(metabolite, 0)) * v
                   for r, v in zip(model.reactions, witness)) == 0
    for reaction, value in zip(model.reactions, witness):
        assert reaction.lower_bound <= value <= reaction.upper_bound
    return model


def state(model):
    return {
        'reactions': [(r.id, sorted((m.id, s) for m, s in r.metabolites.items()), r.bounds)
                      for r in model.reactions],
        'variables': sorted((v.name, v.lb, v.ub, v.type) for v in model.variables),
        'constraints': sorted((c.name, str(c.expression), c.lb, c.ub) for c in model.constraints),
        'objective': (str(model.objective.expression), model.objective.direction),
        'contexts': len(model._contexts),
    }


@pytest.mark.parametrize(('reverse', 'minimum'), [(False, 0), (False, 1), (True, 1)])
def test_exact_signed_path_and_context(reverse, minimum):
    model = signed_path(reverse, minimum)
    # Exercise nested model contexts and preservation of preexisting solver state.
    with model:
        probe = model.problem.Variable('review_probe', lb=0, ub=2)
        gate = model.problem.Constraint(probe, name='review_gate', ub=1)
        model.add_cons_vars([probe, gate])
        before = state(model)
        with warnings.catch_warnings():
            warnings.filterwarnings('error', message='Solver status is')
            consistent = fastcc(model)
        assert {r.id for r in consistent.reactions} == {'u', 'p', 'd'}
        assert state(model) == before
        assert [(r.id, r.bounds) for r in consistent.reactions] == [
            (r.id, r.bounds) for r in model.reactions if r.id in {'u', 'p', 'd'}]


def test_reversible_and_fixed_zero_control():
    model = Model('independent_reversible_cycle')
    model.solver = 'glpk'
    a, b = Metabolite('a'), Metabolite('b')
    forward, back, closed = [Reaction(x) for x in ('f', 'b', 'closed')]
    forward.add_metabolites({a: -1, b: 1})
    back.add_metabolites({a: 1, b: -1})
    closed.add_metabolites({a: -1, b: 1})
    forward.bounds = back.bounds = (-10, 10)
    closed.bounds = (0, 0)
    model.add_reactions([forward, back, closed])
    before = state(model)
    assert {r.id for r in fastcc(model).reactions} == {'f', 'b'}
    assert state(model) == before


def test_exact_closed_positive_path():
    model = signed_path(False, 0)
    model.reactions.p.bounds = (0, 0)
    # The fixed p=0, closed=0, and C balance force all five fluxes to zero.
    before = state(model)
    assert not fastcc(model).reactions
    assert state(model) == before
