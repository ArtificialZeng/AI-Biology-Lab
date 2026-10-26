"""Small custom cases collected through the real public MEMOTE API."""

import pytest
from memote.support.consistency_helpers import get_internals


@pytest.mark.biomass
def test_current_biomass(model, reaction_id):
    assert reaction_id == "branch"


def test_reuse_first(model):
    assert get_internals(model) is get_internals(model)


def test_reuse_second(model):
    assert get_internals(model) is get_internals(model)
