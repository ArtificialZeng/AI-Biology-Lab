"""Real API diagnostics; invoked only by the isolated temporary-copy runner."""

import io
import logging
import sys
from contextlib import redirect_stderr, redirect_stdout
from importlib import metadata
from pathlib import Path


def probe(expected_source):
    import memote
    from cobra import Metabolite, Model, Reaction
    from memote.suite import api
    from memote.support import basic
    from memote.support import consistency_helpers as ch
    from memote.support import helpers

    if not (Path(memote.__file__).resolve().parent == Path(expected_source).resolve()):
        raise RuntimeError('Import or profiling integrity check failed')
    logging.disable(logging.CRITICAL)
    custom = Path(__file__).with_name("custom_cases.py")

    def reset():
        for module in (basic, helpers, ch):
            for fn in vars(module).values():
                if callable(fn) and hasattr(fn, "cache"):
                    fn.clear()

    def build():
        reset()
        model = Model("cache_lifecycle_audit")
        model.solver = "glpk"
        a, b = Metabolite("a_c", compartment="c"), Metabolite("b_c", compartment="c")
        conversion = Reaction("conversion")
        conversion.add_metabolites({a: -1, b: 1})
        model.add_reactions([conversion])
        branch = Reaction("branch")
        branch.add_metabolites({a: -1, b: 3})
        return model, branch

    def run(model, selection="test_stoichiometric_consistency", extra=()):
        calls = [0]
        code_object = ch.get_internals.__wrapped__.__code__

        def profile(frame, event, arg):
            if event == "call" and frame.f_code is code_object:
                calls[0] += 1

        previous = sys.getprofile()
        try:
            sys.setprofile(profile)
            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                code, result = api.test_model(
                    model,
                    results=True,
                    pytest_args=[
                        "-q",
                        "--tb=short",
                        "-o",
                        "addopts=",
                        "-p",
                        "no:cacheprovider",
                        "-k",
                        selection,
                        *map(str, extra),
                    ],
                )
        finally:
            sys.setprofile(previous)
        if not (sys.getprofile() is previous):
            raise RuntimeError('Import or profiling integrity check failed')
        return {
            "returncode": int(code),
            "classification_computations": calls[0],
            "cases": {
                key: {
                    name: value[name]
                    for name in ("result", "data", "metric")
                    if name in value
                }
                for key, value in result.cases.items()
            },
        }

    model, branch = build()
    initial = run(model)
    model.add_reactions([branch])
    added = run(model)
    model, branch = build()
    model.add_reactions([branch])
    remove_before = run(model)
    model.remove_reactions([branch])
    removed = run(model)

    model, branch = build()
    model.add_reactions([branch])
    selection = "test_stoichiometric_consistency or test_current_biomass"
    before = run(model, selection, [custom])
    branch.annotation["sbo"] = "SBO:0000629"
    after = run(model, selection, [custom])
    fresh = run(model.copy(), selection, [custom])

    model, _ = build()
    reuse = run(model, "test_reuse_", [custom])

    class PlannedFailure(Exception):
        pass

    class FailingConfiguration:
        def load(self, current):
            ch.get_internals(current)
            raise PlannedFailure("controlled configuration failure")

    model, branch = build()
    raised = False
    try:
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            api.test_model(
                model, experimental=FailingConfiguration(), pytest_args=["--tb=short"]
            )
    except PlannedFailure:
        raised = True
    model.add_reactions([branch])
    retry = run(model)

    inventory = {
        module.__name__ + "." + name: fn
        for module in (basic, helpers, ch)
        for name, fn in vars(module).items()
        if callable(fn)
        and hasattr(fn, "cache")
        and getattr(fn, "__module__", None) == module.__name__
    }
    try:
        from memote.support.cache import clear_support_caches
    except ImportError:
        clear_support_caches = None
    cleared = None
    if clear_support_caches:
        for fn in inventory.values():
            fn.cache[("cache inventory marker",)] = object()
        clear_support_caches()
        cleared = all(len(fn.cache) == 0 for fn in inventory.values())

    def consistent(row):
        return row["cases"]["test_stoichiometric_consistency"]["data"]

    checks = {
        "initial_exact_oracle": consistent(initial) is True,
        "add_refresh": consistent(added) is False,
        "remove_refresh": consistent(removed) is True,
        "annotation_initial_oracle": consistent(before) is False,
        "annotation_fresh_copy_oracle": consistent(fresh) is True,
        "annotation_refresh": consistent(after) is True,
        "collection_refresh": after["cases"]["test_current_biomass"]["result"]
        == {"branch": "passed"},
        "within_run_reuse": reuse["returncode"] == 0
        and reuse["classification_computations"] == 1,
        "controlled_exception_seen": raised,
        "exception_retry_refresh": consistent(retry) is False,
        "ten_wrappers_retained": len(inventory) == 10
        and all(fn.size() == 2 for fn in inventory.values()),
    }
    if cleared is not None:
        checks["all_current_caches_cleared"] = cleared
    return {
        "imported_temporary_source": True,
        "versions": {
            name: metadata.version(name)
            for name in ("memote", "cobra", "optlang", "swiglpk", "pytest", "pylru")
        },
        "checks": checks,
        "inventory": sorted(inventory),
        "reaction_sequence": [initial, added, remove_before, removed],
        "annotation_sequence": [before, after, fresh],
        "reuse": reuse,
        "retry": retry,
    }
