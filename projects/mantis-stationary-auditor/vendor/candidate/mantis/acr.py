"""Absolute Concentration Robustness (ACR) detection — Shinar & Feinberg, *Science* 2010.

A mass-action system has **absolute concentration robustness** in a species *S* when
the steady-state concentration of *S* is identical in *every* positive steady state,
regardless of the initial totals (the stoichiometric compatibility class).  This is
exactly the robustness property a biosensor or homeostatic module is designed to
have, which makes it directly relevant to diagnostic-circuit design.

Shinar & Feinberg give a purely structural *sufficient* condition:

    If a mass-action system has deficiency one, admits a positive steady state, and
    possesses two **non-terminal** complexes that differ only in species *S*, then the
    system has ACR in *S*.

"Non-terminal" complexes are those that do not lie in a *terminal strong linkage
class* (a strongly connected component of the reaction graph with no reactions
leading out of it).  Two complexes "differ only in *S*" when their species vectors
differ in exactly one coordinate.

This module implements that test on top of the complex graph mantis already builds.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import networkx as nx
import numpy as np

from .parsing import Complex, Reaction


@dataclass
class ACRResult:
    """Outcome of an ACR analysis.

    Attributes
    ----------
    has_acr : bool
        Whether the Shinar–Feinberg sufficient condition detected ACR.
    species : list[str]
        Species for which ACR is certified (empty when ``has_acr`` is False).
    deficiency : int
        Network deficiency (the test applies only when this is 1).
    applies : bool
        Whether the deficiency-one precondition is met.
    acr_values : dict[str, float]
        Steady-state value of each ACR species, when a positive steady state was
        computed (the value is independent of the initial totals).
    nonterminal_complex_pairs : list[tuple[str, str, str]]
        Witnessing ``(complex_a, complex_b, species)`` triples.
    summary_lines : list[str]
    """
    has_acr: bool
    species: list[str] = field(default_factory=list)
    deficiency: int = 0
    applies: bool = False
    acr_values: dict[str, float] = field(default_factory=dict)
    nonterminal_complex_pairs: list[tuple[str, str, str]] = field(default_factory=list)
    summary_lines: list[str] = field(default_factory=list)

    def __str__(self) -> str:
        return "\n".join(self.summary_lines)


def _complex_label(c: Complex) -> str:
    if not c:
        return "0"
    return " + ".join(
        f"{coeff}{name}" if coeff > 1 else name
        for name, coeff in sorted(c, key=lambda x: x[0])
    )


def _complex_vector(c: Complex, species: list[str]) -> np.ndarray:
    idx = {s: i for i, s in enumerate(species)}
    v = np.zeros(len(species))
    for name, coeff in c:
        if name in idx:
            v[idx[name]] = coeff
    return v


def terminal_strong_linkage_classes(graph: nx.DiGraph) -> list[set[Complex]]:
    """Return the terminal strong linkage classes of a complex graph.

    A strong linkage class (SCC) is *terminal* when no edge leaves it for another
    SCC — i.e. it is a sink in the condensation.
    """
    sccs = list(nx.strongly_connected_components(graph))
    condensation = nx.condensation(graph, scc=sccs)
    terminal = []
    for node in condensation.nodes():
        if condensation.out_degree(node) == 0:
            terminal.append(set(condensation.nodes[node]["members"]))
    return terminal


def detect_acr(
    reactions: list[Reaction],
    species: list[str],
    graph: nx.DiGraph,
    deficiency: int,
    rate_values: dict[str, float] | None = None,
    initial_conditions: dict[str, float] | None = None,
    chemostatted_values: dict[str, float] | None = None,
) -> ACRResult:
    """Detect ACR via the Shinar–Feinberg deficiency-one sufficient condition.

    ``graph`` is the (reduced) complex digraph; ``deficiency`` the network deficiency.
    When ``rate_values`` and ``initial_conditions`` are supplied and a positive steady
    state exists, the robust steady-state value of each ACR species is also reported.
    """
    if deficiency != 1:
        return ACRResult(
            has_acr=False,
            deficiency=deficiency,
            applies=False,
            summary_lines=[
                "ACR (Shinar–Feinberg deficiency-one test): NOT applicable "
                f"(δ = {deficiency} ≠ 1).",
                "  → This sufficient condition only covers deficiency-one networks; "
                "absence here does not rule out ACR by other mechanisms.",
            ],
        )

    terminal = terminal_strong_linkage_classes(graph)
    terminal_nodes: set[Complex] = set().union(*terminal) if terminal else set()
    nonterminal = [c for c in graph.nodes() if c not in terminal_nodes]

    acr_species: list[str] = []
    pairs: list[tuple[str, str, str]] = []
    for i in range(len(nonterminal)):
        for j in range(i + 1, len(nonterminal)):
            ca, cb = nonterminal[i], nonterminal[j]
            diff = _complex_vector(ca, species) - _complex_vector(cb, species)
            nz = np.flatnonzero(np.abs(diff) > 1e-12)
            if len(nz) == 1:
                s = species[int(nz[0])]
                if s not in acr_species:
                    acr_species.append(s)
                pairs.append((_complex_label(ca), _complex_label(cb), s))

    has_acr = len(acr_species) > 0

    acr_values: dict[str, float] = {}
    if has_acr and rate_values is not None and initial_conditions is not None:
        from .analysis import find_steady_states
        states = find_steady_states(
            reactions, species, rate_values, initial_conditions,
            n_attempts=12, chemostatted_values=chemostatted_values or None,
        )
        for st in states:
            y = np.array([st.concentrations[s] for s in species])
            if np.all(y >= 0) and st.residual < 1e-4 * (np.max(np.abs(y)) + 1e-30):
                for s in acr_species:
                    acr_values.setdefault(s, float(st.concentrations[s]))
                break

    if has_acr:
        lines = [
            "ACR (Shinar–Feinberg deficiency-one test): DETECTED",
            f"  → Robust species: {', '.join(sorted(acr_species))}",
            "  → Their steady-state concentrations are independent of the initial "
            "totals (within any class admitting a positive steady state).",
        ]
        for a, b, s in pairs:
            lines.append(f"     witness: non-terminal complexes ‘{a}’ and ‘{b}’ differ only in {s}")
        for s, v in sorted(acr_values.items()):
            lines.append(f"     {s}* = {v:.6g} (computed at a positive steady state)")
    else:
        lines = [
            "ACR (Shinar–Feinberg deficiency-one test): NOT detected",
            "  → δ = 1 but no two non-terminal complexes differ in a single species.",
        ]

    return ACRResult(
        has_acr=has_acr,
        species=sorted(acr_species),
        deficiency=deficiency,
        applies=True,
        acr_values=acr_values,
        nonterminal_complex_pairs=pairs,
        summary_lines=lines,
    )
