"""CRNT graph construction, deficiency, weak reversibility, and Feinberg theorem checks."""
from dataclasses import dataclass, field

import numpy as np
import networkx as nx

from .parsing import Complex, Reaction
from .stoichiometry import build_species_list, build_stoichiometry_matrix, matrix_rank


@dataclass
class CRNTResult:
    n_species: int
    n_complexes: int
    n_reactions: int
    n_linkage_classes: int
    rank_N: int
    deficiency: int
    is_weakly_reversible: bool
    deficiency_zero_applies: bool
    deficiency_one_applies: bool
    summary_lines: list[str] = field(default_factory=list)


def build_complex_graph(reactions: list[Reaction], chemostatted: set[str] | None = None) -> nx.DiGraph:
    """
    Directed graph where nodes are Complex objects and edges are reactions.
    Each edge carries the reaction index as attribute 'rxn_idx'.

    Parameters
    ----------
    chemostatted : set[str], optional
        Species held at fixed concentrations; they are removed from each complex
        node before building the graph (reduced CRNT graph).
    """
    from .parsing import reduce_complex
    chem = chemostatted or set()
    G = nx.DiGraph()
    for idx, rxn in enumerate(reactions):
        r = reduce_complex(rxn.reactants, chem)
        p = reduce_complex(rxn.products, chem)
        G.add_node(r)
        G.add_node(p)
        G.add_edge(r, p, rxn_idx=idx)
    return G


def get_linkage_classes(G: nx.DiGraph) -> list[set[Complex]]:
    """Weakly connected components of G (undirected view)."""
    return [set(c) for c in nx.weakly_connected_components(G)]


def check_weak_reversibility(G: nx.DiGraph) -> bool:
    """
    Weakly reversible iff every weakly connected component is strongly connected.
    """
    for wcc in nx.weakly_connected_components(G):
        subgraph = G.subgraph(wcc)
        if not nx.is_strongly_connected(subgraph):
            return False
    return True


def compute_deficiency(n_complexes: int, n_linkage_classes: int, rank_N: int) -> int:
    """δ = n - l - s, always non-negative."""
    return n_complexes - n_linkage_classes - rank_N


def _per_lc_deficiency(
    lc: set[Complex],
    reactions: list[Reaction],
    species: list[str],
    chemostatted: set[str] | None = None,
) -> int:
    """Deficiency of one linkage class."""
    from .parsing import reduce_complex
    chem = chemostatted or set()
    lc_rxn_indices = [
        i for i, r in enumerate(reactions)
        if reduce_complex(r.reactants, chem) in lc
        and reduce_complex(r.products, chem) in lc
    ]
    n_l = len(lc)
    if not lc_rxn_indices:
        return n_l - 1
    N_full = build_stoichiometry_matrix(reactions, species)
    N_l = N_full[:, lc_rxn_indices]
    # Restrict to species that appear in this linkage class
    lc_species: set[str] = set()
    for c in lc:
        for name, _ in c:
            lc_species.add(name)
    sp_idx = [i for i, s in enumerate(species) if s in lc_species]
    N_sub = N_l[np.ix_(sp_idx, list(range(len(lc_rxn_indices))))]
    rank_l = matrix_rank(N_sub)
    return n_l - 1 - rank_l


def check_deficiency_one_theorem(
    deficiency: int,
    lcs: list[set[Complex]],
    reactions: list[Reaction],
    species: list[str],
    is_weakly_reversible: bool,
    chemostatted: set[str] | None = None,
) -> tuple[bool, str]:
    """
    Structural check for Deficiency One Theorem applicability:
      1. δ = 1
      2. Each linkage class has per-LC deficiency 0 or 1
      3. At most one linkage class has per-LC deficiency 1
    Returns (applies, explanation).
    """
    if deficiency != 1:
        return False, "Deficiency One Theorem: NOT applicable (δ ≠ 1)"
    lc_defs = [_per_lc_deficiency(lc, reactions, species, chemostatted) for lc in lcs]
    if any(d > 1 for d in lc_defs):
        return False, (
            "Deficiency One Theorem: NOT applicable "
            "(some linkage class has per-LC deficiency > 1)"
        )
    if sum(1 for d in lc_defs if d == 1) > 1:
        return False, (
            "Deficiency One Theorem: NOT applicable "
            "(more than one linkage class has per-LC deficiency 1)"
        )
    return True, (
        "Deficiency One Theorem: Applicable (δ=1, structural conditions met)\n"
        "  → For mass-action kinetics: admits at most one steady state per "
        "stoichiometry class."
    )


def crnt_analysis(
    reactions: list[Reaction],
    N: np.ndarray,
    species: list[str],
    chemostatted: set[str] | None = None,
) -> CRNTResult:
    chem = chemostatted or set()
    G = build_complex_graph(reactions, chem)
    complexes = list(G.nodes())
    lcs = get_linkage_classes(G)
    n = len(complexes)
    l = len(lcs)
    s = matrix_rank(N)
    delta = compute_deficiency(n, l, s)
    wr = check_weak_reversibility(G)

    dzt_applies = (delta == 0 and wr)
    d1t_applies, d1t_msg = check_deficiency_one_theorem(delta, lcs, reactions, species, wr, chem)

    lines = [
        f"CRNetwork: {len(species)} species, {n} complexes, {l} linkage classes",
        f"Stoichiometry matrix rank: {s}",
        f"Deficiency: δ = {n} - {l} - {s} = {delta}",
        f"Weakly reversible: {'Yes' if wr else 'No'}",
        "",
    ]
    if dzt_applies:
        lines.append(
            "Deficiency Zero Theorem: Applicable\n"
            "  → Unique asymptotically stable steady state per stoichiometry class.\n"
            "  → Bistability and oscillations are impossible."
        )
    else:
        reason = "δ ≠ 0" if delta != 0 else "not weakly reversible"
        lines.append(f"Deficiency Zero Theorem: NOT applicable ({reason})")
    lines.append(d1t_msg)

    return CRNTResult(
        n_species=len(species),
        n_complexes=n,
        n_reactions=len(reactions),
        n_linkage_classes=l,
        rank_N=s,
        deficiency=delta,
        is_weakly_reversible=wr,
        deficiency_zero_applies=dzt_applies,
        deficiency_one_applies=d1t_applies,
        summary_lines=lines,
    )
