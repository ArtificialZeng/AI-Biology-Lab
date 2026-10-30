"""Reaction graph visualization, phase portrait, and bifurcation diagrams."""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import matplotlib.axes
    import networkx as nx


_SUBSCRIPT_DIGITS = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def _prettify_species(name: str) -> str:
    """Turn ASCII species names into typographically nicer forms.

    Trailing digits become LaTeX subscripts (`H1` → `H$_{1}$`); underscores that
    join species into a complex name become the centred dot (`H1H2_CP` →
    `H$_{1}$H$_{2}$$\\cdot$CP`).
    """
    parts = name.split("_")
    pretty_parts = []
    for part in parts:
        # subscript any run of digits that follows an alphabetic prefix
        out, buf, in_digits = [], [], False
        for ch in part:
            if ch.isdigit():
                buf.append(ch)
                in_digits = True
            else:
                if in_digits:
                    out.append(f"$_{{{''.join(buf)}}}$")
                    buf = []
                    in_digits = False
                out.append(ch)
        if in_digits:
            out.append(f"$_{{{''.join(buf)}}}$")
        pretty_parts.append("".join(out))
    return "$\\cdot$".join(pretty_parts)


def _complex_label(c) -> str:
    """Human-readable label for a Complex (frozenset of (name, coeff) pairs)."""
    parts = []
    for name, coeff in sorted(c, key=lambda x: x[0]):
        sp = _prettify_species(name)
        parts.append(f"{coeff}{sp}" if coeff > 1 else sp)
    return " + ".join(parts) if parts else "∅"


def _layout_components(
    G: nx.DiGraph,
    layout_fn,
    seed: int = 42,
    grid_spacing: float = 3.5,
) -> dict:
    """
    Lay out each weakly connected component independently, then arrange
    the components in a grid.  Prevents disconnected LCs from collapsing
    onto each other (spring layout has no inter-component repulsion).
    Components are sorted largest-first so they occupy the top-left cells.
    """
    import numpy as np
    import networkx as nx

    wccs = sorted(nx.weakly_connected_components(G), key=len, reverse=True)
    n = len(wccs)
    cols = max(1, int(np.ceil(np.sqrt(n))))

    pos: dict = {}
    for i, wcc in enumerate(wccs):
        row, col = divmod(i, cols)
        sub = G.subgraph(wcc)

        if len(wcc) == 1:
            sub_pos = {next(iter(wcc)): np.array([0.0, 0.0])}
        else:
            try:
                sub_pos = layout_fn(sub, seed=seed)
            except TypeError:
                sub_pos = layout_fn(sub)
            except Exception:
                sub_pos = nx.spring_layout(sub, seed=seed)

        # Normalise to [-1, 1] so each component occupies the same area
        coords = np.array(list(sub_pos.values()), dtype=float)
        if len(coords) > 1:
            coords -= coords.mean(axis=0)
            scale = np.max(np.abs(coords))
            if scale > 1e-9:
                coords /= scale

        center = np.array([col * grid_spacing, -row * grid_spacing])
        for node, coord in zip(sub_pos.keys(), coords):
            pos[node] = coord + center

    return pos


# A small qualitative palette for linkage classes (colour-blind safe, muted).
_LC_PALETTE = [
    "#4C78A8", "#F58518", "#54A24B", "#E45756",
    "#72B7B2", "#B279A2", "#FF9DA6", "#9D755D",
    "#BAB0AC", "#EECA3B",
]


def draw_reaction_graph(
    G: nx.DiGraph,
    ax=None,
    layout: str = "kamada_kawai",
    font_size: int = 9,
    annotate_stats: bool = True,
    show_legend: bool = True,
) -> matplotlib.axes.Axes:
    """
    Publication-style rendering of the *complex reaction graph*.

    Each node is a chemical complex (multiset of species) — *not* a single
    species. This is the convention required by Feinberg's framework: the
    deficiency δ = n − ℓ − s counts complexes (n) and linkage classes (ℓ)
    of this graph. A species-resolved view (Petri net / SR-graph) would be
    a different object; see ``draw_species_reaction_graph`` for that.

    Reversible reactions are rendered as a single double-headed arrow.
    Each weakly connected component (linkage class) is tinted with its
    own colour so the partition that enters the deficiency formula is
    visually explicit.
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import networkx as nx
    import numpy as np

    if ax is None:
        _, ax = plt.subplots(figsize=(11, 7))

    layouts = {
        "spring":       nx.spring_layout,
        "shell":        nx.shell_layout,
        "spectral":     nx.spectral_layout,
        "planar":       nx.planar_layout,
        "kamada_kawai": nx.kamada_kawai_layout,
    }
    layout_fn = layouts.get(layout, nx.kamada_kawai_layout)
    pos = _layout_components(G, layout_fn)

    # Set data limits from node positions (set_axis_off disables autoscale).
    if pos:
        xs = [p[0] for p in pos.values()]
        ys = [p[1] for p in pos.values()]
        pad_x = max(0.6, 0.18 * (max(xs) - min(xs) if len(xs) > 1 else 1.0))
        pad_y = max(0.6, 0.18 * (max(ys) - min(ys) if len(ys) > 1 else 1.0))
        ax.set_xlim(min(xs) - pad_x, max(xs) + pad_x)
        ax.set_ylim(min(ys) - pad_y, max(ys) + pad_y)

    # Assign each linkage class a colour.
    wccs = sorted(nx.weakly_connected_components(G), key=len, reverse=True)
    lc_color = {}
    for i, wcc in enumerate(wccs):
        c = _LC_PALETTE[i % len(_LC_PALETTE)]
        for node in wcc:
            lc_color[node] = c

    # ── 1. Labels via rounded-box patches (drawn first so edges can clip to them) ──
    label_artists = {}
    for node, (x, y) in pos.items():
        txt = _complex_label(node)
        artist = ax.text(
            x, y, txt,
            ha="center", va="center",
            fontsize=font_size,
            family="STIXGeneral",
            color="#1a1a2e",
            zorder=3,
            bbox=dict(
                boxstyle="round,pad=0.45",
                fc="white",
                ec=lc_color[node],
                lw=1.4,
            ),
        )
        label_artists[node] = artist

    # Force a draw so each text patch has a real bounding box for arrow clipping.
    ax.figure.canvas.draw()

    # ── 2. Edges — collapse reversible pairs into a single double-headed arrow ──
    seen = set()
    for u, v in G.edges():
        if (u, v) in seen or (v, u) in seen:
            continue
        reversible = G.has_edge(v, u)
        seen.add((u, v))

        arrowstyle = "<|-|>" if reversible else "-|>"
        ax.annotate(
            "",
            xy=pos[v], xycoords="data",
            xytext=pos[u], textcoords="data",
            arrowprops=dict(
                arrowstyle=arrowstyle,
                color="#2C3E50",
                lw=1.3,
                shrinkA=22, shrinkB=22,           # leave room for the label boxes
                connectionstyle="arc3,rad=0.0",
                mutation_scale=14,
            ),
            zorder=1,
        )

    # ── 3. Stats annotation (n, ℓ, s, δ if computable) ──
    if annotate_stats:
        n = G.number_of_nodes()
        ell = len(wccs)
        stats_lines = [f"n = {n}  complexes", f"ℓ = {ell}  linkage class{'es' if ell != 1 else ''}"]
        # The graph alone doesn't carry s or δ, so we just print n and ℓ.
        ax.text(
            0.01, 0.98,
            "\n".join(stats_lines),
            transform=ax.transAxes,
            ha="left", va="top",
            fontsize=9, family="STIXGeneral",
            bbox=dict(boxstyle="round,pad=0.4", fc="#FAFAFA",
                      ec="#BDC3C7", lw=0.8),
            zorder=4,
        )

    # ── 4. Legend mapping colour → linkage class ──
    if show_legend and len(wccs) > 1:
        handles = [
            mpatches.Patch(
                facecolor="white",
                edgecolor=_LC_PALETTE[i % len(_LC_PALETTE)],
                linewidth=1.4,
                label=f"Linkage class {i+1}",
            )
            for i in range(len(wccs))
        ]
        ax.legend(
            handles=handles,
            loc="upper right",
            frameon=True,
            framealpha=0.95,
            fontsize=8,
            handlelength=1.4,
            prop={"family": "STIXGeneral"},
        )

    ax.set_axis_off()
    ax.margins(0.18)
    return ax


def plot_bifurcation(
    result: BifurcationResult,
    species: str,
    ax=None,
    stable_color: str = "#2196F3",
    unstable_color: str = "#F44336",
) -> matplotlib.axes.Axes:
    """
    Bifurcation diagram: x = parameter value (log scale), y = steady-state [species].
    Solid lines = stable branches, dashed lines = unstable branches.

    Branches are identified by sorting steady states by concentration at each parameter
    value and connecting same-index entries across consecutive parameter points.
    Stability changes within a branch (e.g. Hopf bifurcation) are handled by splitting
    the branch into segments of consistent stability.
    """
    import matplotlib.pyplot as plt
    from .analysis import BifurcationResult

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))

    # Determine the maximum number of simultaneous branches
    max_branches = max((len(ssl) for ssl in result.steady_states), default=0)
    if max_branches == 0:
        return ax

    # Collect (pval, conc, is_stable) per branch, sorted by concentration for
    # consistent branch ordering across parameter values.
    branch_data: list[list[tuple[float, float, bool]]] = [[] for _ in range(max_branches)]
    for pval, ss_list in zip(result.parameter_values, result.steady_states):
        sorted_ss = sorted(ss_list, key=lambda s: s.concentrations.get(species, 0.0))
        for i, ss in enumerate(sorted_ss):
            branch_data[i].append((
                float(pval),
                ss.concentrations.get(species, float("nan")),
                ss.is_stable,
            ))

    # Plot each branch as one or more line segments, splitting on stability changes.
    legend_added: set[str] = set()

    def _plot_segment(xs, ys, stable):
        color = stable_color if stable else unstable_color
        ls = "-" if stable else "--"
        label_str = "Stable" if stable else "Unstable"
        lbl = label_str if label_str not in legend_added else "_nolegend_"
        legend_added.add(label_str)
        ax.plot(xs, ys, ls, color=color, linewidth=2.0, label=lbl, zorder=3)

    for branch in branch_data:
        if not branch:
            continue
        cur_stable = branch[0][2]
        seg_x: list[float] = [branch[0][0]]
        seg_y: list[float] = [branch[0][1]]

        for pval, conc, is_stable in branch[1:]:
            if is_stable != cur_stable:
                # Close the current segment and start a new one at the transition
                seg_x.append(pval)
                seg_y.append(conc)
                _plot_segment(seg_x, seg_y, cur_stable)
                seg_x, seg_y = [pval], [conc]
                cur_stable = is_stable
            else:
                seg_x.append(pval)
                seg_y.append(conc)

        _plot_segment(seg_x, seg_y, cur_stable)

    ax.set_xscale("log")
    ax.set_xlabel(result.parameter_name)
    ax.set_ylabel(f"[{species}] (M)")
    ax.set_title(f"Bifurcation: {result.parameter_name}")
    ax.legend(framealpha=0.7)
    return ax


# Forward reference fix
try:
    from .analysis import BifurcationResult
except ImportError:
    pass
