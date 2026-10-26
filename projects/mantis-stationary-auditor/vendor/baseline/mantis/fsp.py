"""Finite State Projection — Munsky & Khammash (*J. Chem. Phys.* 2006).

The chemical master equation (CME) is an (often infinite) linear ODE for the
probability of every molecular-count state.  FSP **truncates** the state space to a
finite, reachable subset, assembles the corresponding generator matrix, and solves
the truncated CME directly — giving accurate distributions exactly where the SSA is
noisiest (the low-copy-number biosensor regime), together with a rigorous bound on
the truncation error (the probability that has leaked out of the projection).

This module enumerates the reachable states (BFS under the reaction stoichiometry,
respecting conservation laws, capped at ``max_states``), builds the sparse generator
``Q``, and provides:

* **transient** distributions ``p(t) = expm(Q·t)·p₀`` with the leaked-mass error bound;
* **stationary** distributions (null space of ``Q``) for closed, fully enumerated
  networks — which, for complex-balanced networks, reproduces the exact
  Anderson–Craciun–Kurtz product form (a strong cross-check).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .analysis import AVOGADRO
from .parsing import Reaction


@dataclass
class FSPResult:
    """A finite-state-projection distribution over molecular-count states.

    Attributes
    ----------
    species : list[str]
    states : np.ndarray
        ``(M, n)`` integer array of enumerated states.
    probabilities : np.ndarray
        Length-``M`` probability vector aligned with ``states``.
    truncation_error : float
        Probability that left the projection (0 for a closed, fully enumerated
        network; an upper bound on the distribution error otherwise).
    time : float | None
        Evaluation time for a transient solve, or ``None`` for the stationary one.
    truncated : bool
        True if state enumeration hit ``max_states`` before closing.
    """
    species: list[str]
    states: np.ndarray
    probabilities: np.ndarray
    truncation_error: float
    time: float | None = None
    truncated: bool = False
    _index: dict = field(default=None, repr=False)

    def probability(self, counts: dict[str, float]) -> float:
        """Probability of an exact integer state (0 if outside the projection)."""
        key = tuple(int(counts.get(s, 0)) for s in self.species)
        if self._index is None:
            self._index = {tuple(int(v) for v in row): i
                           for i, row in enumerate(self.states)}
        i = self._index.get(key)
        return float(self.probabilities[i]) if i is not None else 0.0

    def marginal(self, species: str) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(values, probabilities)`` for one species' marginal."""
        i = self.species.index(species)
        col = self.states[:, i]
        vals = np.arange(col.min(), col.max() + 1)
        out = np.zeros(len(vals))
        for v, p in zip(col, self.probabilities):
            out[v - col.min()] += p
        return vals, out

    def expected_counts(self) -> dict[str, float]:
        """Mean molecule counts under the (renormalised) projected distribution."""
        w = self.probabilities / self.probabilities.sum()
        mean = (self.states * w[:, None]).sum(axis=0)
        return {s: float(m) for s, m in zip(self.species, mean)}


def _propensity_fn(c_eff, react_sp, react_co):
    """Return a function state → propensity vector (stochastic mass action)."""
    n_rxn = c_eff.shape[0]
    maxr = react_sp.shape[1]

    def propensities(x):
        a = np.zeros(n_rxn)
        for j in range(n_rxn):
            if c_eff[j] <= 0.0:
                continue
            h = 1.0
            ok = True
            for r in range(maxr):
                idx = react_sp[j, r]
                if idx < 0:
                    break
                coeff = react_co[j, r]
                ni = x[idx]
                if ni < coeff:
                    ok = False
                    break
                ff = 1.0
                for d in range(coeff):
                    ff *= (ni - d)
                h *= ff
            a[j] = c_eff[j] * h if ok else 0.0
        return a

    return propensities


def _enumerate_reachable(n0, change, propensities, max_states):
    """BFS over reachable states; returns (states_array, index_map, truncated)."""
    start = tuple(int(v) for v in n0)
    index = {start: 0}
    order = [start]
    frontier = [np.array(start, dtype=np.int64)]
    truncated = False

    while frontier:
        x = frontier.pop()
        a = propensities(x)
        for j in range(change.shape[0]):
            if a[j] <= 0.0:
                continue
            y = x + change[j]
            if np.any(y < 0):
                continue
            key = tuple(int(v) for v in y)
            if key not in index:
                if len(index) >= max_states:
                    truncated = True
                    continue
                index[key] = len(order)
                order.append(key)
                frontier.append(y)

    states = np.array(order, dtype=np.int64)
    return states, index, truncated


def _build_generator(states, index, change, propensities):
    """Assemble the sparse CME generator ``Q`` and total per-state leak rate."""
    from scipy import sparse

    M = len(states)
    rows, cols, data = [], [], []
    leak = np.zeros(M)
    for i in range(M):
        x = states[i]
        a = propensities(x)
        diag = 0.0
        for j in range(change.shape[0]):
            if a[j] <= 0.0:
                continue
            diag -= a[j]
            key = tuple(int(v) for v in (x + change[j]))
            k = index.get(key)
            if k is not None:
                rows.append(k); cols.append(i); data.append(a[j])
            else:
                leak[i] += a[j]  # probability flux leaving the projection
        rows.append(i); cols.append(i); data.append(diag)
    Q = sparse.csc_matrix((data, (rows, cols)), shape=(M, M))
    return Q, leak


def _kernel_arrays(reactions, species, rate_values, chemostatted_values, volume_L):
    from ._kernels import build_kernel_arrays
    from .stoichiometry import fold_chemostatted_into_rates
    chem_vals = chemostatted_values or {}
    eff = (fold_chemostatted_into_rates(reactions, rate_values, chem_vals)
           if chem_vals else rate_values)
    na_v = AVOGADRO * volume_L
    return build_kernel_arrays(reactions, species, eff, set(chem_vals.keys()), na_v)


def fsp_solve(
    reactions: list[Reaction],
    species: list[str],
    rate_values: dict[str, float],
    initial_conditions: dict[str, float],
    volume_L: float,
    t: float | None = None,
    initial_as: str = "concentration",
    chemostatted_values: dict[str, float] | None = None,
    max_states: int = 100_000,
) -> FSPResult:
    """Solve the CME by finite state projection.

    With ``t`` given, returns the transient distribution ``p(t)``; with ``t=None``,
    returns the stationary distribution (null space of ``Q``) — valid when the
    enumerated state space is closed.
    """
    na_v = AVOGADRO * volume_L
    n0 = np.zeros(len(species), dtype=np.int64)
    for i, s in enumerate(species):
        val = initial_conditions.get(s, 0.0)
        n0[i] = int(round(val)) if initial_as == "count" else int(round(val * na_v))

    c_eff, react_sp, react_co, change = _kernel_arrays(
        reactions, species, rate_values, chemostatted_values, volume_L
    )
    propensities = _propensity_fn(c_eff, react_sp, react_co)
    states, index, truncated = _enumerate_reachable(n0, change, propensities, max_states)
    Q, leak = _build_generator(states, index, change, propensities)
    M = len(states)

    if t is not None:
        from scipy.sparse.linalg import expm_multiply
        p0 = np.zeros(M)
        p0[index[tuple(int(v) for v in n0)]] = 1.0
        p = expm_multiply(Q * float(t), p0)
        p = np.maximum(p, 0.0)
        trunc_err = float(max(0.0, 1.0 - p.sum()))
        return FSPResult(
            species=list(species), states=states, probabilities=p,
            truncation_error=trunc_err, time=float(t), truncated=truncated,
            _index=index,
        )

    # Stationary distribution: solve Qπ = 0 with Σπ = 1 (replace one row with the
    # normalisation constraint).
    A = Q.toarray()
    A[-1, :] = 1.0
    b = np.zeros(M)
    b[-1] = 1.0
    pi = np.linalg.solve(A, b)
    pi = np.maximum(pi, 0.0)
    pi = pi / pi.sum()
    trunc_err = float(leak.max()) if M else 0.0
    return FSPResult(
        species=list(species), states=states, probabilities=pi,
        truncation_error=trunc_err, time=None, truncated=truncated,
        _index=index,
    )
