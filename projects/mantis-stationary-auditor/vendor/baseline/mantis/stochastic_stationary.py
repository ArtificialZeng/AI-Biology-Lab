"""Exact stochastic stationary distributions — Anderson, Craciun & Kurtz (*Bull.
Math. Biol.* 2010).

For a **complex-balanced** mass-action network the stochastic (chemical master
equation) model has an *exact* stationary distribution in closed form — a product
of independent Poissons whose means are the deterministic complex-balanced
equilibrium scaled to molecule counts:

    π(x) = M · Π_i  λ_i^{x_i} / x_i!,        λ_i = (N_A · V) · c_i*,

supported on the reachable state space (the stoichiometric compatibility class).
With no conservation laws the state space is all of ℤ≥0ⁿ and π factorises into
independent ``Poisson(λ_i)``; with conservation laws π is that product *conditioned*
on the conserved totals (for a single moiety this is exactly a Binomial /
multinomial).  This lets you write the stationary distribution down **from network
structure alone — no simulation** — for any weakly reversible deficiency-zero
network, the regime where SSA is noisiest.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .analysis import AVOGADRO
from .parsing import Reaction


def _enumerate_class(cons_vectors, cons_totals, n_species, max_states):
    """Enumerate non-negative integer states satisfying all conservation constraints.

    Returns ``(states, exact)``: an ``(M, n)`` int array and whether enumeration was
    complete (``exact=False`` if a species is unconstrained or the cap was hit).
    """
    V = [np.rint(v).astype(np.int64) for v in cons_vectors]
    T = [int(round(t)) for t in cons_totals]
    # Every species must appear in at least one conservation law, else the class is
    # infinite and cannot be enumerated.
    covered = np.zeros(n_species, dtype=bool)
    for v in V:
        covered |= v > 0
    if not covered.all():
        return None, False

    states: list[list[int]] = []
    partial = [0] * n_species
    remaining = list(T)

    def upper_bound(i, remaining):
        ub = None
        for li, v in enumerate(V):
            if v[i] > 0:
                b = remaining[li] // v[i]
                ub = b if ub is None else min(ub, b)
        return 0 if ub is None else max(ub, 0)

    def dfs(i):
        if len(states) > max_states:
            return False
        if i == n_species:
            if all(r == 0 for r in remaining):
                states.append(list(partial))
            return True
        ub = upper_bound(i, remaining)
        for val in range(ub + 1):
            partial[i] = val
            for li, v in enumerate(V):
                remaining[li] -= v[i] * val
            ok = dfs(i + 1)
            for li, v in enumerate(V):
                remaining[li] += v[i] * val
            if not ok:
                return False
        partial[i] = 0
        return True

    complete = dfs(0)
    if not complete or len(states) > max_states:
        return None, False
    return np.array(states, dtype=np.int64).reshape(-1, n_species), True


@dataclass
class StationaryDistribution:
    """Exact product-form stationary distribution of a complex-balanced network.

    Attributes
    ----------
    species : list[str]
    means : np.ndarray
        The Poisson means ``λ_i`` (in molecule counts) = ``N_A · V · c_i*``.
    volume_L : float
    is_conditioned : bool
        Whether conservation laws restrict the support to a compatibility class.
    is_exact : bool
        True when the distribution is represented exactly (independent Poissons, or a
        fully enumerated finite class); False when the class could not be enumerated
        and the independent-Poisson factorisation is used as an approximation.
    """
    species: list[str]
    means: np.ndarray
    volume_L: float
    is_conditioned: bool
    is_exact: bool
    _states: np.ndarray | None = field(default=None, repr=False)
    _logZ: float | None = field(default=None, repr=False)

    # ── core ────────────────────────────────────────────────────────────────────
    def poisson_means(self) -> dict[str, float]:
        """Poisson means ``λ_i`` keyed by species (molecule counts)."""
        return {s: float(m) for s, m in zip(self.species, self.means)}

    def _log_weight(self, x: np.ndarray) -> float:
        from scipy.special import gammaln
        x = np.asarray(x, dtype=float)
        lam = np.maximum(self.means, 1e-300)
        return float(np.sum(x * np.log(lam) - gammaln(x + 1.0)))

    def probability(self, counts: dict[str, float]) -> float:
        """Normalised stationary probability of an integer state ``counts``."""
        from scipy.stats import poisson
        x = np.array([counts.get(s, 0) for s in self.species], dtype=np.int64)
        if not self.is_conditioned:
            return float(np.prod([poisson.pmf(xi, li)
                                  for xi, li in zip(x, self.means)]))
        if self._states is None:
            # Could not enumerate the class → independent-Poisson approximation.
            return float(np.prod([poisson.pmf(xi, li)
                                  for xi, li in zip(x, self.means)]))
        from scipy.special import logsumexp
        if self._logZ is None:
            self._logZ = float(logsumexp([self._log_weight(s) for s in self._states]))
        return float(np.exp(self._log_weight(x) - self._logZ))

    def expected_counts(self) -> dict[str, float]:
        """Stationary mean molecule counts (conditional means when conditioned)."""
        if not self.is_conditioned or self._states is None:
            return {s: float(m) for s, m in zip(self.species, self.means)}
        probs = np.array([self.probability(dict(zip(self.species, st)))
                          for st in self._states])
        mean = (self._states * probs[:, None]).sum(axis=0)
        return {s: float(m) for s, m in zip(self.species, mean)}

    def marginal(self, species: str) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(values, probabilities)`` for the marginal of one species."""
        i = self.species.index(species)
        if not self.is_conditioned or self._states is None:
            from scipy.stats import poisson
            hi = int(poisson.ppf(0.9999, self.means[i])) + 1
            vals = np.arange(0, hi + 1)
            return vals, poisson.pmf(vals, self.means[i])
        probs = np.array([self.probability(dict(zip(self.species, st)))
                          for st in self._states])
        col = self._states[:, i]
        vals = np.arange(col.min(), col.max() + 1)
        out = np.zeros(len(vals))
        for v, p in zip(col, probs):
            out[v - col.min()] += p
        return vals, out

    def sample(self, size: int = 1, seed: int | None = None) -> dict[str, np.ndarray]:
        """Draw samples from the stationary distribution."""
        rng = np.random.default_rng(seed)
        if not self.is_conditioned or self._states is None:
            draws = np.array([rng.poisson(m, size=size) for m in self.means]).T
        else:
            probs = np.array([self.probability(dict(zip(self.species, st)))
                              for st in self._states])
            probs = probs / probs.sum()
            idx = rng.choice(len(self._states), size=size, p=probs)
            draws = self._states[idx]
        return {s: draws[:, i] for i, s in enumerate(self.species)}


def stationary_distribution(
    reactions: list[Reaction],
    species: list[str],
    rate_values: dict[str, float],
    deficiency: int,
    is_weakly_reversible: bool,
    initial_conditions: dict[str, float],
    volume_L: float,
    initial_as: str = "concentration",
    chemostatted_values: dict[str, float] | None = None,
    max_states: int = 200_000,
) -> StationaryDistribution:
    """Build the exact ACK product-form stationary distribution.

    Raises ``ValueError`` if the network is not complex-balanced (the closed-form
    result does not apply).
    """
    from .stability import is_complex_balanced, complex_balanced_equilibrium
    from .analysis import _conservation_law_vectors

    chem = chemostatted_values or {}
    NA_V = AVOGADRO * volume_L

    # Initial molecule counts and the equivalent concentrations.
    n0 = np.zeros(len(species))
    for i, s in enumerate(species):
        val = initial_conditions.get(s, 0.0)
        n0[i] = round(val) if initial_as == "count" else round(val * NA_V)
    ic_conc = {s: n0[i] / NA_V for i, s in enumerate(species)}

    if not is_complex_balanced(reactions, species, rate_values, deficiency,
                               is_weakly_reversible, ic_conc, chem or None):
        raise ValueError(
            "Network is not complex-balanced at these rate constants, so the "
            "Anderson–Craciun–Kurtz product-form stationary distribution does not "
            "apply. (It holds for every weakly reversible deficiency-zero network, "
            "and for higher-deficiency networks only at complex-balanced rates.)"
        )

    c_star = complex_balanced_equilibrium(
        reactions, species, rate_values, deficiency, is_weakly_reversible,
        ic_conc, chem or None,
    )
    means = np.array([c_star[s] * NA_V for s in species])

    cons_vectors = _conservation_law_vectors(reactions, species)
    is_conditioned = len(cons_vectors) > 0
    states = None
    is_exact = True
    if is_conditioned:
        cons_totals = [float(v @ n0) for v in cons_vectors]
        states, exact = _enumerate_class(cons_vectors, cons_totals, len(species), max_states)
        is_exact = exact

    return StationaryDistribution(
        species=list(species),
        means=means,
        volume_L=volume_L,
        is_conditioned=is_conditioned,
        is_exact=is_exact,
        _states=states,
    )
