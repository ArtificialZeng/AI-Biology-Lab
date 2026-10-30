"""Numba-accelerated stochastic kernels with a transparent pure-Python fallback.

The Gillespie SSA and tau-leap inner loops are pure-Python ``for``-loops in
:mod:`mantis.analysis` — fine for correctness, slow at scale.  This module
re-expresses those loops over flat NumPy arrays so they can be JIT-compiled by
``numba`` (optional dependency, ``pip install mantis-delta[fast]``).  When numba
is not installed, the very same functions run as ordinary interpreted Python, so
behaviour is identical and only speed differs.

System encoding (built once by the callers):

* ``c_eff``      : float64[n_rxn]  — stochastic rate constant with the
  ``1/Π coeff!`` combinatorial factor and the ``1/(N_A V)^(order-1)`` volume
  factor already folded in, so a propensity is simply ``c_eff[j] · Π falling``.
* ``react_sp``   : int64[n_rxn, maxR] — reactant species indices, ``-1`` padded.
* ``react_co``   : int64[n_rxn, maxR] — matching stoichiometric coefficients.
* ``change``     : int64[n_rxn, n_sp] — net state change per reaction firing.

Reproducibility uses ``np.random.seed`` so each backend is deterministic for a
given integer seed (the JIT and interpreted streams are independent — tests pin
behaviour per backend, not across backends).
"""
from __future__ import annotations

import math

import numpy as np

try:  # pragma: no cover - trivial import guard
    from numba import njit as _njit
    HAVE_NUMBA = True
except Exception:  # numba absent → identity decorator, plain-Python execution
    HAVE_NUMBA = False

    def _njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]

        def wrap(func):
            return func

        return wrap


def build_kernel_arrays(reactions, species, eff_rates, chem, na_v):
    """Flatten a reaction system into the dense arrays the kernels consume.

    ``eff_rates`` are the deterministic rate constants (already chemostat-folded);
    ``chem`` is the set of chemostatted species names; ``na_v`` is N_A · volume.
    Returns ``(c_eff, react_sp, react_co, change)`` with NumPy integer/float dtypes.
    """
    sp_idx = {s: i for i, s in enumerate(species)}
    n_sp = len(species)
    n_rxn = len(reactions)

    reactant_terms: list[list[tuple[int, int]]] = []
    c_eff = np.zeros(n_rxn, dtype=np.float64)
    change = np.zeros((n_rxn, n_sp), dtype=np.int64)

    for j, rxn in enumerate(reactions):
        k_det = float(eff_rates.get(rxn.rate_key, 0.0))
        terms = [(sp_idx[name], coeff) for name, coeff in rxn.reactants
                 if name in sp_idx and name not in chem]
        reactant_terms.append(terms)
        order = sum(c for _, c in terms)
        c = k_det / (na_v ** (order - 1)) if order > 1 else k_det
        for _, coeff in terms:
            c /= math.factorial(coeff)
        c_eff[j] = c
        for name, coeff in rxn.reactants:
            if name in sp_idx and name not in chem:
                change[j, sp_idx[name]] -= coeff
        for name, coeff in rxn.products:
            if name in sp_idx and name not in chem:
                change[j, sp_idx[name]] += coeff

    maxr = max((len(t) for t in reactant_terms), default=0)
    maxr = max(maxr, 1)
    react_sp = np.full((n_rxn, maxr), -1, dtype=np.int64)
    react_co = np.zeros((n_rxn, maxr), dtype=np.int64)
    for j, terms in enumerate(reactant_terms):
        for r, (idx, coeff) in enumerate(terms):
            react_sp[j, r] = idx
            react_co[j, r] = coeff

    return c_eff, react_sp, react_co, change


def normalize_seed(seed):
    """Coerce a seed (or None) into a 32-bit integer the RNG backends accept."""
    if seed is None:
        seed = int(np.random.SeedSequence().generate_state(1)[0])
    return int(seed) & 0xFFFFFFFF


@_njit(cache=True)
def _propensities(n, c_eff, react_sp, react_co, out):
    """Fill ``out`` with mass-action stochastic propensities for state ``n``."""
    n_rxn = c_eff.shape[0]
    maxr = react_sp.shape[1]
    for j in range(n_rxn):
        if c_eff[j] <= 0.0:
            out[j] = 0.0
            continue
        h = 1.0
        ok = True
        for r in range(maxr):
            idx = react_sp[j, r]
            if idx < 0:
                break
            coeff = react_co[j, r]
            ni = n[idx]
            if ni < coeff:
                ok = False
                break
            ff = 1.0
            for d in range(coeff):
                ff *= (ni - d)
            h *= ff
        out[j] = c_eff[j] * h if ok else 0.0


@_njit(cache=True)
def ssa_kernel(n0, c_eff, react_sp, react_co, change, t0, tf, max_events, seed):
    """Direct-method Gillespie SSA; records the full event-by-event trajectory.

    Returns ``(times, states, n_events, success)`` where ``states`` is
    ``int64[len(times), n_sp]`` including the initial point and a final padding
    point at ``tf``.
    """
    np.random.seed(seed)
    n_rxn = c_eff.shape[0]
    n_sp = n0.shape[0]

    cap = 1024
    times = np.empty(cap, dtype=np.float64)
    states = np.empty((cap, n_sp), dtype=np.int64)
    n = n0.copy()
    times[0] = t0
    for s in range(n_sp):
        states[0, s] = n[s]
    count = 1

    a = np.empty(n_rxn, dtype=np.float64)
    t = t0
    n_events = 0
    success = True

    while t < tf:
        if n_events >= max_events:
            success = False
            break
        _propensities(n, c_eff, react_sp, react_co, a)
        a0 = 0.0
        for j in range(n_rxn):
            a0 += a[j]
        if a0 <= 0.0:
            break
        r1 = np.random.random()
        r2 = np.random.random()
        tau = -np.log(r1) / a0
        t_new = t + tau
        if t_new > tf:
            break
        threshold = r2 * a0
        cum = 0.0
        chosen = n_rxn - 1
        for j in range(n_rxn):
            cum += a[j]
            if cum >= threshold:
                chosen = j
                break
        for s in range(n_sp):
            n[s] += change[chosen, s]
        t = t_new

        if count >= cap:
            cap *= 2
            new_times = np.empty(cap, dtype=np.float64)
            new_states = np.empty((cap, n_sp), dtype=np.int64)
            for i in range(count):
                new_times[i] = times[i]
                for s in range(n_sp):
                    new_states[i, s] = states[i, s]
            times = new_times
            states = new_states
        times[count] = t
        for s in range(n_sp):
            states[count, s] = n[s]
        count += 1
        n_events += 1

    # Pad to tf so the trajectory always spans the full interval.
    if times[count - 1] < tf:
        if count >= cap:
            cap += 1
            new_times = np.empty(cap, dtype=np.float64)
            new_states = np.empty((cap, n_sp), dtype=np.int64)
            for i in range(count):
                new_times[i] = times[i]
                for s in range(n_sp):
                    new_states[i, s] = states[i, s]
            times = new_times
            states = new_states
        times[count] = tf
        for s in range(n_sp):
            states[count, s] = n[s]
        count += 1

    return times[:count].copy(), states[:count].copy(), n_events, success


@_njit(cache=True)
def tau_leap_kernel(n0, c_eff, react_sp, react_co, change, t0, tf,
                    fixed_tau, epsilon, record_times, seed):
    """Adaptive (or fixed-step) tau-leap; records at ``record_times`` snapshots.

    Returns ``(record_states, success)`` of shape ``int64[len(record_times), n_sp]``.
    ``fixed_tau <= 0`` selects the adaptive Cao-2006 step; otherwise a fixed leap.
    """
    np.random.seed(seed)
    n_rxn = c_eff.shape[0]
    n_sp = n0.shape[0]
    n_rec = record_times.shape[0]

    record_states = np.zeros((n_rec, n_sp), dtype=np.int64)
    n = n0.copy()
    for s in range(n_sp):
        record_states[0, s] = n[s]
    next_record = 1
    success = True

    a = np.empty(n_rxn, dtype=np.float64)
    t = t0
    while t < tf and next_record < n_rec:
        _propensities(n, c_eff, react_sp, react_co, a)
        a0 = 0.0
        for j in range(n_rxn):
            a0 += a[j]
        if a0 <= 0.0:
            break

        if fixed_tau > 0.0:
            dt = fixed_tau
        else:
            # Cao et al. 2006 leap-size bound (g_i = 2 approximation, ≤2nd order).
            dt = np.inf
            for s in range(n_sp):
                mu = 0.0
                sigma2 = 0.0
                for j in range(n_rxn):
                    nu = change[j, s]
                    mu += nu * a[j]
                    sigma2 += nu * nu * a[j]
                xi = n[s]
                if xi < 1:
                    xi = 1
                bound = epsilon * xi / 2.0
                if bound < 1.0:
                    bound = 1.0
                if mu < 0.0:
                    mu = -mu
                if mu > 0.0:
                    cand = bound / mu
                    if cand < dt:
                        dt = cand
                if sigma2 > 0.0:
                    cand = (bound * bound) / sigma2
                    if cand < dt:
                        dt = cand

        if dt > tf - t:
            dt = tf - t
        if dt > record_times[next_record] - t:
            dt = record_times[next_record] - t
        if dt <= 0.0:
            dt = tf - t

        # Provisional leap; reject into one exact SSA step if it goes negative.
        negative = False
        new_n = n.copy()
        for j in range(n_rxn):
            k = np.random.poisson(a[j] * dt)
            if k != 0:
                for s in range(n_sp):
                    new_n[s] += k * change[j, s]
        for s in range(n_sp):
            if new_n[s] < 0:
                negative = True
                break

        if negative:
            r1 = np.random.random()
            r2 = np.random.random()
            ssa_dt = -np.log(r1 if r1 > 1e-300 else 1e-300) / a0
            if t + ssa_dt > tf:
                t = tf
                break
            threshold = r2 * a0
            cum = 0.0
            chosen = n_rxn - 1
            for j in range(n_rxn):
                cum += a[j]
                if cum >= threshold:
                    chosen = j
                    break
            for s in range(n_sp):
                n[s] += change[chosen, s]
            t = t + ssa_dt
        else:
            for s in range(n_sp):
                n[s] = new_n[s]
            t = t + dt

        while next_record < n_rec and record_times[next_record] <= t:
            for s in range(n_sp):
                record_states[next_record, s] = n[s]
            next_record += 1

    for i in range(next_record, n_rec):
        for s in range(n_sp):
            record_states[i, s] = n[s]

    return record_states, success
