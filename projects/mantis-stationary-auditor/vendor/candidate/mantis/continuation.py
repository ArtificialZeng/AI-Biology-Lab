"""Pseudo-arclength continuation of steady-state branches (fold / Hopf detection).

The :meth:`~mantis.network.CRNetwork.bifurcation` log-scan recomputes steady
states *independently* at each parameter value.  Near a saddle-node (fold) the
branch turns back on itself, so two states share one parameter value and one
disappears as the parameter is swept past the fold — a naive scan loses the
branch, mis-pairs points across the fold, and never reports the fold itself.

**Pseudo-arclength continuation** (Keller 1977; Allgower–Georg 2003) fixes this
by parameterising the branch not by the bifurcation parameter ``λ`` but by
*arclength* ``s`` along the solution curve ``(c(s), λ(s))``.  Folds — where
``dλ/ds = 0`` — are then ordinary points of the curve and are traversed
smoothly.  The method is the standard predictor–corrector:

1. **Predictor** — step ``ds`` along the unit tangent of the branch.
2. **Corrector** — Newton on the *augmented* square system

       G(c, λ)                                  = 0     (n steady-state equations)
       ⟨t_c, c − c_prev⟩ + t_λ (λ − λ_prev) − ds = 0     (arclength constraint)

   The arclength row keeps the Jacobian non-singular *through* the fold, where
   the steady-state Jacobian ``G_c`` alone is singular.

``G`` is the conservation-constrained steady-state map: ``rank(N)`` independent
rows of the mass-action ODE field stacked on the ``n − rank(N)`` conservation
laws (totals fixed by the initial conditions), so the branch stays inside one
stoichiometric compatibility class.

**Bifurcation detection** along the branch:

* **Fold (saddle-node)** — ``det G_c`` changes sign (the steady-state Jacobian
  is singular exactly at a fold), the canonical fold test function.
* **Hopf** — the number of eigenvalues with positive real part changes *without*
  a fold (so a complex-conjugate pair, not a real eigenvalue, has crossed the
  imaginary axis).  Reported as a candidate; a Hopf marks the birth of a limit
  cycle.  Eigenvalues are those of the Jacobian restricted to the stoichiometric
  subspace (the ``n − rank(N)`` exact zeros from conservation are dropped).

The classic demonstration is the Schlögl model: continuation traces the full
S-shaped branch and pins both fold points bounding the bistable region, neither
of which the log-scan reports.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .parsing import Reaction


@dataclass
class BifurcationPoint:
    """A detected bifurcation on a continuation branch.

    Attributes
    ----------
    kind : str
        ``"fold"`` (saddle-node) or ``"hopf"`` (candidate Hopf).
    parameter : float
        Bifurcation-parameter value at the point (linearly localised between the
        two bracketing continuation steps).
    state : dict[str, float]
        Species concentrations at the bifurcation.
    eigenvalue : complex
        The critical eigenvalue nearest the imaginary axis at the crossing.
    """
    kind: str
    parameter: float
    state: dict[str, float]
    eigenvalue: complex

    def __str__(self) -> str:
        loc = ", ".join(f"{s}={v:.4g}" for s, v in self.state.items())
        return (f"{self.kind.upper()} at {self.parameter:.6g} "
                f"(λ_crit={self.eigenvalue:.3g}); {loc}")


@dataclass
class ContinuationResult:
    """A steady-state branch traced by pseudo-arclength continuation.

    The branch is ordered by arclength, so ``parameter_values`` is generally
    *non-monotonic* (it reverses at each fold) — that is the whole point.

    Attributes
    ----------
    parameter_name : str
        Rate key that was continued.
    parameter_values : np.ndarray
        Bifurcation-parameter value at each branch point (arclength order).
    branch : list[dict[str, float]]
        Species concentrations at each branch point.
    stable : np.ndarray
        Boolean stability of each branch point (all reduced eigenvalues Re < 0).
    eigenvalues : list[np.ndarray]
        Reduced-Jacobian eigenvalues at each branch point.
    bifurcations : list[BifurcationPoint]
    """
    parameter_name: str
    parameter_values: np.ndarray
    branch: list[dict[str, float]]
    stable: np.ndarray
    eigenvalues: list[np.ndarray] = field(default_factory=list)
    bifurcations: list[BifurcationPoint] = field(default_factory=list)

    def folds(self) -> list[BifurcationPoint]:
        return [b for b in self.bifurcations if b.kind == "fold"]

    def hopfs(self) -> list[BifurcationPoint]:
        return [b for b in self.bifurcations if b.kind == "hopf"]

    def species_branch(self, species: str) -> np.ndarray:
        """The concentration of one species along the branch (arclength order)."""
        return np.array([pt[species] for pt in self.branch])

    def __str__(self) -> str:
        lines = [
            f"Continuation in {self.parameter_name!r}: {len(self.branch)} branch "
            f"points over λ ∈ [{self.parameter_values.min():.4g}, "
            f"{self.parameter_values.max():.4g}]",
        ]
        if self.bifurcations:
            lines.append(f"  {len(self.bifurcations)} bifurcation(s):")
            lines += [f"    • {b}" for b in self.bifurcations]
        else:
            lines.append("  no bifurcations detected on this branch")
        return "\n".join(lines)


class _SteadyStateMap:
    """The conservation-constrained steady-state map ``G(c, λ)`` and its derivatives.

    ``λ`` is the value of one rate constant (``param_key``); every evaluation
    folds it (and any chemostatted concentrations) into the effective rates and
    reuses the tested ODE / Jacobian builders from :mod:`mantis.analysis`.
    """

    def __init__(self, reactions, species, base_rates, param_key,
                 chem_vals, initial_conditions):
        from scipy.linalg import qr
        from .stoichiometry import (
            build_stoichiometry_matrix, matrix_rank, fold_chemostatted_into_rates,
        )
        from .analysis import _conservation_law_vectors

        self.reactions = reactions
        self.species = species
        self.base_rates = dict(base_rates)
        self.param_key = param_key
        self.chem_vals = chem_vals or {}
        self.chem_keys = set(self.chem_vals)
        self._fold = fold_chemostatted_into_rates

        self.N = build_stoichiometry_matrix(reactions, species)
        self.rank = matrix_rank(self.N)
        if self.rank > 0:
            _, _, perm = qr(self.N.T, pivoting=True)
            self.ind = list(perm[: self.rank])
        else:
            self.ind = []

        self.cl_vecs = _conservation_law_vectors(reactions, species)
        y_ref = np.array([initial_conditions.get(s, 0.0) for s in species])
        self.cl_totals = np.array([float(v @ y_ref) for v in self.cl_vecs])

        # Index of the reaction whose rate constant is the continuation parameter.
        self.p = next(j for j, r in enumerate(reactions) if r.rate_key == param_key)
        sp_idx = {s: i for i, s in enumerate(species)}
        self.p_reactants = [
            (sp_idx[name], coeff)
            for name, coeff in reactions[self.p].reactants
            if name not in self.chem_keys
        ]
        self.n = len(species)

    def _eff_rates(self, lam: float) -> dict[str, float]:
        rates = dict(self.base_rates)
        rates[self.param_key] = lam
        if self.chem_vals:
            return self._fold(self.reactions, rates, self.chem_vals)
        return rates

    def G(self, c: np.ndarray, lam: float) -> np.ndarray:
        from .analysis import build_ode_function
        f = build_ode_function(self.reactions, self.species, self._eff_rates(lam),
                               self.chem_keys)
        ode = f(0.0, c)[self.ind] if self.ind else np.zeros(0)
        cons = np.array([v @ c - tot for v, tot in zip(self.cl_vecs, self.cl_totals)])
        return np.concatenate([ode, cons])

    def G_c(self, c: np.ndarray, lam: float) -> np.ndarray:
        """∂G/∂c — independent ODE-Jacobian rows stacked on the conservation rows."""
        from .analysis import _full_jacobian_fn
        J = _full_jacobian_fn(self.reactions, self.species, self._eff_rates(lam),
                              self.chem_keys)(c)
        top = J[self.ind, :] if self.ind else np.zeros((0, self.n))
        bottom = np.array(self.cl_vecs) if self.cl_vecs else np.zeros((0, self.n))
        return np.vstack([top, bottom]) if (len(top) + len(bottom)) else np.zeros((0, self.n))

    def G_lam(self, c: np.ndarray, lam: float) -> np.ndarray:
        """∂G/∂λ — only the chosen reaction's flux depends (linearly) on λ."""
        eff = self._eff_rates(lam)
        flux_p = eff[self.param_key]
        for idx, coeff in self.p_reactants:
            flux_p *= max(c[idx], 0.0) ** coeff
        dode = self.N[self.ind, self.p] * (flux_p / lam) if (self.ind and lam != 0) \
            else np.zeros(len(self.ind))
        return np.concatenate([dode, np.zeros(len(self.cl_vecs))])

    def reduced_eigenvalues(self, c: np.ndarray, lam: float) -> np.ndarray:
        """Eigenvalues of the Jacobian restricted to the stoichiometric subspace.

        The full mass-action Jacobian carries ``n − rank(N)`` exact zero
        eigenvalues (one per conservation law, since ``wᵀN = 0 ⟹ wᵀJ = 0``);
        those are dropped, leaving the ``rank(N)`` dynamically relevant ones.
        """
        from scipy.linalg import eigvals
        from .analysis import _full_jacobian_fn
        J = _full_jacobian_fn(self.reactions, self.species, self._eff_rates(lam),
                              self.chem_keys)(c)
        eigs = eigvals(J)
        n_drop = self.n - self.rank
        if n_drop > 0:
            keep = np.argsort(np.abs(eigs))[n_drop:]
            eigs = eigs[keep]
        return eigs


def _tangent(Gc: np.ndarray, Glam: np.ndarray, prev: np.ndarray | None) -> np.ndarray:
    """Unit tangent of the branch: the null vector of [G_c | G_λ], oriented forward."""
    A = np.hstack([Gc, Glam.reshape(-1, 1)])
    # Right singular vector of the smallest singular value spans the null space.
    _, _, Vt = np.linalg.svd(A)
    t = Vt[-1]
    t = t / np.linalg.norm(t)
    if prev is not None:
        if np.dot(t, prev) < 0:          # keep a consistent travel direction (through folds)
            t = -t
    elif t[-1] < 0:                       # first step: head toward increasing λ
        t = -t
    return t


def _critical_eigenvalue(eigs: np.ndarray) -> complex:
    """The eigenvalue nearest the imaginary axis (smallest |Re|)."""
    return complex(eigs[np.argmin(np.abs(eigs.real))])


def pseudo_arclength_continuation(
    reactions: list[Reaction],
    species: list[str],
    base_rates: dict[str, float],
    parameter: str,
    param_range: tuple[float, float],
    initial_conditions: dict[str, float],
    chemostatted_values: dict[str, float] | None = None,
    ds: float = 0.05,
    max_steps: int = 2000,
    newton_tol: float = 1e-10,
    initial_state: dict[str, float] | None = None,
) -> ContinuationResult:
    """Trace a steady-state branch by pseudo-arclength continuation.

    Parameters
    ----------
    parameter : str
        Rate key to continue (e.g. ``"A + 2X -> 3X"``).
    param_range : (float, float)
        ``(λ_min, λ_max)`` box the branch is traced within; continuation stops
        when it leaves the box.  To capture an S-curve the box must bracket both
        folds.
    initial_conditions : dict[str, float]
        Fixes the conservation totals (the compatibility class) and, unless
        ``initial_state`` is given, seeds the branch (a steady state is solved at
        ``λ_min`` to start from).
    ds : float
        Arclength step (in the non-dimensional ``(c/scale, λ/λ_min)`` metric);
        adapted down on Newton failure and back up on easy steps.
    max_steps : int
        Cap on continuation steps (guards against closed isolas / cycling).
    initial_state : dict[str, float], optional
        Explicit starting steady state on the desired branch; overrides the
        automatic seed (use to select which branch to follow).

    Returns
    -------
    ContinuationResult
    """
    from .parsing import normalize_rate_key
    from .analysis import find_steady_states, classify_steady_state
    from .steady_states import all_steady_states

    try:
        param_key = normalize_rate_key(parameter)
    except ValueError:
        param_key = parameter

    lam_min, lam_max = float(param_range[0]), float(param_range[1])
    chem_vals = chemostatted_values or {}
    smap = _SteadyStateMap(reactions, species, base_rates, param_key,
                           chem_vals, initial_conditions)
    n = smap.n

    # ── Seed: a steady state on the branch at λ_min ──────────────────────────
    if initial_state is not None:
        c = np.array([initial_state.get(s, 0.0) for s in species], float)
    else:
        seed_rates = dict(base_rates)
        seed_rates[param_key] = lam_min
        # Exhaustive enumeration gives a reliable seed (the multi-start solver can
        # converge to spurious points on open systems); fall back to it otherwise.
        ss = []
        try:
            ss = all_steady_states(reactions, species, seed_rates, initial_conditions,
                                   chemostatted_values=chem_vals or None)
        except Exception:
            ss = []
        if not ss:
            ss = find_steady_states(reactions, species, seed_rates, initial_conditions,
                                    n_attempts=20, seed=0,
                                    chemostatted_values=chem_vals or None)
        if not ss:
            raise RuntimeError("Could not find a steady state to seed continuation.")
        c = np.array([ss[0].concentrations[s] for s in species], float)
    lam = lam_min

    # Non-dimensionalising scales so c and λ contribute comparably to arclength.
    cscale = max(float(np.max(np.abs(c))), 1.0)
    lscale = max(abs(lam_min), abs(lam_max), 1e-300)

    def scaled(cv, lv):
        return np.concatenate([cv / cscale, [lv / lscale]])

    def unscaled(x):
        return x[:n] * cscale, x[n] * lscale

    # ── Record helpers ───────────────────────────────────────────────────────
    params: list[float] = []
    branch: list[dict[str, float]] = []
    stable: list[bool] = []
    eig_list: list[np.ndarray] = []
    det_signs: list[float] = []
    n_unstable: list[int] = []
    bifs: list[BifurcationPoint] = []

    xs_scaled: list[np.ndarray] = []      # accepted points in scaled coordinates
    tans: list[np.ndarray] = []           # unit tangent at each accepted point

    def record(cv, lv, tan):
        eigs = smap.reduced_eigenvalues(cv, lv)
        is_stable, _ = classify_steady_state(eigs)
        params.append(lv)
        branch.append(dict(zip(species, cv.tolist())))
        stable.append(is_stable)
        eig_list.append(eigs)
        det = float(np.linalg.det(smap.G_c(cv, lv)))
        det_signs.append(np.sign(det) if det != 0 else (det_signs[-1] if det_signs else 1.0))
        n_unstable.append(int(np.sum(eigs.real > 1e-9 * (np.max(np.abs(eigs)) + 1e-300))))
        xs_scaled.append(scaled(cv, lv))
        tans.append(tan)

    record(c, lam, None)
    prev_tan = None

    for _ in range(max_steps):
        Gc = smap.G_c(c, lam)
        Glam = smap.G_lam(c, lam)
        # Scale derivatives consistently with the (c/cscale, λ/lscale) variables.
        tan = _tangent(Gc * cscale, Glam * lscale, prev_tan)

        x0 = scaled(c, lam)
        x_pred = x0 + ds * tan

        x, converged = _correct(smap, n, cscale, lscale, x_pred, tan, x0, ds, newton_tol)

        if not converged:
            ds *= 0.5
            if ds < 1e-9:
                break
            continue

        cv, lv = unscaled(x)
        if np.any(cv < -1e-7 * (np.max(np.abs(cv)) + 1e-30)):
            break                                   # left the positive orthant
        cv = np.maximum(cv, 0.0)

        if tans[-1] is None:
            tans[-1] = tan                          # tangent leaving the seed point
        record(cv, lv, tan)

        c, lam = cv, lv
        prev_tan = tan
        ds = min(ds * 1.1, 0.1)                      # grow step on success

        if lam < lam_min or lam > lam_max:
            break

    _detect_bifurcations(smap, species, n, cscale, lscale, xs_scaled, tans,
                         params, det_signs, n_unstable, bifs)

    return ContinuationResult(
        parameter_name=param_key,
        parameter_values=np.array(params),
        branch=branch,
        stable=np.array(stable),
        eigenvalues=eig_list,
        bifurcations=bifs,
    )


def _correct(smap, n, cscale, lscale, x_pred, tan, x0, ds, newton_tol, max_iter=50):
    """Newton-correct ``x_pred`` onto the branch under the pseudo-arclength row.

    The constraint ``⟨tan, x − x0⟩ = ds`` keeps the augmented Jacobian regular
    through folds, where the steady-state block ``G_c`` alone is singular.
    """
    x = x_pred.copy()
    for _it in range(max_iter):
        cv, lv = x[:n] * cscale, x[n] * lscale
        g = smap.G(cv, lv)
        arc = float(np.dot(tan, x - x0) - ds)
        res = np.concatenate([g, [arc]])
        if np.linalg.norm(res) < newton_tol * (1 + np.linalg.norm(g) + 1):
            return x, True
        Jc = smap.G_c(cv, lv) * cscale
        Jl = (smap.G_lam(cv, lv) * lscale).reshape(-1, 1)
        A = np.vstack([np.hstack([Jc, Jl]), tan.reshape(1, -1)])
        try:
            dx = np.linalg.solve(A, res)
        except np.linalg.LinAlgError:
            return x, False
        x = x - dx
    return x, False


def _det_scaled(smap, n, cscale, lscale, x):
    """``det G_c`` (the fold test function) at a scaled point ``x``."""
    return float(np.linalg.det(smap.G_c(x[:n] * cscale, x[n] * lscale)))


def _refine_fold(smap, n, cscale, lscale, xa, xb, ta, newton_tol):
    """Locate the fold (``det G_c = 0``) between bracketing scaled points ``xa``/``xb``.

    Bisects the chord parameter; at each trial it re-corrects the predictor *onto
    the branch* (orthogonal to the local tangent ``ta``) before evaluating the
    test function, so the located point is a genuine steady state, not a chord
    interpolant.
    """
    da = _det_scaled(smap, n, cscale, lscale, xa)
    lo, hi = 0.0, 1.0
    x_lo, x_hi = xa, xb
    best = 0.5 * (xa + xb)
    for _ in range(40):
        t = 0.5 * (lo + hi)
        pred = (1 - t) * xa + t * xb
        # Correct transversally onto the branch (anchor the arclength row at pred,
        # ds = 0 ⇒ project orthogonally to the tangent).
        x, ok = _correct(smap, n, cscale, lscale, pred, ta, pred, 0.0, newton_tol)
        if not ok:
            x = pred
        best = x
        dt = _det_scaled(smap, n, cscale, lscale, x)
        if abs(dt) < 1e-12:
            break
        if np.sign(dt) == np.sign(da):
            lo = t
        else:
            hi = t
    return best


def _detect_bifurcations(smap, species, n, cscale, lscale, xs_scaled, tans,
                         params, det_signs, n_unstable, bifs):
    """Scan consecutive branch points for fold (det sign flip) and Hopf crossings."""
    for j in range(1, len(params)):
        i = j - 1
        if det_signs[i] * det_signs[j] < 0:
            xstar = _refine_fold(smap, n, cscale, lscale,
                                 xs_scaled[i], xs_scaled[j], tans[j], 1e-11)
            cstar, pstar = np.maximum(xstar[:n] * cscale, 0.0), xstar[n] * lscale
            eigs = smap.reduced_eigenvalues(cstar, pstar)
            bifs.append(BifurcationPoint("fold", float(pstar),
                                         dict(zip(species, cstar.tolist())),
                                         _critical_eigenvalue(eigs)))
        elif n_unstable[i] != n_unstable[j]:
            # Stability changed without a fold ⇒ a complex pair crossed ⇒ Hopf.
            cstar = 0.5 * (xs_scaled[i][:n] + xs_scaled[j][:n]) * cscale
            pstar = 0.5 * (params[i] + params[j])
            eigs = smap.reduced_eigenvalues(cstar, pstar)
            crit = _critical_eigenvalue(eigs)
            if abs(crit.imag) > 1e-9 * (abs(crit) + 1e-300):
                bifs.append(BifurcationPoint("hopf", float(pstar),
                                             dict(zip(species, cstar.tolist())), crit))
