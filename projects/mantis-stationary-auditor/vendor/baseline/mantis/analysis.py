"""Numerical ODE simulation, steady-state finding, stability analysis, and bifurcation scanning."""
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.optimize import least_squares
from scipy.linalg import eigvals

from .parsing import Reaction
from .stoichiometry import (
    build_stoichiometry_matrix,
    conservation_laws_sympy,
    matrix_rank,
)


@dataclass
class SteadyState:
    concentrations: dict[str, float]
    eigenvalues: np.ndarray
    is_stable: bool
    is_oscillatory: bool
    residual: float


@dataclass
class BifurcationResult:
    parameter_name: str
    parameter_values: np.ndarray
    steady_states: list[list[SteadyState]]


@dataclass
class SimulationResult:
    times: np.ndarray
    concentrations: dict[str, np.ndarray]
    success: bool

    def at(self, t: float) -> dict[str, float]:
        """Return interpolated concentrations at time t."""
        idx = int(np.searchsorted(self.times, t, side="right")) - 1
        idx = max(0, min(idx, len(self.times) - 1))
        return {sp: float(arr[idx]) for sp, arr in self.concentrations.items()}

    def final(self) -> dict[str, float]:
        """Return concentrations at the last time point."""
        return {sp: float(arr[-1]) for sp, arr in self.concentrations.items()}


# Avogadro's number — converts between molarity and molecule count.
AVOGADRO = 6.02214076e23


@dataclass
class StochasticResult:
    """Trajectory from a single Gillespie SSA realization.

    Attributes
    ----------
    times          : 1-D ndarray of reaction times (length n_events + 1)
    counts         : dict species → 1-D ndarray of integer molecule counts
    concentrations : dict species → 1-D ndarray of concentrations (M)
                     (counts / (volume * Avogadro))
    n_events       : number of reaction firings recorded
    success        : True if integration finished before exhausting events
    volume_L       : reaction volume in liters
    """
    times: np.ndarray
    counts: dict[str, np.ndarray]
    concentrations: dict[str, np.ndarray]
    n_events: int
    success: bool
    volume_L: float

    def at(self, t: float) -> dict[str, float]:
        """Step-function lookup of concentrations at time t."""
        idx = int(np.searchsorted(self.times, t, side="right")) - 1
        idx = max(0, min(idx, len(self.times) - 1))
        return {sp: float(arr[idx]) for sp, arr in self.concentrations.items()}

    def final(self) -> dict[str, float]:
        return {sp: float(arr[-1]) for sp, arr in self.concentrations.items()}


def gillespie_simulate(
    reactions: list[Reaction],
    species: list[str],
    rate_values: dict[str, float],
    initial_conditions: dict[str, float],
    t_span: tuple[float, float],
    volume_L: float,
    *,
    initial_as: str = "concentration",
    max_events: int = 1_000_000,
    seed: int | None = None,
    chemostatted_values: dict[str, float] | None = None,
) -> StochasticResult:
    """
    Single-trajectory Gillespie direct-method SSA.

    Parameters
    ----------
    reactions, species, rate_values
        Same shape as :func:`simulate_ode` — deterministic mass-action rates.
    initial_conditions
        Either molecule counts (set ``initial_as='count'``) or molar
        concentrations (default).
    t_span : (t0, tf)
        Simulation interval.
    volume_L : float
        Reaction volume in liters.  Required because stochastic propensities
        depend on absolute molecule counts.  For 100 µL physiological volume,
        pass ``1e-4``.
    initial_as : 'concentration' or 'count'
    max_events : safety cap on reaction firings.
    seed : optional RNG seed for reproducibility.
    chemostatted_values : species kept at constant concentration; folded into
        the propensities of reactions that consume them.

    Returns
    -------
    StochasticResult
        Full trajectory of times, counts and concentrations.
    """
    from .stoichiometry import fold_chemostatted_into_rates
    from ._kernels import build_kernel_arrays, normalize_seed, ssa_kernel

    chem_vals = chemostatted_values or {}
    eff_rates = (
        fold_chemostatted_into_rates(reactions, rate_values, chem_vals)
        if chem_vals else rate_values
    )
    chem = set(chem_vals.keys())

    sp_idx = {s: i for i, s in enumerate(species)}
    n_sp = len(species)
    NA_V = AVOGADRO * volume_L

    # Initial counts
    n0 = np.zeros(n_sp, dtype=np.int64)
    for s, val in initial_conditions.items():
        if s not in sp_idx:
            continue
        if initial_as == "count":
            n0[sp_idx[s]] = int(round(val))
        else:
            n0[sp_idx[s]] = int(round(val * NA_V))

    c_eff, react_sp, react_co, change = build_kernel_arrays(
        reactions, species, eff_rates, chem, NA_V
    )

    # The direct-method SSA loop runs in a numba-JIT kernel when available,
    # else as an equivalent pure-Python loop (see mantis._kernels).
    times_arr, history_arr, n_events, success = ssa_kernel(
        n0, c_eff, react_sp, react_co, change,
        float(t_span[0]), float(t_span[1]), int(max_events), normalize_seed(seed),
    )

    counts_dict = {sp: history_arr[:, i] for i, sp in enumerate(species)}
    conc_dict = {sp: history_arr[:, i] / NA_V for i, sp in enumerate(species)}

    return StochasticResult(
        times=times_arr,
        counts=counts_dict,
        concentrations=conc_dict,
        n_events=int(n_events),
        success=bool(success),
        volume_L=volume_L,
    )


def tau_leap_simulate(
    reactions: list[Reaction],
    species: list[str],
    rate_values: dict[str, float],
    initial_conditions: dict[str, float],
    t_span: tuple[float, float],
    volume_L: float,
    *,
    initial_as: str = "concentration",
    tau: float | None = None,
    epsilon: float = 0.03,
    n_record: int = 200,
    seed: int | None = None,
    chemostatted_values: dict[str, float] | None = None,
) -> StochasticResult:
    """
    Approximate Gillespie via τ-leap (Gillespie 2001; Cao, Gillespie, Petzold 2006).

    Instead of one reaction per step, fire all reactions in parallel over a
    timestep τ, sampling each firing count from Poisson(a_j · τ).  Much faster
    than direct SSA when populations are large and propensities change slowly,
    at the cost of accuracy: τ-leap is exact only in the limit of frequent
    reactions per leap (large molecule counts).

    Parameters
    ----------
    reactions, species, rate_values, initial_conditions, t_span, volume_L,
    initial_as, seed, chemostatted_values
        Same as :func:`gillespie_simulate`.
    tau : optional fixed leap size (seconds).  If None, an adaptive τ is
        chosen each step bounded by ``epsilon`` (Cao 2006 algorithm).
    epsilon : adaptive-τ tolerance.  Controls the maximum allowed fractional
        change in any propensity per leap; smaller = more accurate but slower.
        Ignored when ``tau`` is given.
    n_record : approximate number of evenly-spaced time points to record in
        the returned trajectory (the simulator's internal step is independent
        of this — recording is via linear interpolation of step boundaries).

    Returns
    -------
    StochasticResult
        Same shape as :func:`gillespie_simulate`.

    Notes
    -----
    To preserve non-negativity in stiff systems, this implementation falls
    back to a single direct-SSA step when the chosen τ would drive any
    species count below zero (a simple "leap rejection" — not the binomial
    Tian/Burrage refinement, but adequate for most kinetic-design contexts).
    """
    from .stoichiometry import fold_chemostatted_into_rates
    from ._kernels import build_kernel_arrays, normalize_seed, tau_leap_kernel

    chem_vals = chemostatted_values or {}
    eff_rates = (
        fold_chemostatted_into_rates(reactions, rate_values, chem_vals)
        if chem_vals else rate_values
    )
    chem = set(chem_vals.keys())

    sp_idx = {s: i for i, s in enumerate(species)}
    n_sp = len(species)
    NA_V = AVOGADRO * volume_L

    n0 = np.zeros(n_sp, dtype=np.int64)
    for s, val in initial_conditions.items():
        if s not in sp_idx:
            continue
        if initial_as == "count":
            n0[sp_idx[s]] = int(round(val))
        else:
            n0[sp_idx[s]] = int(round(val * NA_V))

    c_eff, react_sp, react_co, change_matrix = build_kernel_arrays(
        reactions, species, eff_rates, chem, NA_V
    )

    t0 = float(t_span[0])
    tf = float(t_span[1])
    record_times = np.linspace(t0, tf, max(2, n_record))

    # Adaptive (fixed_tau<=0) or fixed-step tau-leap, JIT-compiled when numba is
    # present, otherwise an equivalent pure-Python loop (see mantis._kernels).
    record_states, success = tau_leap_kernel(
        n0, c_eff, react_sp, react_co, change_matrix, t0, tf,
        float(tau) if tau is not None else -1.0, float(epsilon),
        record_times, normalize_seed(seed),
    )

    counts_dict = {sp: record_states[:, i] for i, sp in enumerate(species)}
    conc_dict = {sp: record_states[:, i] / NA_V for i, sp in enumerate(species)}

    return StochasticResult(
        times=record_times,
        counts=counts_dict,
        concentrations=conc_dict,
        n_events=int(len(record_times)),
        success=success,
        volume_L=volume_L,
    )


def build_ode_function(
    reactions: list[Reaction],
    species: list[str],
    rate_values: dict[str, float],
    chemostatted: set[str] | None = None,
) -> Any:
    """
    Return a callable f(t, y) → dydt using pure numpy.
    Pre-computes N matrix and ordered rate array for speed.

    Chemostatted species are excluded from reactants_info (their contribution
    has already been folded into rate_values via fold_chemostatted_into_rates).
    """
    chem = chemostatted or set()
    sp_idx = {s: i for i, s in enumerate(species)}
    N = build_stoichiometry_matrix(reactions, species)
    rates_arr = np.array([rate_values.get(r.rate_key, 0.0) for r in reactions])

    # Pre-compute reactant info: list of (species_index, coeff) per reaction
    # Skip chemostatted species — their concentrations are already in rates_arr
    reactants_info = [
        [(sp_idx[name], coeff) for name, coeff in rxn.reactants if name not in chem]
        for rxn in reactions
    ]

    def f(t, y):
        y = np.maximum(y, 0.0)
        fluxes = np.empty(len(reactions))
        for j, (rate, ri) in enumerate(zip(rates_arr, reactants_info)):
            flux = rate
            for idx, coeff in ri:
                flux *= y[idx] ** coeff
            fluxes[j] = flux
        return N @ fluxes

    return f


def _conservation_law_vectors(
    reactions: list[Reaction],
    species: list[str],
) -> list[np.ndarray]:
    """
    Return list of non-negative integer numpy vectors spanning the left null space of N.
    Uses the same sign-normalized basis as conservation_laws_sympy() so that
    each vector has all non-negative entries (physical moiety totals).
    """
    from .stoichiometry import conservation_laws_sympy
    import sympy
    N = build_stoichiometry_matrix(reactions, species)
    laws = conservation_laws_sympy(N, species)
    sp_sym_map = {s: sympy.Symbol(s) for s in species}
    result = []
    for law in laws:
        v = np.array([
            float(law.coeff(sp_sym_map[s])) for s in species
        ], dtype=float)
        result.append(v)
    return result


def _build_augmented_system(
    reactions: list[Reaction],
    species: list[str],
    rate_values: dict[str, float],
    cl_vectors: list[np.ndarray],
    cl_totals: list[float],
    chemostatted: set[str] | None = None,
) -> Any:
    """
    Build f(y) = 0 system combining independent ODE equations + conservation constraints.
    Selects rank(N) linearly independent rows of N using QR pivoting.
    """
    from scipy.linalg import qr
    chem = chemostatted or set()
    sp_idx = {s: i for i, s in enumerate(species)}
    N = build_stoichiometry_matrix(reactions, species)
    rates_arr = np.array([rate_values.get(r.rate_key, 0.0) for r in reactions])
    reactants_info = [
        [(sp_idx[name], coeff) for name, coeff in rxn.reactants if name not in chem]
        for rxn in reactions
    ]

    rank = matrix_rank(N)
    if rank > 0:
        _, _, perm = qr(N.T, pivoting=True)
        independent_rows = perm[:rank]
        N_ind = N[independent_rows, :]
    else:
        N_ind = np.zeros((0, len(reactions)))

    def system(y):
        y_pos = np.maximum(y, 0.0)
        fluxes = np.empty(len(reactions))
        for j, (rate, ri) in enumerate(zip(rates_arr, reactants_info)):
            flux = rate
            for idx, coeff in ri:
                flux *= y_pos[idx] ** coeff
            fluxes[j] = flux
        ode_eqs = N_ind @ fluxes
        cl_eqs = np.array([
            v @ y_pos - tot
            for v, tot in zip(cl_vectors, cl_totals)
        ])
        return np.concatenate([ode_eqs, cl_eqs])

    return system


def _make_feasible_initial(
    cl_vectors: list[np.ndarray],
    cl_totals: list[float],
    n_sp: int,
    rng: np.random.Generator,
) -> np.ndarray | None:
    """
    Generate a random initial point that approximately satisfies conservation constraints.
    Uses a simple proportional allocation strategy.
    """
    y0 = rng.uniform(1e-12, 1e-7, size=n_sp)
    # Scale each conservation law to match its total
    for _ in range(20):
        for v, tot in zip(cl_vectors, cl_totals):
            mass = v @ y0
            if mass > 1e-30:
                mask = v > 0.5  # species in this law
                y0[mask] *= tot / mass
        y0 = np.maximum(y0, 0.0)
    return y0


def _full_jacobian_fn(
    reactions: list[Reaction],
    species: list[str],
    rate_values: dict[str, float],
    chemostatted: set[str] | None = None,
):
    """Return a function y → J(y) for the full ODE Jacobian dF/dy."""
    from .stoichiometry import build_stoichiometry_matrix
    chem = chemostatted or set()
    sp_idx = {s: i for i, s in enumerate(species)}
    N_full = build_stoichiometry_matrix(reactions, species)
    n_sp = len(species)
    rates_arr = np.array([rate_values.get(r.rate_key, 0.0) for r in reactions])
    reactants_info = [
        [(sp_idx[name], coeff) for name, coeff in rxn.reactants if name not in chem]
        for rxn in reactions
    ]

    def full_jac(y):
        y_pos = np.maximum(y, 0.0)
        J = np.zeros((n_sp, n_sp))
        for j, (rate, ri) in enumerate(zip(rates_arr, reactants_info)):
            flux = rate
            for idx, coeff in ri:
                flux *= y_pos[idx] ** coeff
            for k_idx, k_coeff in ri:
                if y_pos[k_idx] > 1e-30:
                    d_flux = flux * k_coeff / y_pos[k_idx]
                else:
                    d_flux = rate * k_coeff * (y_pos[k_idx] ** max(k_coeff - 1, 0))
                    for other_idx, other_coeff in ri:
                        if other_idx != k_idx:
                            d_flux *= y_pos[other_idx] ** other_coeff
                J[:, k_idx] += N_full[:, j] * d_flux
        return J

    return full_jac


def simulate_ode(
    reactions: list[Reaction],
    species: list[str],
    rate_values: dict[str, float],
    initial_conditions: dict[str, float],
    t_span: tuple[float, float],
    t_eval: np.ndarray | None = None,
    chemostatted_values: dict[str, float] | None = None,
    rtol: float = 1e-8,
    atol: float = 1e-12,
) -> SimulationResult:
    """
    Integrate the ODE system forward in time and return the full trajectory.

    Parameters
    ----------
    t_span : (t0, tf)
        Start and end times in seconds.
    t_eval : array-like, optional
        Times at which to store the solution.  Defaults to 200 log-spaced
        points across t_span (or linear if t_span[0] == 0 and t_span[1] <= 0).
    chemostatted_values : dict, optional
        Fixed concentrations folded into rate constants (same semantics as
        ``find_steady_states``).

    Returns
    -------
    SimulationResult
        ``.times`` (1-D array), ``.concentrations`` (dict species → 1-D array),
        ``.success`` (bool).  Use ``.final()`` to get the last time-point dict
        or ``.at(t)`` for a specific time.
    """
    from scipy.integrate import solve_ivp
    from .stoichiometry import fold_chemostatted_into_rates

    chem_vals = chemostatted_values or {}
    eff_rates = fold_chemostatted_into_rates(reactions, rate_values, chem_vals) if chem_vals else rate_values
    chem_keys = set(chem_vals.keys())

    f_ode = build_ode_function(reactions, species, eff_rates, chem_keys)
    y0 = np.array([initial_conditions.get(s, 0.0) for s in species])

    t0, tf = float(t_span[0]), float(t_span[1])
    if t_eval is None:
        n_pts = 200
        if t0 <= 0 and tf > 0:
            t_eval = np.logspace(np.log10(max(tf * 1e-6, 1e-6)), np.log10(tf), n_pts - 1)
            t_eval[-1] = tf
            t_eval = np.concatenate([[t0], t_eval])
        else:
            t_eval = np.linspace(t0, tf, n_pts)

    try:
        sol = solve_ivp(
            f_ode,
            [t0, tf],
            y0,
            method="Radau",
            t_eval=t_eval,
            rtol=rtol,
            atol=atol,
            dense_output=False,
        )
        times = sol.t
        conc = {sp: np.maximum(sol.y[i], 0.0) for i, sp in enumerate(species)}
        if not sol.success:
            print("solve_ivp failed:", getattr(sol, 'message', 'No message'))
        return SimulationResult(times=times, concentrations=conc, success=sol.success)
    except Exception as e:
        import traceback
        traceback.print_exc()
        times = np.array([t0, tf])
        conc = {sp: np.full(2, y0[i]) for i, sp in enumerate(species)}
        return SimulationResult(times=times, concentrations=conc, success=False)


def _integrate_to_ss(
    f_ode: Any,
    y0: np.ndarray,
    t_end: float = 1e4,
    rtol: float = 1e-8,
    atol: float = 1e-12,
) -> tuple[np.ndarray | None, float]:
    """
    Integrate ODE from y0 to t_end and return the final state.
    Returns (y_final, residual_norm) or (None, inf) on failure.

    Note: t_end is intentionally short (1e4 default) so that limit-cycle systems
    (unstable fixed points / oscillators) do not cause the solver to hang.
    The residual check downstream will reject non-converged trajectories, and the
    algebraic least_squares fallback handles those cases.

    Only call this for systems with conservation laws (closed/semi-closed networks).
    For fully open systems with no CLs, skip integration and use least_squares directly.
    """
    from scipy.integrate import solve_ivp
    try:
        sol = solve_ivp(
            f_ode,
            [0, t_end],
            y0,
            method="Radau",
            rtol=rtol,
            atol=atol,
            dense_output=False,
        )
        if not sol.success:
            return None, float("inf")
        y_final = np.maximum(sol.y[:, -1], 0.0)
        residual = float(np.linalg.norm(f_ode(0, y_final)))
        return y_final, residual
    except Exception:
        return None, float("inf")


def find_steady_states(
    reactions: list[Reaction],
    species: list[str],
    rate_values: dict[str, float],
    initial_conditions: dict[str, float],
    n_attempts: int = 50,
    tol: float = 1e-8,
    seed: int | None = None,
    chemostatted_values: dict[str, float] | None = None,
    t_end: float = 1e4,
) -> list[SteadyState]:
    """
    Steady-state finder using three strategies:
    1. Direct least_squares from user IC (finds both stable and unstable fixed points).
    2. ODE integration from user IC to t_end (conserves CLs; reaches stable attractors).
    3. Multi-start with scale-aware random ICs → ODE then least_squares fallback.

    Chemostatted species are folded into effective rate constants so the ODE
    system only tracks dynamic species.

    Parameters
    ----------
    t_end : float
        Integration horizon for ODE-based strategies (seconds).  The default
        (1e4 s) is intentionally short so that oscillatory / limit-cycle
        systems do not hang; the algebraic least_squares fallback handles those
        cases.  For systems with very slow reactions (e.g. leakage pathways
        with τ >> 1e4 s) increase this value so that ODE integration reaches
        the true attractor.  A safe upper bound is 10–100× the slowest
        relevant timescale (τ = 1 / (k_slow × [reactant])).

    Returns a list of SteadyState objects, sorted by true ODE residual (lowest first),
    filtered to remove duplicate states and high-residual algebraic artifacts.
    """
    from .stoichiometry import fold_chemostatted_into_rates
    chem_vals = chemostatted_values or {}
    chem_keys = set(chem_vals.keys())
    eff_rates = fold_chemostatted_into_rates(reactions, rate_values, chem_vals) if chem_vals else rate_values

    n_sp = len(species)
    rng = np.random.default_rng(seed)

    cl_vecs = _conservation_law_vectors(reactions, species)
    y_ref = np.array([initial_conditions.get(s, 0.0) for s in species])
    cl_totals = [float(v @ y_ref) for v in cl_vecs]

    f_ode = build_ode_function(reactions, species, eff_rates, chem_keys)
    full_jac = _full_jacobian_fn(reactions, species, eff_rates, chem_keys)
    system_fn = _build_augmented_system(reactions, species, eff_rates, cl_vecs, cl_totals, chem_keys)

    collected: list[SteadyState] = []

    def _try_add(y_sol: np.ndarray, _residual_hint: float) -> None:
        y_sol = np.maximum(y_sol, 0.0)
        # Recompute the true ODE residual |dy/dt|.  We store it for post-filtering
        # (relative to the best-found solution) rather than applying a fixed absolute
        # threshold, which fails on stiff systems where all flux magnitudes are tiny.
        abs_res = float(np.linalg.norm(f_ode(0, y_sol)))
        if np.any(y_sol < -tol):
            return
        for existing in collected:
            y_ex = np.array([existing.concentrations[s] for s in species])
            if np.linalg.norm(y_sol - y_ex) / (np.linalg.norm(y_ex) + 1e-30) < 0.10:
                return
        J_num = full_jac(y_sol)
        eigs = eigvals(J_num)
        is_stable, is_osc = classify_steady_state(eigs)
        ss = SteadyState(
            concentrations=dict(zip(species, y_sol.tolist())),
            eigenvalues=eigs,
            is_stable=is_stable,
            is_oscillatory=is_osc,
            residual=abs_res,
        )
        collected.append(ss)

    def _try_least_squares(y0: np.ndarray) -> None:
        """Try algebraic least_squares from y0; can find unstable fixed points."""
        try:
            result = least_squares(
                system_fn,
                y0,
                bounds=(0.0, np.inf),
                method="trf",
                ftol=tol,
                xtol=tol,
                gtol=tol,
                max_nfev=5000,
            )
            _try_add(result.x, float(np.linalg.norm(result.fun)))
        except Exception:
            pass

    has_cl = len(cl_vecs) > 0

    if has_cl:
        # ── Closed / semi-closed system (has conservation laws) ──────────────
        # Strategy 1a: ODE integration from user IC — respects CLs, reaches the
        # physically relevant attractor (unique within each conservation class).
        y_final, res = _integrate_to_ss(f_ode, y_ref.copy(), t_end=t_end)
        if y_final is not None:
            _try_add(y_final, res)

        # Strategy 1b: algebraic solve from user IC — finds any fixed point
        # including those that are unstable (integration diverges from them).
        _try_least_squares(y_ref.copy())

        # Strategy 2: multi-start on the constraint manifold
        scale = max(float(np.max(np.abs(y_ref))), 1.0)
        for _ in range(n_attempts - 2):
            y0 = _make_feasible_initial(cl_vecs, cl_totals, n_sp, rng)
            if y0 is None:
                continue
            y_ivp, res_ivp = _integrate_to_ss(f_ode, y0, t_end=t_end)
            if y_ivp is not None and res_ivp < 1e-4:
                _try_add(y_ivp, res_ivp)
                continue
            _try_least_squares(y0)

    else:
        # ── Open / chemostatted system (no conservation laws) ─────────────────
        # ODE integration is useless here: limit-cycle attractors make it hang,
        # and there is no conservation manifold to project onto.  Use algebraic
        # least_squares exclusively.

        # Strategy 1: algebraic solve from user IC (finds both stable & unstable FPs)
        _try_least_squares(y_ref.copy())

        # Strategy 2: scale-aware multi-start algebraic solve
        scale = max(float(np.max(np.abs(y_ref))), 1.0)
        for _ in range(n_attempts - 1):
            y0 = np.abs(rng.normal(loc=scale, scale=scale * 0.5, size=n_sp))
            y0 = np.maximum(y0, 1e-12)
            _try_least_squares(y0)

    collected.sort(key=lambda s: s.residual)
    # Post-filter: drop candidates whose ODE residual is more than 1e4× worse than
    # the best (lowest-residual) solution.  This correctly handles stiff systems
    # where all fluxes are tiny in absolute terms — a relative comparison is the
    # only meaningful discriminant between "converged" and "stuck-at-IC" solutions.
    if collected:
        best_res = collected[0].residual
        collected = [s for s in collected if s.residual <= best_res * 1e4 + 1e-30]
    return collected


def classify_steady_state(
    eigenvalues: np.ndarray,
    tol_zero: float = 1e-4,
) -> tuple[bool, bool]:
    """
    Returns (is_stable, is_oscillatory).

    ``is_stable``     — True iff all significant eigenvalues have Re < 0.
    ``is_oscillatory`` — True iff any significant eigenvalue has a non-trivial
                         imaginary part (|Im| > 1e-10 * |λ|).  This covers both
                         stable spirals (Re<0) *and* unstable foci/spirals (Re>0),
                         i.e. any fixed point near a Hopf bifurcation.

    Filters near-zero eigenvalues (from conservation law dimensions) before classifying.
    """
    if len(eigenvalues) == 0:
        return True, False
    max_abs = np.max(np.abs(eigenvalues))
    if max_abs < 1e-30:
        return True, False
    # Filter eigenvalues that are effectively zero relative to the largest
    significant = eigenvalues[np.abs(eigenvalues) > tol_zero * max_abs]
    if len(significant) == 0:
        return True, False
    is_stable = bool(np.all(np.real(significant) < 0))
    # Oscillatory = complex eigenvalues regardless of stability
    # (captures both stable spirals and unstable foci / Hopf-born limit cycles)
    is_oscillatory = bool(
        np.any(np.abs(np.imag(significant)) > 1e-10 * np.abs(significant))
    )
    return is_stable, is_oscillatory


def compute_eigenvalues(J_numeric: np.ndarray) -> np.ndarray:
    return eigvals(J_numeric)


def scan_bifurcation(
    reactions: list[Reaction],
    species: list[str],
    base_rate_values: dict[str, float],
    parameter: str,
    param_range: tuple[float, float],
    n_points: int,
    initial_conditions: dict[str, float],
    n_attempts: int = 20,
    chemostatted_values: dict[str, float] | None = None,
    t_end: float = 1e4,
) -> BifurcationResult:
    """
    Vary one rate constant over a log-spaced range and collect steady states.

    Parameters
    ----------
    t_end : float
        ODE integration horizon passed to ``find_steady_states`` at each point.
        See that function's ``t_end`` documentation for guidance on choosing a
        value appropriate for the slowest reaction timescale in your network.
    """
    from .parsing import normalize_rate_key
    try:
        norm_param = normalize_rate_key(parameter)
    except ValueError:
        norm_param = parameter

    param_values = np.logspace(
        np.log10(param_range[0]), np.log10(param_range[1]), n_points
    )
    all_ss: list[list[SteadyState]] = []
    for pval in param_values:
        rates = dict(base_rate_values)
        rates[norm_param] = pval
        ss_list = find_steady_states(
            reactions, species, rates, initial_conditions,
            n_attempts=n_attempts, seed=42,
            chemostatted_values=chemostatted_values,
            t_end=t_end,
        )
        all_ss.append(ss_list)

    return BifurcationResult(
        parameter_name=norm_param,
        parameter_values=param_values,
        steady_states=all_ss,
    )
