"""CRNetwork: unified interface for CRNT structural analysis, symbolic ODEs, ODE simulation, and steady-state finding."""
from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np

from .parsing import Reaction, Complex, parse_reactions, normalize_rate_key
from .stoichiometry import (
    build_species_list,
    build_stoichiometry_matrix,
    matrix_rank,
    conservation_laws_sympy,
)
from .crnt import (
    build_complex_graph,
    get_linkage_classes,
    check_weak_reversibility,
    compute_deficiency,
    crnt_analysis,
    CRNTResult,
)

if TYPE_CHECKING:
    import sympy
    import networkx as nx
    import matplotlib.axes


class CRNetwork:
    """
    Unified interface for Chemical Reaction Network Theory (CRNT) analysis.

    Construct with `CRNetwork.from_string(reaction_strings, rates=...)`.
    Structural properties (deficiency, conservation laws, weak reversibility) are
    computed lazily and cached; they do not require rate constants.  Symbolic and
    numerical methods require rate constants to be supplied.

    Parameters
    ----------
    reactions : list[Reaction]
        Directed Reaction objects.  Each reversible string `A <-> B` produces two.
    rates : dict[str, float], optional
        Rate constants keyed by reaction strings, e.g. ``{"A -> B": 1.0}``.
        Keys are automatically normalized to canonical form (species sorted
        alphabetically); use ``rate_keys()`` to inspect the expected form.
    """

    def __init__(
        self,
        reactions: list[Reaction],
        rates: dict[str, float] | None = None,
        chemostatted: dict[str, float] | None = None,
    ) -> None:
        self._reactions = reactions
        self._rates: dict[str, float] = {}
        if rates:
            # Normalize user-supplied keys to canonical form
            for key, val in rates.items():
                try:
                    norm = normalize_rate_key(key)
                except ValueError:
                    norm = key
                self._rates[norm] = val
        self._chemostatted: dict[str, float] = dict(chemostatted) if chemostatted else {}

    @classmethod
    def from_string(
        cls,
        reaction_strings: list[str],
        rates: dict[str, float] | None = None,
        chemostatted: dict[str, float] | None = None,
    ) -> CRNetwork:
        """
        Construct from human-readable reaction strings.

        Parameters
        ----------
        reaction_strings : list[str]
            Each string is one reaction, e.g. ``"A + 2B <-> C"`` or ``"ES -> E + P"``.
        rates : dict[str, float], optional
            Maps reaction strings to rate constants.  Order of species within each
            side does not matter; keys are normalized before lookup.
        chemostatted : dict[str, float], optional
            Species held at fixed concentrations by an external reservoir.
            They are excluded from the ODE system and stoichiometry matrix rows,
            but their concentrations are folded into flux expressions.
        """
        reactions = parse_reactions(reaction_strings)
        return cls(reactions, rates, chemostatted)

    # ── Structural properties (no rates needed) ──────────────────────────────

    @cached_property
    def species(self) -> list[str]:
        return build_species_list(self._reactions, set(self._chemostatted.keys()))

    @cached_property
    def complexes(self) -> list[Complex]:
        return list(self._crnt_graph.nodes())

    @cached_property
    def stoichiometry_matrix(self) -> np.ndarray:
        return build_stoichiometry_matrix(self._reactions, self.species)

    @cached_property
    def n_species(self) -> int:
        return len(self.species)

    @cached_property
    def n_complexes(self) -> int:
        return len(self.complexes)

    @cached_property
    def n_reactions(self) -> int:
        return len(self._reactions)

    @cached_property
    def _crnt_graph(self) -> nx.DiGraph:
        return build_complex_graph(self._reactions, set(self._chemostatted.keys()))

    @cached_property
    def n_linkage_classes(self) -> int:
        return len(get_linkage_classes(self._crnt_graph))

    @cached_property
    def deficiency(self) -> int:
        return self._crnt_result.deficiency

    @cached_property
    def is_weakly_reversible(self) -> bool:
        return self._crnt_result.is_weakly_reversible

    @cached_property
    def conservation_laws(self) -> list[sympy.Expr]:
        return conservation_laws_sympy(self.stoichiometry_matrix, self.species)

    @cached_property
    def _crnt_result(self) -> CRNTResult:
        return crnt_analysis(self._reactions, self.stoichiometry_matrix, self.species, set(self._chemostatted.keys()))

    @property
    def chemostatted(self) -> dict[str, float]:
        """Chemostatted species and their fixed concentrations."""
        return dict(self._chemostatted)

    def rate_keys(self) -> list[str]:
        """Return canonical rate key strings for all reactions (for debugging)."""
        return [r.rate_key for r in self._reactions]

    def crnt_summary(self) -> str:
        """Return a human-readable CRNT analysis summary."""
        r = self._crnt_result
        lines = list(r.summary_lines)
        if self._chemostatted:
            chem_str = ", ".join(f"{s}={v}" for s, v in sorted(self._chemostatted.items()))
            lines.insert(0, f"Chemostatted species: {chem_str}")
            lines.insert(1, "")
        lines.append("")
        laws = self.conservation_laws
        lines.append(f"Conservation laws ({len(laws)}):")
        for i, law in enumerate(laws, 1):
            lines.append(f"  [{i}] {law} = const")
        return "\n".join(lines)

    # ── Symbolic ──────────────────────────────────────────────────────────────

    def odes(self, numeric_rates: bool = True) -> dict[str, sympy.Expr]:
        """
        Return mass-action ODEs as SymPy expressions keyed by species name.

        Parameters
        ----------
        numeric_rates : bool
            If True (default) and rates were supplied, substitute numerical values
            so that the returned expressions contain floats rather than ``k_i`` symbols.
        """
        from .symbolic import (
            make_species_symbols,
            make_rate_symbols,
            build_odes,
            substitute_rates,
        )
        sp_syms = make_species_symbols(self.species)
        rate_syms, key_to_sym = make_rate_symbols(self._reactions)
        odes_sym = build_odes(
            self._reactions, self.species, sp_syms, rate_syms, key_to_sym,
            chemostatted_values=self._chemostatted or None,
        )
        if numeric_rates and self._rates:
            odes_sym = {
                sp: substitute_rates(expr, key_to_sym, self._rates)
                for sp, expr in odes_sym.items()
            }
        return odes_sym

    def jacobian(self) -> sympy.Matrix:
        """Return the symbolic Jacobian matrix (n_species × n_species), cached after first call."""
        if "_jacobian_cache" not in self.__dict__:
            from .symbolic import make_species_symbols, build_jacobian
            sp_syms = make_species_symbols(self.species)
            J = build_jacobian(self.odes(numeric_rates=False), self.species, sp_syms)
            self.__dict__["_jacobian_cache"] = J
        return self.__dict__["_jacobian_cache"]

    def jacobian_at(self, ss: dict[str, float]) -> sympy.Matrix:
        from .symbolic import make_species_symbols, substitute_rates, substitute_steady_state
        J_sym = self.jacobian()
        sp_syms = make_species_symbols(self.species)
        from .symbolic import make_rate_symbols
        rate_syms, key_to_sym = make_rate_symbols(self._reactions)
        J_rates = substitute_rates(J_sym, key_to_sym, self._rates)
        return substitute_steady_state(J_rates, sp_syms, ss)

    # ── Numerical ─────────────────────────────────────────────────────────────

    def steady_states(
        self,
        initial_conditions: dict[str, float],
        n_attempts: int = 50,
        seed: int | None = None,
        t_end: float = 1e4,
    ) -> list:
        """
        Find steady states from the given initial conditions.

        Uses ODE integration (Radau) as the primary strategy, which exactly
        preserves conservation laws, then falls back to multi-start least_squares
        for additional steady states. Duplicate solutions (relative L2 distance
        < 10%) and non-physical solutions (any concentration < −tol) are discarded.
        Results are sorted by residual ascending and filtered by relative residual
        against the best solution found to remove spurious stuck-at-IC artifacts.

        Parameters
        ----------
        initial_conditions : dict[str, float]
            Initial concentrations; missing species default to 0.
        n_attempts : int
            Total number of solver attempts (including the integration from IC).
        seed : int, optional
            RNG seed for reproducibility.
        t_end : float
            ODE integration horizon in seconds (default 1e4).  Increase for
            systems with slow reactions whose timescale τ = 1/(k·[X]) >> 1e4 s.

        Returns
        -------
        list[SteadyState]
            Each entry has ``.concentrations``, ``.eigenvalues``, ``.is_stable``,
            ``.is_oscillatory``, and ``.residual``.
        """
        from .analysis import find_steady_states
        return find_steady_states(
            self._reactions,
            self.species,
            self._rates,
            initial_conditions,
            n_attempts=n_attempts,
            seed=seed,
            chemostatted_values=self._chemostatted or None,
            t_end=t_end,
        )

    def all_steady_states(
        self,
        initial_conditions: dict[str, float],
        backend: str = "auto",
        scale: float | None = None,
    ) -> list:
        """
        Exhaustively enumerate *all* positive steady states (with completeness).

        Unlike :meth:`steady_states` (ODE integration + multi-start least-squares,
        which can silently miss unstable fixed points), this solves the mass-action
        steady-state polynomial system for *every* complex root and returns the real
        non-negative ones — so it certifies the full steady-state set rather than
        sampling it.  The classic payoff: on the Schlögl model it returns all three
        equilibria, including the unstable middle branch that integration never
        reaches.

        Parameters
        ----------
        initial_conditions : dict[str, float]
            Used only to fix the conservation-law totals (the stoichiometric
            compatibility class) for closed networks.
        backend : {'auto', 'sympy', 'phcpy'}
            ``'sympy'`` is the pure-Python Gröbner-basis enumerator (default via
            ``'auto'``).  ``'phcpy'`` uses optional polynomial homotopy continuation
            (PHCpack) for larger / stiffer systems — ``pip install mantis-delta[homotopy]``.
        scale : float, optional
            Characteristic concentration used to non-dimensionalise the variables.
            Inferred from the conservation totals / initial conditions when omitted;
            override only if the solver struggles on a badly-scaled network.

        Returns
        -------
        list[SteadyState]
            One entry per distinct positive steady state, sorted by residual, each
            with ``.eigenvalues``, ``.is_stable``, ``.is_oscillatory``.
        """
        from .steady_states import all_steady_states
        return all_steady_states(
            self._reactions,
            self.species,
            self._rates,
            initial_conditions,
            chemostatted_values=self._chemostatted or None,
            backend=backend,
            scale=scale,
        )

    def bifurcation(
        self,
        parameter: str,
        range: tuple[float, float],
        n_points: int = 100,
        initial_conditions: dict[str, float] | None = None,
        plot: bool = False,
        t_end: float = 1e4,
    ) -> object:
        """
        Scan one rate constant over a log-spaced range and collect steady states.

        Parameters
        ----------
        parameter : str
            Rate key to vary, e.g. ``"miR21 + H1 -> miR21_H1"``.
        range : tuple[float, float]
            (min, max) values for the parameter (log scale).
        n_points : int
            Number of points in the scan.
        initial_conditions : dict[str, float], optional
            Initial concentrations for each scan point (defaults to all-zero).
        plot : bool
            If True, display a bifurcation diagram for the first species.
        t_end : float
            ODE integration horizon in seconds (default 1e4).  Increase for
            systems with slow reactions whose timescale τ = 1/(k·[X]) >> 1e4 s.

        Returns
        -------
        BifurcationResult
        """
        from .analysis import scan_bifurcation
        result = scan_bifurcation(
            self._reactions,
            self.species,
            self._rates,
            parameter=parameter,
            param_range=range,
            n_points=n_points,
            initial_conditions=initial_conditions or {},
            n_attempts=20,
            chemostatted_values=self._chemostatted or None,
            t_end=t_end,
        )
        if plot:
            from .plot import plot_bifurcation
            plot_bifurcation(result, species=self.species[0])
        return result

    def continuation(
        self,
        parameter: str,
        range: tuple[float, float],
        initial_conditions: dict[str, float],
        ds: float = 0.05,
        max_steps: int = 2000,
        initial_state: dict[str, float] | None = None,
    ) -> object:
        """
        Trace a steady-state branch by **pseudo-arclength continuation**.

        Unlike :meth:`bifurcation` (an independent log-scan that loses the branch
        at folds and never reports them), this parameterises the branch by
        arclength along the solution curve, so saddle-node folds are traversed
        smoothly and detected.  The classic payoff: on the Schlögl model it traces
        the full S-curve and pins both folds that bound the bistable region.

        Parameters
        ----------
        parameter : str
            Rate key to continue, e.g. ``"A + 2X -> 3X"``.
        range : tuple[float, float]
            ``(λ_min, λ_max)`` box the branch is traced within; to capture an
            S-curve the box must bracket both folds.
        initial_conditions : dict[str, float]
            Fixes the conservation totals (the compatibility class) and seeds the
            branch (a steady state is solved at ``λ_min`` unless ``initial_state``
            is given).
        ds : float
            Arclength step (adapted automatically on Newton failure / success).
        max_steps : int
            Cap on continuation steps (guards against closed isolas).
        initial_state : dict[str, float], optional
            Explicit starting steady state — use to select which branch to follow.

        Returns
        -------
        ContinuationResult
            With ``.parameter_values`` (arclength order — non-monotonic across
            folds), ``.branch``, ``.stable``, ``.bifurcations`` (``.folds()`` /
            ``.hopfs()``), and ``.species_branch(species)``.  ``str()`` summarises.
        """
        from .continuation import pseudo_arclength_continuation
        return pseudo_arclength_continuation(
            self._reactions,
            self.species,
            self._rates,
            parameter=parameter,
            param_range=range,
            initial_conditions=initial_conditions,
            chemostatted_values=self._chemostatted or None,
            ds=ds,
            max_steps=max_steps,
            initial_state=initial_state,
        )

    def simulate(
        self,
        initial_conditions: dict[str, float],
        t_span: tuple[float, float],
        t_eval=None,
        rtol: float = 1e-8,
        atol: float = 1e-12,
    ):
        """
        Integrate the ODE system forward in time and return the full trajectory.

        Parameters
        ----------
        initial_conditions : dict[str, float]
            Initial concentrations; missing species default to 0.
        t_span : (t0, tf)
            Start and end times in seconds.
        t_eval : array-like, optional
            Times at which to store the solution.  Defaults to 200 log-spaced
            points across t_span.
        rtol, atol : float
            Solver tolerances (passed to scipy Radau).

        Returns
        -------
        SimulationResult
            ``.times``, ``.concentrations`` (dict species → 1-D array),
            ``.success``.  Use ``.final()`` for the last time-point dict
            or ``.at(t)`` for a specific time.
        """
        from .analysis import simulate_ode
        return simulate_ode(
            self._reactions,
            self.species,
            self._rates,
            initial_conditions,
            t_span=t_span,
            t_eval=t_eval,
            chemostatted_values=self._chemostatted or None,
            rtol=rtol,
            atol=atol,
        )

    def stochastic_simulate(
        self,
        initial_conditions: dict[str, float],
        t_span: tuple[float, float],
        volume_L: float,
        *,
        initial_as: str = "concentration",
        max_events: int = 1_000_000,
        seed: int | None = None,
    ):
        """
        Single-trajectory Gillespie SSA realization.

        Use when molecule counts are low (∼ ≤ 10³) and the deterministic ODE
        gives the wrong answer — e.g., a CHA cascade at single-cell
        concentrations or stochastic switching in a bistable circuit.

        Parameters
        ----------
        initial_conditions : species → initial count or concentration.
        t_span             : (t0, tf) in seconds.
        volume_L           : reaction volume in liters (e.g. 1e-4 = 100 µL).
        initial_as         : 'concentration' (default) or 'count'.
        max_events         : safety cap on reaction firings.
        seed               : RNG seed for reproducibility.

        Returns
        -------
        StochasticResult with ``.times``, ``.counts``, ``.concentrations``,
        ``.at(t)``, ``.final()``.
        """
        from .analysis import gillespie_simulate
        return gillespie_simulate(
            self._reactions,
            self.species,
            self._rates,
            initial_conditions,
            t_span=t_span,
            volume_L=volume_L,
            initial_as=initial_as,
            max_events=max_events,
            seed=seed,
            chemostatted_values=self._chemostatted or None,
        )

    def stationary_distribution(
        self,
        initial_conditions: dict[str, float],
        volume_L: float,
        *,
        initial_as: str = "concentration",
        max_states: int = 200_000,
    ):
        """
        Exact stochastic stationary distribution (Anderson–Craciun–Kurtz 2010).

        For a complex-balanced network the chemical master equation has a closed-form
        product-of-Poissons stationary distribution with means ``λ_i = N_A·V·c_i*``
        (``c*`` the complex-balanced equilibrium) — written down *from structure
        alone, with no simulation*. With conservation laws it is that product
        conditioned on the conserved totals (e.g. a Binomial for a single moiety).

        Parameters
        ----------
        initial_conditions : dict[str, float]
            Fixes the compatibility class; interpreted per ``initial_as``.
        volume_L : float
            Reaction volume in litres (sets the absolute molecule counts).
        initial_as : 'concentration' or 'count'
        max_states : int
            Cap on enumerated states for the conditioned (closed-network) case.

        Returns
        -------
        StationaryDistribution
            With ``.poisson_means()``, ``.probability(state)``, ``.expected_counts()``,
            ``.marginal(species)`` and ``.sample(...)``.

        Raises
        ------
        ValueError
            If the network is not complex-balanced (the closed form does not apply).
        """
        from .stochastic_stationary import stationary_distribution
        return stationary_distribution(
            self._reactions, self.species, self._rates,
            self.deficiency, self.is_weakly_reversible,
            initial_conditions, volume_L,
            initial_as=initial_as,
            chemostatted_values=self._chemostatted or None,
            max_states=max_states,
        )

    def tau_leap_simulate(
        self,
        initial_conditions: dict[str, float],
        t_span: tuple[float, float],
        volume_L: float,
        *,
        initial_as: str = "concentration",
        tau: float | None = None,
        epsilon: float = 0.03,
        n_record: int = 200,
        seed: int | None = None,
    ):
        """
        τ-leap stochastic simulation (approximate Gillespie).

        Same interface as :meth:`stochastic_simulate` but fires reactions in
        Poisson-distributed bursts over each leap rather than one at a time.
        Roughly N× faster than direct SSA when populations are large (N is
        the mean number of firings per leap), at the cost of asymptotic
        accuracy as populations shrink.

        Parameters
        ----------
        tau : optional fixed leap size (s); if None, adaptive τ per Cao 2006.
        epsilon : adaptive-τ tolerance (max fractional propensity change per
                  leap).  Smaller = more accurate, slower.
        n_record : number of evenly-spaced time points to record.

        See :func:`mantis.analysis.tau_leap_simulate` for full details.
        """
        from .analysis import tau_leap_simulate
        return tau_leap_simulate(
            self._reactions,
            self.species,
            self._rates,
            initial_conditions,
            t_span=t_span,
            volume_L=volume_L,
            initial_as=initial_as,
            tau=tau,
            epsilon=epsilon,
            n_record=n_record,
            seed=seed,
            chemostatted_values=self._chemostatted or None,
        )

    def fsp(
        self,
        initial_conditions: dict[str, float],
        volume_L: float,
        t: float | None = None,
        *,
        initial_as: str = "concentration",
        max_states: int = 100_000,
    ):
        """
        Solve the chemical master equation by Finite State Projection (Munsky–Khammash 2006).

        Truncates the count-state space to the reachable set, builds the CME generator,
        and solves it directly — accurate where the SSA is noisy (low molecule numbers),
        with a rigorous truncation-error bound.

        Parameters
        ----------
        initial_conditions : dict[str, float]
            Starting state; interpreted per ``initial_as``.
        volume_L : float
            Reaction volume in litres.
        t : float, optional
            If given, returns the transient distribution ``p(t)``; if ``None``, returns
            the stationary distribution (valid for closed, fully enumerated networks).
        initial_as : 'concentration' or 'count'
        max_states : int
            Cap on the number of enumerated states.

        Returns
        -------
        FSPResult
            With ``.probabilities``, ``.states``, ``.marginal(species)``,
            ``.expected_counts()``, and ``.truncation_error``.
        """
        from .fsp import fsp_solve
        return fsp_solve(
            self._reactions, self.species, self._rates,
            initial_conditions, volume_L, t=t, initial_as=initial_as,
            chemostatted_values=self._chemostatted or None,
            max_states=max_states,
        )

    def draw(self, ax=None, layout: str = "spring") -> matplotlib.axes.Axes:
        from .plot import draw_reaction_graph
        return draw_reaction_graph(self._crnt_graph, ax=ax, layout=layout)

    # ── Modern structural analysis ──────────────────────────────────────────────

    def detect_acr(
        self, initial_conditions: dict[str, float] | None = None
    ):
        """
        Detect Absolute Concentration Robustness (Shinar–Feinberg, *Science* 2010).

        A deficiency-one network with two non-terminal complexes differing in a single
        species *S* has ACR in *S*: its steady-state concentration is independent of the
        initial totals — the robustness property at the heart of biosensor design.

        Parameters
        ----------
        initial_conditions : dict[str, float], optional
            When supplied (with rate constants set), the robust steady-state value of
            each ACR species is computed at a positive steady state.

        Returns
        -------
        ACRResult
            With ``.has_acr``, ``.species``, ``.acr_values``, and a readable summary
            via ``str()``.
        """
        from .acr import detect_acr
        return detect_acr(
            self._reactions, self.species, self._crnt_graph, self.deficiency,
            rate_values=self._rates or None,
            initial_conditions=initial_conditions,
            chemostatted_values=self._chemostatted or None,
        )

    def is_injective(self):
        """
        Test mass-action injectivity (Craciun–Feinberg, *SIAM J. Appl. Math.* 2005/06).

        An injective network has **at most one** positive steady state in each
        stoichiometric compatibility class — ruling out multistationarity. Unlike the
        Deficiency One Theorem this works for higher-deficiency networks too. The test
        is the sign-definiteness of the steady-state Jacobian determinant over all
        positive concentrations and rate constants (a sufficient condition).

        Returns
        -------
        InjectivityResult
            With ``.injective`` (True ⇒ certified injective; False ⇒ inconclusive),
            ``.determinant``, and a readable summary via ``str()``.
        """
        from .injectivity import test_injectivity
        return test_injectivity(
            self._reactions, self.species,
            chemostatted_values=self._chemostatted or None,
        )

    def multistationarity_region(self):
        """
        Compute the multistationarity parameter region (Conradi–Feliu–Mincheva–Wiuf,
        *PLoS Comput. Biol.* 2017).

        Where :meth:`is_injective` answers *whether* multistationarity is possible,
        this maps *where in parameter space* it occurs.  It forms the critical
        function ``φ = (−1)^s·det J`` (``s = rank(N)``) and reads it as a polynomial
        in the steady-state concentrations: if every coefficient is positive for all
        positive rates the network is monostationary for all parameters; otherwise
        each coefficient that can go negative is a face of the multistationarity
        region.  When rate constants are set, also reports whether *those* rates lie
        inside the region.

        Returns
        -------
        MultistationarityResult
            With ``.monostationary``, ``.multistationary_possible``,
            ``.critical_function``, ``.region_conditions``,
            ``.multistationary_at_rates``, and a readable summary via ``str()``.
        """
        from .multistationarity import multistationarity_region
        return multistationarity_region(
            self._reactions, self.species,
            rate_values=self._rates or None,
            chemostatted_values=self._chemostatted or None,
        )

    def is_multistationary(
        self, initial_conditions: dict[str, float], backend: str = "auto"
    ) -> bool:
        """
        Exactly decide multistationarity in one compatibility class by enumeration.

        Whereas :meth:`multistationarity_region` gives the rate-free *structural*
        verdict (and the symbolic region), this is the concrete numeric answer for
        the supplied rate constants and totals: it enumerates *all* positive steady
        states (:meth:`all_steady_states`) in the class fixed by ``initial_conditions``
        and returns ``True`` iff there are at least two.  Definitive for any network —
        no parametrisation or sign condition required.

        Parameters
        ----------
        initial_conditions : dict[str, float]
            Fixes the stoichiometric compatibility class (the conservation totals).
        backend : {'auto', 'sympy', 'phcpy'}
            Steady-state enumeration backend (see :meth:`all_steady_states`).

        Returns
        -------
        bool
            ``True`` ⇔ ≥ 2 positive steady states in this class.
        """
        states = self.all_steady_states(initial_conditions, backend=backend)
        return len(states) >= 2

    # ── Global stability (Horn–Jackson Lyapunov) ────────────────────────────────

    def is_complex_balanced(
        self, initial_conditions: dict[str, float] | None = None
    ) -> bool:
        """
        Whether the network is complex-balanced.

        Every weakly reversible deficiency-zero network is complex-balanced for
        *all* rate constants (Horn–Jackson). For weakly reversible higher-deficiency
        networks the property is rate-dependent and is verified numerically at the
        positive equilibrium — pass ``initial_conditions`` to enable that check.
        """
        from .stability import is_complex_balanced
        return is_complex_balanced(
            self._reactions, self.species, self._rates,
            self.deficiency, self.is_weakly_reversible,
            initial_conditions=initial_conditions,
            chemostatted_values=self._chemostatted or None,
        )

    def complex_balanced_equilibrium(
        self, initial_conditions: dict[str, float]
    ) -> dict[str, float]:
        """
        Return the **Birch point** — the unique positive complex-balanced equilibrium
        in the stoichiometric compatibility class fixed by ``initial_conditions``.

        Defined only for complex-balanced networks (every weakly reversible
        deficiency-zero network qualifies, for all rates). This equilibrium sets the
        means of the exact stochastic stationary distribution (see
        :meth:`stationary_distribution`). Raises ``ValueError`` otherwise.
        """
        from .stability import complex_balanced_equilibrium
        return complex_balanced_equilibrium(
            self._reactions, self.species, self._rates,
            self.deficiency, self.is_weakly_reversible, initial_conditions,
            chemostatted_values=self._chemostatted or None,
        )

    def certify_global_stability(
        self,
        initial_conditions: dict[str, float],
        n_samples: int = 200,
        seed: int = 0,
    ):
        """
        Certify *global* asymptotic stability via the Horn–Jackson Lyapunov function.

        For complex-balanced networks the pseudo-Helmholtz function
        ``V(c) = Σ_i [c_i(ln(c_i/c_i*) − 1) + c_i*]`` is a strict Lyapunov function,
        so the positive equilibrium ``c*`` is globally asymptotically stable within
        its stoichiometric compatibility class — a strictly stronger statement than
        the local eigenvalue test.

        Parameters
        ----------
        initial_conditions : dict[str, float]
            Fixes the compatibility class (conservation totals) and the equilibrium.
        n_samples : int
            Number of random points in the class at which ``dV/dt ≤ 0`` is verified.
        seed : int
            RNG seed for the sampling.

        Returns
        -------
        StabilityCertificate
            With ``.globally_stable``, ``.equilibrium``, ``.lyapunov_function`` (V),
            and ``.lyapunov_derivative`` (dV/dt).  ``str()`` gives a readable summary.
        """
        from .stability import certify_global_stability
        return certify_global_stability(
            self._reactions, self.species, self._rates,
            self.deficiency, self.is_weakly_reversible,
            self.odes(numeric_rates=True), initial_conditions,
            chemostatted_values=self._chemostatted or None,
            n_samples=n_samples, seed=seed,
        )

    # ── Interoperability (SBML) ─────────────────────────────────────────────────

    def to_sbml(self, path: str | None = None) -> str:
        """
        Export this network to SBML (Level 3 Version 2).

        Requires the optional ``python-libsbml`` package
        (``pip install mantis-delta[sbml]``).

        Parameters
        ----------
        path : str, optional
            If given, the SBML document is written to this file path.  The XML
            string is returned in all cases.

        Returns
        -------
        str
            The SBML document as an XML string.
        """
        from .sbml import network_to_sbml, write_sbml
        if path is not None:
            write_sbml(self, path)
        return network_to_sbml(self)

    @classmethod
    def from_sbml(cls, source: str) -> CRNetwork:
        """
        Construct a network from an SBML model.

        Requires the optional ``python-libsbml`` package
        (``pip install mantis-delta[sbml]``).

        Parameters
        ----------
        source : str
            Either a path to an ``.xml`` SBML file or a raw SBML/XML string.
            Mass-action kinetic laws are parsed to recover rate constants where
            possible; SBML boundary/constant species become chemostatted species.

        Returns
        -------
        CRNetwork
        """
        from .sbml import network_from_sbml
        return network_from_sbml(source)
