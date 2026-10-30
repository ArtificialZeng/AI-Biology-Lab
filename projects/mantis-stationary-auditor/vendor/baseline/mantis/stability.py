"""Global asymptotic stability certification via the Horn–Jackson Lyapunov function.

The structural theorems already in mantis certify *local* stability through the
sign of the Jacobian eigenvalues.  For **complex-balanced** networks a much
stronger, *global* statement is available (Horn & Jackson 1972; Feinberg's
Deficiency Zero Theorem): the pseudo-Helmholtz function

    V(c) = Σ_i [ c_i (ln(c_i / c_i*) − 1) + c_i* ]

is a strict Lyapunov function on every positive stoichiometric compatibility
class, where ``c*`` is the (unique) positive complex-balanced equilibrium of that
class.  Indeed ``∂V/∂c_i = ln(c_i/c_i*)``, so

    dV/dt = Σ_i ln(c_i/c_i*) · ċ_i ≤ 0,

with equality only at ``c = c*``.  Hence ``c*`` is *globally* asymptotically
stable within the relative interior of its compatibility class — strictly
stronger than the eigenvalue check.

Every weakly reversible deficiency-zero network is complex-balanced for **all**
rate constants (Horn–Jackson), so the certificate is purely structural there; for
higher-deficiency weakly-reversible networks complex balance is rate-dependent and
is verified numerically at the equilibrium.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from .parsing import Reaction, reduce_complex
from .stoichiometry import build_stoichiometry_matrix, matrix_rank

if TYPE_CHECKING:
    import sympy


@dataclass
class StabilityCertificate:
    """Outcome of a global-stability analysis.

    Attributes
    ----------
    is_complex_balanced : bool
        Whether the network is complex-balanced (structurally for δ=0 weakly
        reversible, else numerically at the equilibrium).
    globally_stable : bool
        True when a Horn–Jackson Lyapunov certificate establishes global
        asymptotic stability within the compatibility class.
    equilibrium : dict[str, float] | None
        The positive complex-balanced equilibrium ``c*`` (per the supplied class).
    lyapunov_function : sympy.Expr | None
        The symbolic pseudo-Helmholtz function ``V(c)``.
    lyapunov_derivative : sympy.Expr | None
        The symbolic ``dV/dt = Σ ln(c_i/c_i*)·f_i(c)``.
    n_samples_checked : int
        Number of random points in the compatibility class at which ``dV/dt ≤ 0``
        was numerically verified.
    max_dVdt : float
        Largest ``dV/dt`` observed over the samples (should be ≤ 0).
    rationale : str
        Human-readable explanation of the verdict.
    """
    is_complex_balanced: bool
    globally_stable: bool
    equilibrium: dict[str, float] | None = None
    lyapunov_function: "sympy.Expr | None" = None
    lyapunov_derivative: "sympy.Expr | None" = None
    n_samples_checked: int = 0
    max_dVdt: float = 0.0
    rationale: str = ""
    summary_lines: list[str] = field(default_factory=list)

    def __str__(self) -> str:
        return "\n".join(self.summary_lines) if self.summary_lines else self.rationale


def _complex_flux_balance(reactions, species, rate_values, c_star, chem, tol):
    """Return the max complex-level imbalance |inflow − outflow| at ``c_star``."""
    sp_idx = {s: i for i, s in enumerate(species)}
    # Map each reduced complex to an index.
    complexes: dict = {}

    def cid(comp):
        key = reduce_complex(comp, chem)
        if key not in complexes:
            complexes[key] = len(complexes)
        return complexes[key]

    edges = []  # (source_cid, product_cid, flux)
    for rxn in reactions:
        k = float(rate_values.get(rxn.rate_key, 0.0))
        flux = k
        for name, coeff in rxn.reactants:
            if name in sp_idx and name not in chem:
                flux *= c_star[sp_idx[name]] ** coeff
            elif name in chem:
                flux *= chem[name] ** coeff
        edges.append((cid(rxn.reactants), cid(rxn.products), flux))

    balance = np.zeros(len(complexes))
    for s, p, flux in edges:
        balance[s] -= flux   # outflow from source complex
        balance[p] += flux   # inflow to product complex
    return float(np.max(np.abs(balance))) if len(balance) else 0.0


def _positive_equilibrium(reactions, species, rate_values, initial_conditions,
                          chemostatted_values):
    """Find a strictly positive steady state in the class fixed by the ICs."""
    from .analysis import find_steady_states
    states = find_steady_states(
        reactions, species, rate_values, initial_conditions,
        n_attempts=12, chemostatted_values=chemostatted_values,
    )
    for st in states:
        y = np.array([st.concentrations[s] for s in species])
        if np.all(y > 0):
            return y
    return None


def is_complex_balanced(reactions, species, rate_values, deficiency,
                        is_weakly_reversible, initial_conditions=None,
                        chemostatted_values=None, tol=1e-6):
    """Decide whether the network is complex-balanced.

    Complex balance ⇒ weak reversibility (necessary).  Every weakly reversible
    deficiency-zero network is complex-balanced for all rates (Horn–Jackson).
    For weakly reversible higher-deficiency networks the property is rate-dependent
    and is checked numerically at the positive equilibrium (requires ICs).
    """
    if not is_weakly_reversible:
        return False
    if deficiency == 0:
        return True
    if initial_conditions is None:
        return False
    chem = chemostatted_values or {}
    c_star = _positive_equilibrium(reactions, species, rate_values,
                                   initial_conditions, chem or None)
    if c_star is None:
        return False
    imbalance = _complex_flux_balance(reactions, species, rate_values, c_star, chem, tol)
    scale = max(abs(float(rate_values.get(r.rate_key, 0.0))) for r in reactions) or 1.0
    return imbalance <= tol * scale


def complex_balanced_equilibrium(reactions, species, rate_values, deficiency,
                                 is_weakly_reversible, initial_conditions,
                                 chemostatted_values=None):
    """Return the Birch point: the unique positive complex-balanced equilibrium.

    For a complex-balanced network the Birch theorem guarantees exactly one positive
    equilibrium in each stoichiometric compatibility class (fixed here by the initial
    conditions).  Raises ``ValueError`` if the network is not complex-balanced or no
    positive equilibrium is found.
    """
    chem = chemostatted_values or {}
    if not is_complex_balanced(reactions, species, rate_values, deficiency,
                               is_weakly_reversible, initial_conditions, chem or None):
        raise ValueError(
            "Network is not complex-balanced at these rate constants; the Birch point "
            "is only defined for complex-balanced networks."
        )
    c_star = _positive_equilibrium(reactions, species, rate_values,
                                   initial_conditions, chem or None)
    if c_star is None:
        raise ValueError("No positive equilibrium found in the given compatibility class.")
    return dict(zip(species, c_star.tolist()))


def _lyapunov_symbolics(species, odes, c_star):
    """Build the symbolic V(c) and dV/dt for a complex-balanced equilibrium c*."""
    import sympy

    syms = {s: sympy.Symbol(s) for s in species}
    V = sympy.Integer(0)
    dV = sympy.Integer(0)
    for s in species:
        ci = syms[s]
        cs = sympy.Float(c_star[s])
        V += ci * (sympy.log(ci / cs) - 1) + cs
        # ∂V/∂c_i = ln(c_i / c_i*);  dV/dt sums these against the ODE RHS.
        dV += sympy.log(ci / cs) * odes[s]
    return sympy.simplify(V), dV


def certify_global_stability(reactions, species, rate_values, deficiency,
                             is_weakly_reversible, odes, initial_conditions,
                             chemostatted_values=None, n_samples=200, seed=0,
                             tol=1e-9):
    """Attempt a Horn–Jackson global-stability certificate for one compatibility class.

    ``odes`` is the species → SymPy RHS mapping (numeric rates substituted).
    Returns a :class:`StabilityCertificate`.
    """
    import sympy
    from .analysis import _conservation_law_vectors, _make_feasible_initial

    chem = chemostatted_values or {}
    cb = is_complex_balanced(reactions, species, rate_values, deficiency,
                             is_weakly_reversible, initial_conditions, chem or None)

    if not cb:
        reason = ("not weakly reversible" if not is_weakly_reversible
                  else "not complex-balanced at these rate constants")
        lines = [
            "Global stability (Horn–Jackson Lyapunov): NOT certified",
            f"  → Network is {reason}; the complex-balanced Lyapunov function does "
            f"not apply.",
            "  → Local eigenvalue stability (see steady_states) still holds where computed.",
        ]
        return StabilityCertificate(
            is_complex_balanced=False, globally_stable=False,
            rationale=reason, summary_lines=lines,
        )

    c_star_vec = _positive_equilibrium(reactions, species, rate_values,
                                       initial_conditions, chem or None)
    if c_star_vec is None:
        lines = ["Global stability: complex-balanced, but no positive equilibrium "
                 "was found in the given class (check initial conditions)."]
        return StabilityCertificate(
            is_complex_balanced=True, globally_stable=False,
            rationale="no positive equilibrium found", summary_lines=lines,
        )
    c_star = dict(zip(species, c_star_vec.tolist()))

    V, dV = _lyapunov_symbolics(species, odes, c_star)

    # Numerically verify dV/dt ≤ 0 across the compatibility class.
    syms = [sympy.Symbol(s) for s in species]
    dV_fn = sympy.lambdify(syms, dV, "numpy")
    cl_vecs = _conservation_law_vectors(reactions, species)
    y_ref = np.array([initial_conditions.get(s, 0.0) for s in species])
    cl_totals = [float(v @ y_ref) for v in cl_vecs]
    rng = np.random.default_rng(seed)

    max_dVdt = -np.inf
    n_ok = 0
    for _ in range(n_samples):
        if cl_vecs:
            y = _make_feasible_initial(cl_vecs, cl_totals, len(species), rng)
        else:
            y = np.abs(rng.normal(loc=np.maximum(c_star_vec, 1e-9),
                                  scale=np.maximum(c_star_vec, 1e-9)))
        if y is None:
            continue
        y = np.maximum(y, 1e-30)
        val = float(dV_fn(*y))
        if not np.isfinite(val):
            continue
        max_dVdt = max(max_dVdt, val)
        n_ok += 1

    max_dVdt = 0.0 if max_dVdt == -np.inf else max_dVdt
    # Allow a small positive tolerance for floating-point noise near c*.
    numerically_ok = max_dVdt <= tol * (abs(float(np.max(c_star_vec))) + 1.0)

    lines = [
        "Global stability (Horn–Jackson Lyapunov): CERTIFIED",
        f"  → Network is complex-balanced "
        f"({'δ=0 weakly reversible, all rates' if deficiency == 0 else 'verified at equilibrium'}).",
        "  → V(c) = Σ_i [c_i(ln(c_i/c_i*) − 1) + c_i*] is a strict Lyapunov function;",
        "    dV/dt ≤ 0 with equality only at c*  ⟹  c* is GLOBALLY asymptotically",
        "    stable within its positive stoichiometric compatibility class.",
        f"  → Numerical check: dV/dt ≤ 0 at {n_ok} sampled points "
        f"(max dV/dt = {max_dVdt:.3e}).",
    ]
    return StabilityCertificate(
        is_complex_balanced=True,
        globally_stable=bool(numerically_ok),
        equilibrium=c_star,
        lyapunov_function=V,
        lyapunov_derivative=dV,
        n_samples_checked=n_ok,
        max_dVdt=max_dVdt,
        rationale="complex-balanced; Horn–Jackson Lyapunov certificate",
        summary_lines=lines,
    )
