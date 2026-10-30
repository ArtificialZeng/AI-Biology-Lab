from .network import CRNetwork
from .parsing import Reaction, Complex
from .analysis import (
    SteadyState, BifurcationResult, SimulationResult,
    StochasticResult, gillespie_simulate, tau_leap_simulate,
)
from .crnt import CRNTResult
from .steady_states import all_steady_states
from .stability import (
    StabilityCertificate, is_complex_balanced, certify_global_stability,
    complex_balanced_equilibrium,
)
from .acr import ACRResult, detect_acr
from .injectivity import InjectivityResult, test_injectivity
from .stochastic_stationary import StationaryDistribution, stationary_distribution
from .fsp import FSPResult, fsp_solve
from .continuation import (
    ContinuationResult, BifurcationPoint, pseudo_arclength_continuation,
)
from .multistationarity import MultistationarityResult, multistationarity_region

__all__ = [
    "CRNetwork", "Reaction", "Complex",
    "SteadyState", "BifurcationResult", "SimulationResult",
    "StochasticResult", "gillespie_simulate", "tau_leap_simulate",
    "CRNTResult", "all_steady_states",
    "StabilityCertificate", "is_complex_balanced", "certify_global_stability",
    "complex_balanced_equilibrium",
    "ACRResult", "detect_acr",
    "InjectivityResult", "test_injectivity",
    "StationaryDistribution", "stationary_distribution",
    "FSPResult", "fsp_solve",
    "ContinuationResult", "BifurcationPoint", "pseudo_arclength_continuation",
    "MultistationarityResult", "multistationarity_region",
]
from importlib.metadata import PackageNotFoundError, version as _pkg_version

try:
    # Derived from git tags at build time via setuptools-scm.
    __version__ = _pkg_version("mantis-delta")
except PackageNotFoundError:  # source tree that was never installed
    __version__ = "0.0.0+unknown"
