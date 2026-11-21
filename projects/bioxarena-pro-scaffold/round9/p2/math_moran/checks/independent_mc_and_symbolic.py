"""Checker 2: Independent SymPy symbolic proof check + forward-time Monte Carlo."""
import sys
sys.dont_write_bytecode = True

import math
import random
import sympy as sp

from fractions import Fraction
from common import formula_bd, formula_db, formula_mixed, require


def run_symbolic_checks() -> dict:
    """Verify all algebraic recurrence identities in Theorem 1 and Corollaries 1-3 via SymPy."""
    n, k, r, lam = sp.symbols("n k r lam", positive=True)
    H = sp.Function("H")

    def sym_t1(kk):
        return n * H(kk) + n * r * (n * kk - 1) / (kk + 1)

    def sym_t0(kk):
        return n * H(kk - 1) + (n + n**2 * r * (kk - 1)) / kk

    p_in = 1 / (n * (1 + r))
    p_out = r / (n * (1 + r))

    # 1. (k, 0) for general k >= 2
    res_k0 = sp.simplify(sym_t0(k) - (1 + (1 - k * p_in) * sym_t0(k) + k * p_in * sym_t1(k - 1)))
    require(res_k0 == 0, f"Symbolic residual (k,0) != 0: {res_k0}")

    # 2. (k, 1) for general 2 <= k <= n - 1
    rhs_k1 = (
        1
        + (n - k) * p_in * sym_t1(k)
        + k * p_in * sym_t1(k - 1)
        + k * p_out * sym_t0(k)
        + (n - k) * p_out * sym_t0(k + 1)
    )
    res_k1 = sp.simplify((sym_t1(k) - rhs_k1).subs(H(k), H(k - 1) + 1 / k))
    require(res_k1 == 0, f"Symbolic residual (k,1) != 0: {res_k1}")

    # 3. Boundary (n, 1) where (n - k) = 0
    rhs_n1 = 1 + n * p_in * sym_t1(n - 1) + n * p_out * sym_t0(n)
    res_n1 = sp.simplify((sym_t1(n) - rhs_n1).subs(H(n), H(n - 1) + 1 / n))
    require(res_n1 == 0, f"Symbolic residual (n,1) != 0: {res_n1}")

    # 4. Boundary (1, 1) for general n >= 2 (t1(0)=0, t0(1)=0, H(1)=1)
    rhs_11 = 1 + (n - 1) * p_in * sym_t1(1) + (n - 1) * p_out * sym_t0(2)
    res_11 = sp.simplify((sym_t1(1) - rhs_11).subs(H(1), 1))
    require(res_11 == 0, f"Symbolic residual (1,1) != 0: {res_11}")

    # 5. Boundary n = 1
    res_n1_base = sp.simplify(sym_t1(1).subs({n: 1, H(1): 1}) - 1)
    require(res_n1_base == 0, f"Symbolic residual n=1 != 0: {res_n1_base}")

    # 6. Telescoping difference identity
    diff_formula = sp.simplify((sym_t1(k) - sym_t1(k - 1)).subs(H(k), H(k - 1) + 1 / k))
    expected_diff = n / k + n * (n + 1) * r / (k * (k + 1))
    res_diff = sp.simplify(diff_formula - expected_diff)
    require(res_diff == 0, f"Telescoping difference mismatch: {res_diff}")

    # 7. Specializations: bd (r=n), db (r=1/n), mixed
    t_master = sp.simplify(sym_t1(n))
    t_bd = sp.simplify(t_master.subs(r, n))
    require(sp.simplify(t_bd - (n**3 - n**2 + n * H(n))) == 0, "bd specialization mismatch")

    t_db = sp.simplify(t_master.subs(r, 1 / n))
    require(sp.simplify(t_db - (n * H(n) + n - 1)) == 0, "db specialization mismatch")

    r_lam = (n * lam + 1 - lam) / (lam + n * (1 - lam))
    t_mix = sp.simplify(t_master.subs(r, r_lam))
    require(sp.simplify(t_mix.subs(lam, 1) - (n**3 - n**2 + n * H(n))) == 0, "mixed lam=1 mismatch")
    require(sp.simplify(t_mix.subs(lam, 0) - (n * H(n) + n - 1)) == 0, "mixed lam=0 mismatch")
    require(sp.simplify(t_mix.subs(lam, sp.Rational(1, 2)) - (n**2 - n + n * H(n))) == 0, "mixed lam=1/2 mismatch")

    # 8. Identity t_bd(n, 0) == t_bd(n, 1) for n >= 2
    res_n0_eq_n1 = sp.simplify((sym_t0(n) - sym_t1(n)).subs(r, n).subs(H(n), H(n - 1) + 1 / n))
    require(res_n0_eq_n1 == 0, f"t_bd(n,0) != t_bd(n,1): {res_n0_eq_n1}")

    return {
        "sympy_version": sp.__version__,
        "all_symbolic_residuals_zero": True,
        "verified_identities": [
            "(k,0) general k>=2",
            "(k,1) interior 2<=k<=n-1",
            "(n,1) upper boundary",
            "(1,1) lower boundary n>=2",
            "n=1 base case",
            "telescoping difference t(k,1)-t(k-1,1)",
            "bd, db, lambda-mixed specializations",
            "t_bd(n,0) == t_bd(n,1) for n>=2",
        ],
    }


def simulate_forward_star(n: int, lam: float, trials: int, seed: int) -> tuple[float, float]:
    """
    Literal forward-time multi-type Moran simulator on S_n (center=0, leaves=1..n).
    Zero duality, zero state compression, zero formula knowledge.
    Returns (sample_mean, sample_se).
    """
    rng = random.Random(seed)
    N = n + 1
    sum_t = 0.0
    sum_t2 = 0.0
    for _ in range(trials):
        types = list(range(N))
        freq = [1] * N
        distinct = N
        steps = 0
        while distinct > 1:
            steps += 1
            if lam == 1.0 or (lam > 0.0 and rng.random() < lam):
                # bd step: uniform birth, uniform neighbor death
                b = rng.randrange(N)
                if b == 0:
                    d = rng.randrange(1, N)
                else:
                    d = 0
            else:
                # db step: uniform death, uniform neighbor birth
                d = rng.randrange(N)
                if d == 0:
                    b = rng.randrange(1, N)
                else:
                    b = 0
            old_t = types[d]
            new_t = types[b]
            if old_t != new_t:
                types[d] = new_t
                freq[old_t] -= 1
                if freq[old_t] == 0:
                    distinct -= 1
                freq[new_t] += 1
        sum_t += steps
        sum_t2 += steps * steps
    mean = sum_t / trials
    var = max(0.0, (sum_t2 / trials) - mean * mean) * trials / max(1, trials - 1)
    se = math.sqrt(var / trials)
    return mean, se


def run_monte_carlo_checks(z_threshold: float = 4.0) -> list[dict]:
    """Run independent forward Monte Carlo checks across bd, db, and mixed updating."""
    cases = [
        # (rule_label, lam_frac, n, trials, seed)
        ("bd", Fraction(1), 2, 15000, 1002),
        ("bd", Fraction(1), 3, 15000, 1003),
        ("bd", Fraction(1), 4, 10000, 1004),
        ("bd", Fraction(1), 5, 10000, 1005),
        ("bd", Fraction(1), 8, 5000, 1008),
        ("bd", Fraction(1), 12, 3000, 1012),
        ("db", Fraction(0), 2, 15000, 2002),
        ("db", Fraction(0), 3, 15000, 2003),
        ("db", Fraction(0), 4, 15000, 2004),
        ("db", Fraction(0), 5, 15000, 2005),
        ("db", Fraction(0), 8, 10000, 2008),
        ("db", Fraction(0), 12, 10000, 2012),
        ("db", Fraction(0), 20, 10000, 2020),
        ("db", Fraction(0), 30, 10000, 2030),
        ("mixed_half", Fraction(1, 2), 3, 10000, 3003),
        ("mixed_half", Fraction(1, 2), 5, 10000, 3005),
        ("mixed_half", Fraction(1, 2), 10, 5000, 3010),
    ]
    results = []
    for label, lam, n, trials, seed in cases:
        mean, se = simulate_forward_star(n, float(lam), trials, seed)
        if lam == 1:
            exact = float(formula_bd(n))
        elif lam == 0:
            exact = float(formula_db(n))
        else:
            exact = float(formula_mixed(n, lam))
        z = (mean - exact) / se if se > 0 else 0.0
        rel_err = abs(mean - exact) / exact
        require(abs(z) < z_threshold, f"MC mismatch {label} n={n}: mean={mean:.4f}, exact={exact:.4f}, z={z:.2f}")
        results.append({
            "rule": label,
            "lam": str(lam),
            "n": n,
            "N": n + 1,
            "trials": trials,
            "seed": seed,
            "mc_mean": round(mean, 4),
            "mc_se": round(se, 4),
            "exact": round(exact, 6),
            "z_score": round(z, 3),
            "rel_err": round(rel_err, 5),
        })
    return results


def run_all_checker2() -> dict:
    sym = run_symbolic_checks()
    mc = run_monte_carlo_checks()
    return {"symbolic": sym, "monte_carlo": mc}


if __name__ == "__main__":
    out = run_all_checker2()
    max_z = max(abs(r["z_score"]) for r in out["monte_carlo"])
    print(f"Checker 2 PASS: SymPy symbolic residuals == 0; {len(out['monte_carlo'])} MC suites passed (max |z|={max_z:.3f}).")
