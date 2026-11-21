"""Adversarial / corrupted negative controls to verify both checkers detect errors."""
import sys
sys.dont_write_bytecode = True

import math
from fractions import Fraction
import sympy as sp

from common import formula_bd, formula_db, harmonic, require
from exact_markov_checker import exact_lumped_kc_check, exact_unreduced_forward_time
from independent_mc_and_symbolic import simulate_forward_star


def run_negative_controls() -> dict:
    controls = []

    # 1. Wrong formulas caught by Checker 1 (exact rational), SymPy, and Checker 2 (MC)
    wrong_formulas = [
        (
            "off_by_one_N_instead_of_n",
            lambda n: Fraction((n + 1) ** 3 - (n + 1) ** 2) + (n + 1) * harmonic(n + 1),
            lambda nn, HH: (nn + 1) ** 3 - (nn + 1) ** 2 + (nn + 1) * HH(nn + 1),
            5,
            5000,
        ),
        (
            "truncated_harmonic_H_{n-1}",
            lambda n: Fraction(n**3 - n**2) + n * harmonic(n - 1),
            lambda nn, HH: nn**3 - nn**2 + nn * HH(nn - 1),
            3,
            15000,
        ),
        (
            "plus_one_additive_corruption",
            lambda n: formula_bd(n) + 1,
            lambda nn, HH: nn**3 - nn**2 + nn * HH(nn) + 1,
            2,
            20000,
        ),
        (
            "log_instead_of_harmonic",
            None,
            None,
            3,
            20000,
        ),
    ]

    _, exact_n3 = exact_unreduced_forward_time(3, Fraction(1))
    _, exact_n2 = exact_unreduced_forward_time(2, Fraction(1))
    exact_by_n = {2: exact_n2, 3: exact_n3}

    n_sym = sp.symbols("n", positive=True)
    H = sp.Function("H")
    for name, q_fn, sym_fn, mc_n, mc_trials in wrong_formulas:
        # Check exact rational rejection at n=3 (and mc_n if in exact_by_n)
        if q_fn is not None:
            wrong_q3 = q_fn(3)
            require(wrong_q3 != exact_n3, f"Control {name} unexpectedly matched exact n=3")
            diff_n3 = str(wrong_q3 - exact_n3)
        else:
            wrong_f = mc_n**3 - mc_n**2 + mc_n * math.log(mc_n)
            exact_f = float(formula_bd(mc_n))
            require(abs(wrong_f - exact_f) > 2.0, "log control unexpectedly close")
            diff_n3 = f"{wrong_f - exact_f:.4f} (at n={mc_n})"

        # Check SymPy symbolic residual rejection (check whether expression matches formula_bd at n=3)
        if sym_fn is not None:
            sym_diff = sp.simplify(sym_fn(n_sym, H) - (n_sym**3 - n_sym**2 + n_sym * H(n_sym)))
            require(sym_diff != 0, f"Symbolic control {name} failed to trigger nonzero diff")

        # Check MC rejection (|z| > 4.0)
        mean, se = simulate_forward_star(mc_n, 1.0, mc_trials, seed=9000 + mc_n)
        if q_fn is not None:
            wrong_val = float(q_fn(mc_n))
        else:
            wrong_val = mc_n**3 - mc_n**2 + mc_n * math.log(mc_n)
        z = (mean - wrong_val) / se
        require(abs(z) > 4.0, f"MC failed to reject {name} at n={mc_n}: z={z:.2f}")
        controls.append({
            "control": name,
            "detected_by": ["exact_markov", "monte_carlo"] + (["sympy"] if sym_fn else []),
            "exact_discrepancy": diff_n3,
            "mc_n": mc_n,
            "mc_z_rejection": round(z, 2),
        })

    # 2. Self-loop model mismatch (adding self-loops violates Eq. (9))
    for loop_mode in ["center", "all"]:
        _, loop_val = exact_unreduced_forward_time(3, Fraction(1), self_loop_mode=loop_mode)
        require(loop_val > exact_n3, f"Self-loop mode {loop_mode} did not alter T_4")
        controls.append({
            "control": f"self_loop_{loop_mode}_n3",
            "detected_by": ["exact_unreduced_forward_markov"],
            "T_with_loops": str(loop_val),
            "T_simple_graph": str(exact_n3),
            "exact_discrepancy": str(loop_val - exact_n3),
        })

    # 3. Corrupted (k, c) transition matrix at n=4 (row sums preserved == 1)
    perturbed = exact_lumped_kc_check(4, Fraction(1), perturb_eps=Fraction(1, 100))
    exact_n4 = formula_bd(4)
    require(Fraction(perturbed["T_kc"]) != exact_n4, "Perturbed (k,c) kernel failed to alter T_5")
    controls.append({
        "control": "corrupted_kc_kernel_n4_eps_1_100",
        "detected_by": ["exact_lumped_kc_markov"],
        "T_perturbed": perturbed["T_kc"],
        "T_exact": str(exact_n4),
        "exact_discrepancy": str(Fraction(perturbed["T_kc"]) - exact_n4),
    })

    # 4. Swapped rule control (test db simulation against bd formula at n=4)
    mean_db4, se_db4 = simulate_forward_star(4, 0.0, 10000, seed=9104)
    z_swap = (mean_db4 - float(formula_bd(4))) / se_db4
    require(abs(z_swap) > 100.0, f"Swapped rule check did not strongly reject: z={z_swap}")
    require(formula_db(4) != formula_bd(4), "db and bd should differ at n=4")
    controls.append({
        "control": "swapped_bd_vs_db_n4",
        "detected_by": ["exact_markov", "monte_carlo"],
        "mc_db_mean": round(mean_db4, 4),
        "wrong_bd_exact": str(formula_bd(4)),
        "mc_z_rejection": round(z_swap, 2),
    })

    return {"negative_controls_passed": len(controls), "controls": controls}


if __name__ == "__main__":
    out = run_negative_controls()
    print(f"Negative controls PASS: {out['negative_controls_passed']} adversarial controls detected and rejected.")
