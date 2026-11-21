"""Checker 1: Independent exact rational Markov chain solvers on S_n."""
import sys
sys.dont_write_bytecode = True

from collections import deque
from fractions import Fraction

from common import (
    formula_bd,
    formula_db,
    formula_mixed,
    formula_t0,
    formula_t1,
    p_in_out_mixed,
    require,
    solve_linear_fraction,
)


def solve_absorption_from_adj(
    transient_states: list,
    adj: dict,
    absorbing_check,
) -> dict:
    """Given adj[s] = {next_s: prob}, solve (I - P_trans) tau = 1 exactly over Q."""
    idx = {s: i for i, s in enumerate(transient_states)}
    m = len(transient_states)
    A = [[Fraction(0)] * m for _ in range(m)]
    b = [Fraction(1)] * m
    for i, s in enumerate(transient_states):
        row_sum = sum(adj[s].values(), Fraction(0))
        require(row_sum == 1, f"Row sum != 1 at state {s}: {row_sum}")
        require(all(p >= 0 for p in adj[s].values()), f"Negative prob at state {s}")
        A[i][i] = Fraction(1)
        for ns, p in adj[s].items():
            if not absorbing_check(ns):
                require(ns in idx, f"Unknown non-absorbing state {ns}")
                A[i][idx[ns]] -= p
    sol = solve_linear_fraction(A, b)
    return {s: sol[i] for i, s in enumerate(transient_states)}


def verify_reachability(all_states: set, adj: dict, absorbing_check) -> None:
    """Verify every reachable state has a positive-probability path to absorption."""
    absorbing = {s for s in all_states if absorbing_check(s)}
    rev = {s: set() for s in all_states}
    for s, row in adj.items():
        for ns, p in row.items():
            if p > 0:
                rev[ns].add(s)
    reached = set(absorbing)
    dq = deque(absorbing)
    while dq:
        u = dq.popleft()
        for prev in rev[u]:
            if prev not in reached:
                reached.add(prev)
                dq.append(prev)
    require(reached == all_states, "Some states cannot reach absorption")


def exact_unreduced_forward_time(n: int, lam: Fraction, self_loop_mode: str = "none") -> tuple[int, Fraction]:
    """
    Literal N^N vertex-labeled forward Markov chain directly from Brewster et al. Eq. (15).
    Zero type-relabeling and zero symmetry reduction.
    self_loop_mode in {'none', 'center', 'all'} for normal checks and negative controls.
    """
    N = n + 1
    init = tuple(range(N))
    p_in, p_out = p_in_out_mixed(n, lam)
    adj = {}
    seen = {init}
    dq = deque([init])
    while dq:
        s = dq.popleft()
        if len(set(s)) == 1:
            adj[s] = {s: Fraction(1)}
            continue
        trans = {}
        if self_loop_mode == "none":
            for j in range(1, N):
                ns = list(s)
                ns[j] = s[0]
                ns_t = tuple(ns)
                trans[ns_t] = trans.get(ns_t, Fraction(0)) + p_in
                if ns_t not in seen:
                    seen.add(ns_t)
                    dq.append(ns_t)
            for i in range(1, N):
                ns = list(s)
                ns[0] = s[i]
                ns_t = tuple(ns)
                trans[ns_t] = trans.get(ns_t, Fraction(0)) + p_out
                if ns_t not in seen:
                    seen.add(ns_t)
                    dq.append(ns_t)
        elif self_loop_mode == "center":
            # Add a self-loop at center only (center degree becomes n + 1)
            require(lam == 1, "center self-loop mode defined for bd")
            p_in_loop = Fraction(1, N * (n + 1))
            trans[s] = p_in_loop
            for j in range(1, N):
                ns = list(s)
                ns[j] = s[0]
                ns_t = tuple(ns)
                trans[ns_t] = trans.get(ns_t, Fraction(0)) + p_in_loop
                if ns_t not in seen:
                    seen.add(ns_t)
                    dq.append(ns_t)
            for i in range(1, N):
                ns = list(s)
                ns[0] = s[i]
                ns_t = tuple(ns)
                trans[ns_t] = trans.get(ns_t, Fraction(0)) + Fraction(1, N)
                if ns_t not in seen:
                    seen.add(ns_t)
                    dq.append(ns_t)
        elif self_loop_mode == "all":
            # Add self-loop at every vertex (center deg n+1, leaf deg 2)
            require(lam == 1, "all self-loop mode defined for bd")
            p_in_loop = Fraction(1, N * (n + 1))
            p_out_loop = Fraction(1, N * 2)
            trans[s] = p_in_loop + n * p_out_loop
            for j in range(1, N):
                ns = list(s)
                ns[j] = s[0]
                ns_t = tuple(ns)
                trans[ns_t] = trans.get(ns_t, Fraction(0)) + p_in_loop
                if ns_t not in seen:
                    seen.add(ns_t)
                    dq.append(ns_t)
            for i in range(1, N):
                ns = list(s)
                ns[0] = s[i]
                ns_t = tuple(ns)
                trans[ns_t] = trans.get(ns_t, Fraction(0)) + p_out_loop
                if ns_t not in seen:
                    seen.add(ns_t)
                    dq.append(ns_t)
        else:
            raise ValueError(f"Unknown self_loop_mode: {self_loop_mode}")
        adj[s] = trans

    is_abs = lambda st: len(set(st)) == 1
    verify_reachability(seen, adj, is_abs)
    transient = sorted(s for s in seen if not is_abs(s))
    tau = solve_absorption_from_adj(transient, adj, is_abs)
    return len(transient), tau[init]


def exact_canonical_forward_time(n: int, lam: Fraction) -> tuple[int, Fraction]:
    """
    Forward Markov chain reduced ONLY by type-label isomorphism (order of first appearance).
    Leaves 1..n are still tracked individually with zero leaf-symmetry lumping.
    """
    N = n + 1

    def canon(s: list[int] | tuple[int, ...]) -> tuple[int, ...]:
        mapping = {}
        out = []
        for x in s:
            if x not in mapping:
                mapping[x] = len(mapping)
            out.append(mapping[x])
        return tuple(out)

    init = tuple(range(N))
    p_in, p_out = p_in_out_mixed(n, lam)
    adj = {}
    seen = {init}
    dq = deque([init])
    while dq:
        s = dq.popleft()
        if len(set(s)) == 1:
            adj[s] = {s: Fraction(1)}
            continue
        trans = {}
        for j in range(1, N):
            ns = list(s)
            ns[j] = s[0]
            cs = canon(ns)
            trans[cs] = trans.get(cs, Fraction(0)) + p_in
            if cs not in seen:
                seen.add(cs)
                dq.append(cs)
        for i in range(1, N):
            ns = list(s)
            ns[0] = s[i]
            cs = canon(ns)
            trans[cs] = trans.get(cs, Fraction(0)) + p_out
            if cs not in seen:
                seen.add(cs)
                dq.append(cs)
        adj[s] = trans

    is_abs = lambda st: len(set(st)) == 1
    verify_reachability(seen, adj, is_abs)
    transient = sorted(s for s in seen if not is_abs(s))
    tau = solve_absorption_from_adj(transient, adj, is_abs)
    return len(transient), tau[init]


def exact_subset_coalescing_check(n: int, lam: Fraction) -> dict:
    """
    Build full (2^(n+1)-1)-state subset coalescing walk S_m on V = {0..n},
    verify exact strong lumpability under pi(S) = (|S intersect [n]|, 1[0 in S]),
    solve the unlumped 2^N-1 linear system, and check agreement at every subset S.
    """
    N = n + 1
    p_in, p_out = p_in_out_mixed(n, lam)
    r = p_out / p_in
    all_subsets = list(range(1, 1 << N))  # bitmask; bit 0 is center, bits 1..n are leaves

    def pi(mask: int) -> tuple[int, int]:
        c = mask & 1
        k = (mask >> 1).bit_count()
        return (k, c)

    def kc_expected_row(kc: tuple[int, int]) -> dict[tuple[int, int], Fraction]:
        k, c = kc
        if k + c == 1:
            return {kc: Fraction(1)}
        row = {}
        if c == 1:
            if k > 0:
                row[(k - 1, 1)] = k * p_in
            if n - k > 0:
                row[(k, 1)] = (n - k) * p_in
                row[(k + 1, 0)] = (n - k) * p_out
            row[(k, 0)] = k * p_out
        else:
            row[(k - 1, 1)] = k * p_in
            row[(k, 0)] = 1 - k * p_in
        return {dest: prob for dest, prob in row.items() if prob > 0}

    adj = {}
    for mask in all_subsets:
        if mask.bit_count() == 1:
            adj[mask] = {mask: Fraction(1)}
            continue
        trans = {}
        # Edge 0 -> i for i in 1..n: if i in S, replace i by 0 in S
        for i in range(1, N):
            bit_i = 1 << i
            if mask & bit_i:
                nmask = (mask ^ bit_i) | 1
            else:
                nmask = mask
            trans[nmask] = trans.get(nmask, Fraction(0)) + p_in
        # Edge i -> 0 for i in 1..n: if 0 in S, replace 0 by i in S
        for i in range(1, N):
            bit_i = 1 << i
            if mask & 1:
                nmask = (mask ^ 1) | bit_i
            else:
                nmask = mask
            trans[nmask] = trans.get(nmask, Fraction(0)) + p_out
        adj[mask] = trans

        # Strong lumpability check at this subset mask
        lumped_agg = {}
        for nmask, prob in trans.items():
            dest_kc = pi(nmask)
            lumped_agg[dest_kc] = lumped_agg.get(dest_kc, Fraction(0)) + prob
        expected_kc = kc_expected_row(pi(mask))
        require(lumped_agg == expected_kc, f"Lumpability failed at mask={bin(mask)}: {lumped_agg} vs {expected_kc}")

    is_abs = lambda m: m.bit_count() == 1
    verify_reachability(set(all_subsets), adj, is_abs)
    transient = [m for m in all_subsets if not is_abs(m)]
    tau = solve_absorption_from_adj(transient, adj, is_abs)
    for m in transient:
        k, c = pi(m)
        expected = formula_t1(n, k, r) if c == 1 else formula_t0(n, k, r)
        require(tau[m] == expected, f"Subset mismatch n={n}, mask={bin(m)}: {tau[m]} vs {expected}")
    full_mask = (1 << N) - 1
    return {
        "n": n,
        "lam": str(lam),
        "subset_states": len(transient),
        "T_subset": str(tau[full_mask]),
    }


def exact_lumped_kc_check(n: int, lam: Fraction, perturb_eps: Fraction = Fraction(0)) -> dict:
    """
    Build the (2n-1)-state (k, c) chain from Lemma 2, solve (I - P)t = 1 over Q,
    and verify every single entry against formula_t1 / formula_t0.
    If perturb_eps != 0, inject a row-sum-preserving transition perturbation (for negative controls).
    """
    p_in, p_out = p_in_out_mixed(n, lam)
    r = p_out / p_in
    transient = [(k, 1) for k in range(1, n + 1)] + [(k, 0) for k in range(2, n + 1)]
    all_kc = set(transient) | {(0, 1), (1, 0)}
    adj = {(0, 1): {(0, 1): Fraction(1)}, (1, 0): {(1, 0): Fraction(1)}}

    for k in range(1, n + 1):
        row = {}
        row[(k - 1, 1)] = k * p_in
        if n - k > 0:
            row[(k, 1)] = (n - k) * p_in
            row[(k + 1, 0)] = (n - k) * p_out
        row[(k, 0)] = k * p_out
        adj[(k, 1)] = {st: pr for st, pr in row.items() if pr > 0}

    for k in range(2, n + 1):
        adj[(k, 0)] = {
            (k - 1, 1): k * p_in,
            (k, 0): 1 - k * p_in,
        }

    if perturb_eps != 0:
        require(n >= 3, "Perturbation requires n >= 3")
        # Move perturb_eps from (n, 0)->(n, 0) to (n, 0)->(n-1, 1), keeping row sum == 1
        adj[(n, 0)][(n, 0)] -= perturb_eps
        adj[(n, 0)][(n - 1, 1)] += perturb_eps

    is_abs = lambda kc: kc[0] + kc[1] == 1
    verify_reachability(all_kc, adj, is_abs)
    tau = solve_absorption_from_adj(transient, adj, is_abs)

    if perturb_eps == 0:
        for k in range(1, n + 1):
            require(tau[(k, 1)] == formula_t1(n, k, r), f"t({k},1) mismatch for n={n}, lam={lam}")
        for k in range(2, n + 1):
            require(tau[(k, 0)] == formula_t0(n, k, r), f"t({k},0) mismatch for n={n}, lam={lam}")
        if lam == 1 and n >= 2:
            require(tau[(n, 0)] == tau[(n, 1)], f"t_bd(n,0) != t_bd(n,1) at n={n}")

    return {
        "n": n,
        "lam": str(lam),
        "kc_transient": len(transient),
        "T_kc": str(tau[(n, 1)]),
    }


def run_all_exact_checks() -> dict:
    """Run the complete Checker 1 suite and return structured results."""
    unreduced_results = []
    for lam, fn in [(Fraction(1), formula_bd), (Fraction(0), formula_db), (Fraction(1, 2), lambda n: formula_mixed(n, Fraction(1, 2)))]:
        for n in [1, 2, 3]:
            st, val = exact_unreduced_forward_time(n, lam)
            exp = fn(n)
            require(val == exp, f"Unreduced forward mismatch n={n}, lam={lam}: {val} vs {exp}")
            unreduced_results.append({"n": n, "lam": str(lam), "transient_states": st, "T_exact": str(val)})

    canonical_results = []
    for lam, fn in [
        (Fraction(1), formula_bd),
        (Fraction(0), formula_db),
        (Fraction(1, 3), lambda n: formula_mixed(n, Fraction(1, 3))),
        (Fraction(1, 2), lambda n: formula_mixed(n, Fraction(1, 2))),
    ]:
        for n in [1, 2, 3, 4, 5]:
            st, val = exact_canonical_forward_time(n, lam)
            exp = fn(n)
            require(val == exp, f"Canonical forward mismatch n={n}, lam={lam}: {val} vs {exp}")
            canonical_results.append({"n": n, "lam": str(lam), "transient_states": st, "T_exact": str(val)})

    subset_results = []
    for lam in [Fraction(1), Fraction(0), Fraction(1, 2)]:
        for n in range(1, 7):
            subset_results.append(exact_subset_coalescing_check(n, lam))

    kc_results = []
    for lam in [Fraction(1), Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(3, 4)]:
        for n in range(1, 21):
            kc_results.append(exact_lumped_kc_check(n, lam))

    return {
        "unreduced_forward_checks": unreduced_results,
        "canonical_forward_checks": canonical_results,
        "subset_coalescing_checks": subset_results,
        "lumped_kc_checks_count": len(kc_results),
        "lumped_kc_bd_sample": [r for r in kc_results if r["lam"] == "1" and r["n"] in (1, 2, 3, 4, 5, 10, 15, 20)],
        "lumped_kc_db_sample": [r for r in kc_results if r["lam"] == "0" and r["n"] in (1, 2, 3, 4, 5, 10, 15, 20)],
    }


if __name__ == "__main__":
    out = run_all_exact_checks()
    print(f"Checker 1 PASS: {len(out['unreduced_forward_checks'])} unreduced, "
          f"{len(out['canonical_forward_checks'])} canonical, "
          f"{len(out['subset_coalescing_checks'])} subset/lumpability, "
          f"{out['lumped_kc_checks_count']} lumped (k,c) solves.")
