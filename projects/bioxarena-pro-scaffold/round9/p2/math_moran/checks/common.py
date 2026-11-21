"""Shared exact-arithmetic helpers and analytic formulas for star Moran checks."""
import sys
sys.dont_write_bytecode = True

from fractions import Fraction


def require(condition: bool, message: str) -> None:
    """Explicit runtime assertion retained even under python -O."""
    if not condition:
        raise AssertionError(message)


def harmonic(n: int) -> Fraction:
    """Exact nth harmonic number H_n = sum_{k=1}^n 1/k."""
    require(n >= 0, f"harmonic requires n >= 0, got {n}")
    return sum((Fraction(1, k) for k in range(1, n + 1)), Fraction(0))


def solve_linear_fraction(A: list[list[Fraction]], b: list[Fraction]) -> list[Fraction]:
    """Solve A x = b exactly over Q via Gaussian elimination with partial pivoting."""
    m = len(b)
    require(len(A) == m, "Row count mismatch")
    if m == 0:
        return []
    aug = [list(A[i]) + [b[i]] for i in range(m)]
    for col in range(m):
        piv = None
        for r in range(col, m):
            if aug[r][col] != 0:
                piv = r
                break
        require(piv is not None, f"Singular matrix at column {col}")
        if piv != col:
            aug[col], aug[piv] = aug[piv], aug[col]
        inv = Fraction(1) / aug[col][col]
        row_c = aug[col]
        for j in range(col, m + 1):
            row_c[j] *= inv
        for r in range(m):
            if r != col and aug[r][col] != 0:
                fac = aug[r][col]
                row_r = aug[r]
                for j in range(col, m + 1):
                    row_r[j] -= fac * row_c[j]
    return [aug[i][m] for i in range(m)]


def p_in_out_mixed(n: int, lam: Fraction) -> tuple[Fraction, Fraction]:
    """Return (p_in, p_out) for lambda-mixed updating on S_n (N = n + 1)."""
    require(n >= 1, f"n must be >= 1, got {n}")
    lam = Fraction(lam)
    require(0 <= lam <= 1, f"lam must be in [0, 1], got {lam}")
    N = n + 1
    p_in = lam * Fraction(1, N * n) + (1 - lam) * Fraction(1, N)
    p_out = lam * Fraction(1, N) + (1 - lam) * Fraction(1, N * n)
    return p_in, p_out


def formula_t1(n: int, k: int, r: Fraction) -> Fraction:
    """Exact expected coalescing time t(k, 1) for 0 <= k <= n."""
    require(0 <= k <= n, f"Invalid k={k} for n={n}")
    if k == 0:
        return Fraction(0)
    return n * harmonic(k) + Fraction(n * r * (n * k - 1), k + 1)


def formula_t0(n: int, k: int, r: Fraction) -> Fraction:
    """Exact expected coalescing time t(k, 0) for 1 <= k <= n."""
    require(1 <= k <= n, f"Invalid k={k} for n={n}")
    if k == 1:
        return Fraction(0)
    return n * harmonic(k - 1) + Fraction(n + n * n * r * (k - 1), k)


def formula_bd(n: int) -> Fraction:
    """Conjecture Eq. (9): T_N^bd = n^3 - n^2 + n * H_n."""
    require(n >= 1, f"n must be >= 1, got {n}")
    return Fraction(n**3 - n**2) + n * harmonic(n)


def formula_db(n: int) -> Fraction:
    """Exact db consensus time: T_N^db = n * H_n + n - 1."""
    require(n >= 1, f"n must be >= 1, got {n}")
    return n * harmonic(n) + Fraction(n - 1)


def formula_mixed(n: int, lam: Fraction) -> Fraction:
    """Exact lambda-mixed consensus time on S_n."""
    p_in, p_out = p_in_out_mixed(n, lam)
    r = p_out / p_in
    return n * harmonic(n) + n * (n - 1) * r
