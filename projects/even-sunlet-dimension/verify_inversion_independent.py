"""Fresh referee diagnostic; imports no owner code. Python 3 + python-flint.

build_certificate() reconstructs the two source monomials on ALL raw group
labels. Exact ranks over Q are finite diagnostics; the report supplies the
all-q argument. No numerical tolerance and no floating point are used.
"""
import hashlib
import json
from flint import fmpq_mat

def build_certificate():
    cases = []
    for q in (3, 4, 5, 6, 7, 8, 11, 12, 17, 18, 23, 24):
        m = q // 2
        k = m + 1
        n = 6 * k
        fold = lambda v: min(v % q, q - v % q)
        weight = [1] * k + [1] * k + [1] * k + [4 * t + 1 for t in range(k)] + [2] * k + [1] + [3] * m
        require(min(weight) > 0, 'Independent invariant, original line 21')
        point = [2 + (i * i + 13 * i + 7 * q) % 17 for i in range(n)]
        selected_rows, ordinary_rows, serialized = ([], [], [])
        margins, branches = ([], [0, 0])
        for x in range(q):
            for y in range(q):
                a, b, c = (fold(x), fold(y), fold(-x - y))
                e = [0] * n
                f = [0] * n
                for pos in (a, k + b, 2 * k + c, 3 * k + fold(x + y), 4 * k + a):
                    e[pos] += 1
                for pos in (a, k + b, 2 * k + c, 3 * k + b, 5 * k + a):
                    f[pos] += 1
                v1 = sum((z * t for z, t in zip(e, weight)))
                v2 = sum((z * t for z, t in zip(f, weight)))
                require(v1 != v2, 'Independent invariant, original line 37')
                branch = int(v2 < v1)
                require(branch == 1 if x == 0 else branch == int(c > b), 'Independent invariant, original line 39')
                require(v1 - v2 == 4 * (c - b) - (-1 if a == 0 else 1), 'Independent invariant, original line 40')
                selected_rows.append((e, f)[branch])
                margins.append(abs(v1 - v2))
                branches[branch] += 1
                derivative = [0] * n
                for powers in (e, f):
                    value = 1
                    for j, power in enumerate(powers):
                        value *= point[j] ** power
                    for j, power in enumerate(powers):
                        if power:
                            derivative[j] += power * value // point[j]
                ordinary_rows.append(derivative)
                serialized.append([x, y, e, f, v1, v2, branch])
        kernel = []
        for free in range(k + 4):
            coeff = [0] * n
            A = [int(free == i) for i in range(k)]
            B, C, D = (int(free == k + j) for j in range(3))
            K = -B - C - D
            E0 = int(free == k + 3)
            for i in range(k):
                coeff[i] = A[i]
                coeff[k + i] = B
                coeff[2 * k + i] = C
                coeff[3 * k + i] = D
                coeff[4 * k + i] = K - A[i] if i else E0
                coeff[5 * k + i] = K - A[i]
            kernel.append(coeff)
        matrix = fmpq_mat(selected_rows)
        null_vectors = fmpq_mat(kernel)
        require(matrix * null_vectors.transpose() == fmpq_mat(q * q, k + 4), 'Independent invariant, original line 73')
        require(null_vectors.rank() == k + 4, 'Independent invariant, original line 74')
        selected_rank = matrix.rank()
        jacobian_rank = fmpq_mat(ordinary_rows).rank()
        expected = 5 * k - 4 if q >= 5 else {3: 5, 4: 10}[q]
        require(selected_rank == jacobian_rank == expected, 'Independent invariant, original line 78')
        cases.append({'q': q, 'k': k, 'raw_labels': q * q, 'selected_shape': [q * q, n], 'weights': weight, 'smallest_strict_gap': min(margins), 'first_second_column_counts': branches, 'selected_rank_Q': selected_rank, 'integer_kernel_dimension_Q': null_vectors.rank(), 'integer_kernel_annihilation_exact': True, 'ordinary_jacobian_point': point, 'ordinary_jacobian_rank_Q': jacobian_rank, 'proof_upper_dimension': 5 * k - 4 if q >= 5 else None, 'fresh_held_out': q in (23, 24), 'literal_data_sha256': hashlib.sha256(json.dumps(serialized, separators=(',', ':')).encode()).hexdigest()})
    return {'checker': 'independently written fresh referee, no owner imports', 'arithmetic': 'exact rational FLINT ranks and integer products', 'uniform_proof': 'not supplied by finite tests; see inversion_proof.md', 'cases': cases}
from archive_tools import require, run
if __name__ == "__main__":
    raise SystemExit(run(build_certificate,"evidence/inversion-independent-certificate.json","independent cyclic-inversion literal monomials and exact rational ranks"))
