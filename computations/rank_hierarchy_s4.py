"""
Rank hierarchy in H_q(S_4) via the left regular representation.

H_q(S_4) has dimension 24 (one basis element T_w per permutation w in S_4).
Generators T_1, T_2, T_3 satisfy:
  T_i^2 = (q-1)T_i + q
  braid relations
  far-commutativity

We represent each T_i as a 24x24 matrix acting on the left regular module.
Left action: T_i * T_w = T_{s_i w}  if s_i w > w  (ell goes up)
             T_i * T_w = (q-1)T_w + q T_{s_i w}    if s_i w < w  (ell goes down)

This is the standard Kazhdan-Lusztig convention (length-compatible left action).

Basis: all permutations of S_4, ordered by length then lex on one-line notation.
"""

import numpy as np
from fractions import Fraction
import itertools

# ---------------------------------------------------------------------------
# S_4 permutation utilities
# ---------------------------------------------------------------------------

def all_perms():
    """All 24 permutations of {0,1,2,3} as tuples."""
    return list(itertools.permutations(range(4)))

def length(w):
    """Number of inversions."""
    w = list(w)
    inv = 0
    for i in range(len(w)):
        for j in range(i+1, len(w)):
            if w[i] > w[j]:
                inv += 1
    return inv

def apply_si(i, w):
    """Apply simple transposition s_i (swap positions i-1 and i in one-line notation).
    s_1 swaps positions 0,1; s_2 swaps 1,2; s_3 swaps 2,3.
    """
    w = list(w)
    w[i-1], w[i] = w[i], w[i-1]
    return tuple(w)

# Build ordered basis: sort by length, then lex
PERMS = sorted(all_perms(), key=lambda w: (length(w), w))
IDX = {w: i for i, w in enumerate(PERMS)}
N = 24

print("Basis (index, permutation, length):")
for i, w in enumerate(PERMS):
    print(f"  {i:2d}: {w}  len={length(w)}")

# ---------------------------------------------------------------------------
# Build T_i matrices over Fraction (exact arithmetic)
# ---------------------------------------------------------------------------

def build_Ti_frac(i, q_val):
    """
    Build the 24x24 matrix of left multiplication by T_i,
    with q = q_val (a Fraction or int).

    Left action on basis element T_w:
      T_i * T_w = T_{s_i w}           if len(s_i w) = len(w) + 1
      T_i * T_w = (q-1)T_w + q T_{s_i w}   if len(s_i w) = len(w) - 1
    """
    q = Fraction(q_val)
    M = [[Fraction(0)]*N for _ in range(N)]
    for col, w in enumerate(PERMS):
        sw = apply_si(i, w)
        row_sw = IDX[sw]
        if length(sw) > length(w):
            # T_i T_w = T_{sw}
            M[row_sw][col] = Fraction(1)
        else:
            # T_i T_w = (q-1) T_w + q T_{sw}
            M[col][col] += (q - 1)
            M[row_sw][col] += q
    return M

def mat_mul_frac(A, B):
    """Matrix multiply two NxN Fraction matrices."""
    C = [[Fraction(0)]*N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            s = Fraction(0)
            for k in range(N):
                if A[i][k] != 0 and B[k][j] != 0:
                    s += A[i][k] * B[k][j]
            C[i][j] = s
    return C

def rank_frac(M):
    """Compute rank of an NxN Fraction matrix via Gaussian elimination."""
    A = [row[:] for row in M]  # copy
    rank = 0
    pivot_row = 0
    for col in range(N):
        # find pivot
        found = -1
        for row in range(pivot_row, N):
            if A[row][col] != 0:
                found = row
                break
        if found == -1:
            continue
        A[pivot_row], A[found] = A[found], A[pivot_row]
        pivot = A[pivot_row][col]
        # normalize pivot row
        for j in range(N):
            A[pivot_row][j] /= pivot
        # eliminate
        for row in range(N):
            if row != pivot_row and A[row][col] != 0:
                factor = A[row][col]
                for j in range(N):
                    A[row][j] -= factor * A[pivot_row][j]
        pivot_row += 1
        rank += 1
    return rank

def identity_frac():
    M = [[Fraction(0)]*N for _ in range(N)]
    for i in range(N):
        M[i][i] = Fraction(1)
    return M

def add_mat(A, B):
    return [[A[i][j] + B[i][j] for j in range(N)] for i in range(N)]

def scalar_mul(c, M):
    c = Fraction(c)
    return [[c * M[i][j] for j in range(N)] for i in range(N)]

def print_rank_results(label, M):
    r = rank_frac(M)
    print(f"  rank({label}) = {r}")
    return r

# ---------------------------------------------------------------------------
# Helper: Ti+1 matrix (i.e., T_i + identity)
# ---------------------------------------------------------------------------

def build_TiPlus1(i, q_val):
    """T_i + 1 matrix."""
    Ti = build_Ti_frac(i, q_val)
    I = identity_frac()
    return add_mat(Ti, I)

# ---------------------------------------------------------------------------
# Main computation for each q value
# ---------------------------------------------------------------------------

def compute_for_q(q_val):
    q = Fraction(q_val)
    print(f"\n{'='*60}")
    print(f"q = {q_val}")
    print(f"{'='*60}")

    T1 = build_Ti_frac(1, q)
    T2 = build_Ti_frac(2, q)
    T3 = build_Ti_frac(3, q)

    # T_{w_0} = T_1 T_2 T_1 T_3 T_2 T_1  (reduced word s1 s2 s1 s3 s2 s1)
    # Apply right to left: first T1, then T2, then T3, then T1, then T2, then T1
    Tw0 = mat_mul_frac(T1, mat_mul_frac(T2, mat_mul_frac(T1, mat_mul_frac(T3, mat_mul_frac(T2, T1)))))

    # Alternative reduced word: T_1 T_2 T_3 T_1 T_2 T_1 = s1s2s3s1s2s1
    # (another reduced expression for w_0 in S_4)
    Tw0_alt = mat_mul_frac(T1, mat_mul_frac(T2, mat_mul_frac(T3, mat_mul_frac(T1, mat_mul_frac(T2, T1)))))

    # Cactus midpoint: Π = (T_1+1)(T_2+1)(T_1+1)(T_3+1)(T_2+1)(T_1+1)
    # = product in same order as T_{w_0}
    Tp1 = build_TiPlus1(1, q)
    Tp2 = build_TiPlus1(2, q)
    Tp3 = build_TiPlus1(3, q)

    Pi = mat_mul_frac(Tp1, mat_mul_frac(Tp2, mat_mul_frac(Tp1, mat_mul_frac(Tp3, mat_mul_frac(Tp2, Tp1)))))

    # q-symmetrizer: sum of all T_w
    # Build T_w for all w via BFS/reduced word
    # We'll compute them iteratively: T_w = T_i * T_{w'} for a reduced factorization w = s_i w'
    # Simpler: store all T_w matrices
    Tw_matrices = {}
    Tw_matrices[PERMS[0]] = identity_frac()  # T_e = identity

    # BFS by length
    for ell in range(1, 7):
        for w in PERMS:
            if length(w) != ell:
                continue
            if w in Tw_matrices:
                continue
            # find a left descent: i such that s_i w < w
            for i in [1,2,3]:
                sw = apply_si(i, w)
                if length(sw) < length(w) and sw in Tw_matrices:
                    Ti = build_Ti_frac(i, q)
                    Tw_matrices[w] = mat_mul_frac(Ti, Tw_matrices[sw])
                    break

    # Sum all T_w
    Sym = [[Fraction(0)]*N for _ in range(N)]
    for w, Mw in Tw_matrices.items():
        Sym = add_mat(Sym, Mw)

    print(f"\nRanks:")
    r_Pi = print_rank_results(f"Pi (cactus midpoint)", Pi)
    r_Tw0 = print_rank_results(f"T_w0 (s1s2s1s3s2s1)", Tw0)
    r_Tw0_alt = print_rank_results(f"T_w0 alt (s1s2s3s1s2s1)", Tw0_alt)
    r_Sym = print_rank_results(f"q-symmetrizer (sum T_w)", Sym)

    # Check if the two reduced words give the same matrix
    diff = add_mat(Tw0, scalar_mul(-1, Tw0_alt))
    diff_rank = rank_frac(diff)
    diff_nonzero = any(diff[i][j] != 0 for i in range(N) for j in range(N))
    print(f"\n  Tw0 == Tw0_alt: {not diff_nonzero}  (diff matrix rank = {diff_rank})")

    # For q=0: also check Pi is rank 1
    if q == 0:
        print(f"\n  [q=0 check] Pi rank={r_Pi} (expected 1 if Pi is a 1-dim projector up to scalar)")
        # Print the first nonzero column of Pi
        for col in range(N):
            col_nonzero = [Pi[row][col] for row in range(N) if Pi[row][col] != 0]
            if col_nonzero:
                print(f"  First nonzero column (col {col}): {col_nonzero[:5]}...")
                break

    return {
        'q': q_val,
        'rank_Pi': r_Pi,
        'rank_Tw0': r_Tw0,
        'rank_Tw0_alt': r_Tw0_alt,
        'rank_Sym': r_Sym,
        'Tw0_equals_Tw0_alt': not diff_nonzero,
    }

# ---------------------------------------------------------------------------
# Run for q = 0, 1/2, 1, 2
# ---------------------------------------------------------------------------

print("="*60)
print("RANK HIERARCHY IN H_q(S_4) — LEFT REGULAR REPRESENTATION")
print("="*60)

results = []
for q_val in [0, Fraction(1,2), 1, 2]:
    r = compute_for_q(q_val)
    results.append(r)

print("\n\n" + "="*60)
print("SUMMARY TABLE")
print("="*60)
print(f"{'q':<8} {'rank(Pi)':<12} {'rank(Tw0)':<12} {'rank(Tw0_alt)':<14} {'rank(Sym)':<12} {'Tw0=Tw0_alt'}")
print("-"*70)
for r in results:
    print(f"{str(r['q']):<8} {r['rank_Pi']:<12} {r['rank_Tw0']:<12} {r['rank_Tw0_alt']:<14} {r['rank_Sym']:<12} {r['Tw0_equals_Tw0_alt']}")

print("\nDone.")
