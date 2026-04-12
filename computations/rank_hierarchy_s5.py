"""
rank_hierarchy_s5.py

Compute the rank of the cactus midpoint element prod(T_i + I) in the
left regular representation of H_q(S_5) at specific numerical q values.

H_q(S_5): dim = 120, generators T_1,...,T_4 (adjacent transpositions s_i).

Hecke multiplication rules (left action on basis {T_w}):
  T_i * T_w = T_{s_i w}           if ell(s_i w) > ell(w)
  T_i * T_w = (q-1)*T_w + q*T_{s_i w}  if ell(s_i w) < ell(w)

pi_i = T_i + I

Conjecture: rank(prod pi over reduced word of w_0) = ell(w_0) = 10 at generic q,
            rank = 1 at q = 0.
"""

import numpy as np
from itertools import permutations

# -------------------------------------------------------------------
# 1. Enumerate S_5 and build permutation index
# -------------------------------------------------------------------
S5 = list(permutations(range(5)))          # 120 permutations as tuples
perm_index = {p: i for i, p in enumerate(S5)}
N = len(S5)  # 120

def inversion_count(perm):
    """Length (number of inversions) of permutation."""
    count = 0
    n = len(perm)
    for i in range(n):
        for j in range(i+1, n):
            if perm[i] > perm[j]:
                count += 1
    return count

ell = {p: inversion_count(p) for p in S5}

def left_multiply_si(p, i):
    """Apply simple transposition s_i (swap positions i-1, i) on the LEFT of p.
    s_i acts on values: swap the values at positions i-1 and i.
    Actually s_i as a permutation swaps i-1 <-> i (0-indexed: i-1 and i).
    Left multiplication: (s_i * p)(x) = s_i(p(x)).
    As a tuple, s_i * p swaps the two VALUES i-1 and i in the image of p.
    Wait — we use one-line notation. Let's be careful.

    Permutations in one-line notation: p = (p[0], p[1], ..., p[n-1]).
    s_i (1-indexed generator, swap positions i and i+1 in 1-indexed)
    = swap values i-1 and i in 0-indexed notation.

    Left multiply: s_i * p swaps positions (i-1) and i in the tuple
    (since s_i as a function swaps the ELEMENTS at positions i-1 and i).
    """
    lst = list(p)
    # s_i (1-indexed) swaps positions i-1 and i in 0-indexed
    lst[i-1], lst[i] = lst[i], lst[i-1]
    return tuple(lst)

# Verify: s_1 * (0,1,2,3,4) should give (1,0,2,3,4)
assert left_multiply_si((0,1,2,3,4), 1) == (1,0,2,3,4), "s_1 check failed"
# s_2 * (0,1,2,3,4) = (0,2,1,3,4)
assert left_multiply_si((0,1,2,3,4), 2) == (0,2,1,3,4), "s_2 check failed"

# -------------------------------------------------------------------
# 2. Build T_i matrices for a given q
# -------------------------------------------------------------------

def build_Ti(i, q):
    """Build 120x120 matrix of left action of T_i on H_q(S_5)."""
    M = np.zeros((N, N), dtype=np.float64)
    for p in S5:
        col = perm_index[p]           # p is the basis element T_p
        sp = left_multiply_si(p, i)   # s_i * p
        row_sp = perm_index[sp]
        if ell[sp] > ell[p]:
            # T_i * T_p = T_{s_i p}
            M[row_sp, col] += 1.0
        else:
            # T_i * T_p = (q-1)*T_p + q*T_{s_i p}
            M[col, col] += (q - 1.0)
            M[row_sp, col] += q
    return M

# -------------------------------------------------------------------
# 3. Reduced words for w_0 in S_5
#    w_0 = [4,3,2,1,0] (0-indexed), length = 10
#    Word 1: s1 s2 s1 s3 s2 s1 s4 s3 s2 s1  (standard)
#    Word 2: s4 s3 s4 s2 s3 s4 s1 s2 s3 s4  (reverse standard)
# -------------------------------------------------------------------

word1 = [1, 2, 1, 3, 2, 1, 4, 3, 2, 1]
word2 = [4, 3, 4, 2, 3, 4, 1, 2, 3, 4]

# Verify word1 gives w_0
def apply_word(word, start=None):
    """Apply a sequence of simple transpositions to the identity."""
    if start is None:
        p = tuple(range(5))
    else:
        p = start
    for i in word:
        p = left_multiply_si(p, i)
    return p

w0_from_word1 = apply_word(word1)
w0_from_word2 = apply_word(word2)
w0_expected = (4, 3, 2, 1, 0)

print(f"w_0 from word1: {w0_from_word1}  (expected {w0_expected})")
print(f"w_0 from word2: {w0_from_word2}  (expected {w0_expected})")
assert w0_from_word1 == w0_expected, f"word1 does not give w_0: got {w0_from_word1}"
assert w0_from_word2 == w0_expected, f"word2 does not give w_0: got {w0_from_word2}"
print("Both reduced words verified to give w_0 = (4,3,2,1,0).\n")

# -------------------------------------------------------------------
# 4. Main computation
# -------------------------------------------------------------------

q_values = [0.0, 0.5, 1.0, 2.0]

print(f"{'q':>6}  {'rank(word1)':>12}  {'rank(word2)':>12}  {'same?':>6}  {'note'}")
print("-" * 60)

for q in q_values:
    # Build T_i and pi_i = T_i + I for each generator
    I = np.eye(N, dtype=np.float64)
    Pi = {}
    for i in range(1, 5):
        Ti = build_Ti(i, q)
        Pi[i] = Ti + I

    # Product for word1: pi_1 * pi_2 * pi_1 * pi_3 * pi_2 * pi_1 * pi_4 * pi_3 * pi_2 * pi_1
    Prod1 = np.eye(N, dtype=np.float64)
    for i in word1:
        Prod1 = Prod1 @ Pi[i]

    # Product for word2
    Prod2 = np.eye(N, dtype=np.float64)
    for i in word2:
        Prod2 = Prod2 @ Pi[i]

    r1 = np.linalg.matrix_rank(Prod1)
    r2 = np.linalg.matrix_rank(Prod2)
    same = "YES" if r1 == r2 else "NO"

    # Annotate
    if q == 0.0:
        note = "<-- q=0 (coboundary/crystal)"
    elif q == 0.5:
        note = "<-- midpoint q=1/2"
    elif q == 1.0:
        note = "<-- q=1 (group algebra)"
    else:
        note = "<-- generic q"

    print(f"{q:>6.2f}  {r1:>12}  {r2:>12}  {same:>6}  {note}")

print()

# -------------------------------------------------------------------
# 5. Detailed check at q=0: verify rank=1 and what the image is
# -------------------------------------------------------------------
print("=== Detailed check at q=0 ===")
q = 0.0
I = np.eye(N, dtype=np.float64)
Pi0 = {}
for i in range(1, 5):
    Ti = build_Ti(i, q)
    Pi0[i] = Ti + I

Prod0 = np.eye(N, dtype=np.float64)
for i in word1:
    Prod0 = Prod0 @ Pi0[i]

# At q=0, T_i is a projection: T_i^2 = -T_i (since T_i^2=(q-1)T_i+q*I -> -(T_i))
# pi_i = T_i + I is an idempotent (pi_i^2 = pi_i at q=0)
# Check idempotency of pi_1 at q=0
Pi1_sq = Pi0[1] @ Pi0[1]
diff = np.max(np.abs(Pi1_sq - Pi0[1]))
print(f"pi_1 idempotent at q=0? max|pi_1^2 - pi_1| = {diff:.2e}")

# The image should be spanned by sum of T_w over all w
# Check: Prod0 should have all columns proportional
col0 = Prod0[:, 0]
nz = np.where(np.abs(col0) > 1e-10)[0]
print(f"Nonzero entries in column 0 of Prod (q=0): {len(nz)}")
if len(nz) > 0:
    print(f"Values: {np.unique(np.round(col0[nz], 8))}")

# -------------------------------------------------------------------
# 6. Check at q=1: compare to S_5 group algebra
# -------------------------------------------------------------------
print("\n=== Check at q=1 (group algebra) ===")
q = 1.0
I = np.eye(N, dtype=np.float64)
Pi1_mats = {}
for i in range(1, 5):
    Ti = build_Ti(i, q)
    Pi1_mats[i] = Ti + I  # = T_i + I = s_i + e (in group algebra)

# At q=1, T_i = permutation matrix for s_i, so pi_i = s_i + e
# Product is (s1+e)(s2+e)...(s1+e) = sum of all subwords = sum over subsets of word1 of T_{product}
# This is related to the "left-to-right" Gaussian binomial / descent polynomial

Prod1_q1 = np.eye(N, dtype=np.float64)
for i in word1:
    Prod1_q1 = Prod1_q1 @ Pi1_mats[i]

r_q1 = np.linalg.matrix_rank(Prod1_q1)
print(f"Rank at q=1: {r_q1}")
# Also check Frobenius norm
print(f"||Prod(q=1)||_F = {np.linalg.norm(Prod1_q1):.4f}")
# In the group algebra of S_n, this product should be n! * (sum over S_n) / something
# The all-ones-like matrix... let's check if all rows sum to same value
row_sums = Prod1_q1.sum(axis=1)
print(f"Row sum range at q=1: min={row_sums.min():.4f}, max={row_sums.max():.4f}")

print("\nDone.")
