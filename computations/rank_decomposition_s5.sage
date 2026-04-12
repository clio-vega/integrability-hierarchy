#!/usr/bin/env python3
"""
rank_decomposition_s5.sage
(Implemented in Python since SageMath binary is unavailable; .sage extension per task spec)

For each irreducible representation V_lambda of S_5, compute:
  - dim(V_lambda)
  - rank of the cactus midpoint element product = prod_{i in word} (I + rho(s_i))
    where word = [1,2,1,3,2,1,4,3,2,1] is a reduced word for w_0 in S_5

At q=1 (group algebra) the element becomes prod (e + s_i) in C[S_5].

Method: Young Orthogonal Representation (YOR)
  - Basis: standard Young tableaux for partition lambda
  - Generators: adjacent transpositions s_i = (i,i+1) act via YOR formula

YOR formula for s_i acting on standard tableau T:
  Let d = axial distance between i and i+1 in T (= content(i+1) - content(i))
  where content(box in row r, col c) = c - r (0-indexed)

  If i and i+1 are in the same row: rho(s_i) T = T  (eigenvalue +1)
  If i and i+1 are in the same col: rho(s_i) T = -T  (eigenvalue -1)
  Otherwise: rho(s_i) T = (1/d) T + sqrt(1 - 1/d^2) T'
              rho(s_i) T' = sqrt(1 - 1/d^2) T - (1/d) T'
  where T' is the tableau obtained by swapping i and i+1 in T.

At q=1, the YOR matrices have entries in Q(sqrt(various)), but for our purposes
we use exact rational arithmetic: the eigenvalues 1/d are rational, and the
off-diagonal entries sqrt(1-1/d^2) can be handled exactly.

For rank computations at q=1, we can also use the simpler approach:
  - The matrices for s_i in ANY unitary representation have eigenvalues ±1
  - (I + s_i) has rank = dim(+1 eigenspace of s_i in V_lambda)
  - But we need the rank of the PRODUCT, which is more subtle.

We use the YOR matrices directly with exact rational arithmetic via sympy.
"""

import sys
from itertools import permutations
from fractions import Fraction
from sympy import Rational, Matrix, sqrt, zeros, eye, simplify
from sympy import nsimplify

# ============================================================
# 1. Standard Young Tableaux enumeration
# ============================================================

def partitions_of(n):
    """All partitions of n, in lexicographic descending order."""
    result = []
    def helper(remaining, max_part, current):
        if remaining == 0:
            result.append(tuple(current))
            return
        for p in range(min(remaining, max_part), 0, -1):
            helper(remaining - p, p, current + [p])
    helper(n, n, [])
    return result

def standard_young_tableaux(lam):
    """
    Enumerate all standard Young tableaux of shape lam.
    Returns list of tableaux, each represented as a dict: value -> (row, col) 0-indexed.
    """
    n = sum(lam)
    shape = []
    for r, row_len in enumerate(lam):
        for c in range(row_len):
            shape.append((r, c))

    # Build incrementally: fill boxes of shape with 1..n in increasing order
    # ensuring rows and columns are increasing
    result = []

    def backtrack(pos, tableau, row_counts, col_counts):
        """Fill position pos (= value pos+1) into the Young diagram."""
        if pos == n:
            result.append(dict(tableau))
            return

        val = pos + 1  # value to place (1-indexed)
        # A box (r,c) is "addable" at this step if:
        # - row_counts[r] == c  (next free position in row r)
        # - c == 0 or col_counts[c-1] > col_counts[c]  (no box to the left that is unfilled... wait)
        # Actually: place val at (r,c) if:
        # - It's the leftmost unfilled box in row r
        # - The box above (r-1, c) is already filled (or r=0)

        # Outer corners of current partial tableau:
        # A box (r,c) can be added if row_counts[r] == c and (r==0 or col_counts[c] >= r)
        # where col_counts[c] = number of entries placed in column c so far

        num_cols = lam[0]
        for r in range(len(lam)):
            c = row_counts[r]
            if c >= lam[r]:
                continue  # row r is full
            # Check: box above must be filled
            if r > 0 and row_counts[r-1] <= c:
                continue  # box (r-1, c) not yet filled
            # Place val here
            tableau[val] = (r, c)
            row_counts[r] += 1
            old_cc = col_counts.get(c, 0)
            col_counts[c] = old_cc + 1
            backtrack(pos + 1, tableau, row_counts, col_counts)
            del tableau[val]
            row_counts[r] -= 1
            col_counts[c] = old_cc

    row_counts = [0] * len(lam)
    col_counts = {}
    backtrack(0, {}, row_counts, col_counts)
    return result

# ============================================================
# 2. Young Orthogonal Representation matrices
# ============================================================

def content(r, c):
    """Content of box at (row r, col c), 0-indexed."""
    return c - r

def axial_distance(tableau, i, j):
    """
    Axial distance from i to j in tableau: content(j) - content(i).
    i, j are values (1-indexed).
    """
    ri, ci = tableau[i]
    rj, cj = tableau[j]
    return content(rj, cj) - content(ri, ci)

def yor_matrix(lam, i, SYT):
    """
    Build the YOR matrix for simple transposition s_i = (i, i+1)
    in the representation V_lam, with basis = SYT (standard Young tableaux).

    Returns a sympy Matrix with exact entries.

    YOR formula:
      For basis tableau T (indexed by its position in SYT list):
        d = axial distance from i to i+1 in T
        If i, i+1 in same row (d = 1... wait, same row means c_{i+1} = c_i + something):

    Actually the correct YOR formula:
      d = content(i+1) - content(i) in T

      If i+1 is directly to the right of i (same row, adjacent):
        Actually we use d = axial_distance(T, i, i+1)

      Same row: i and i+1 in same row => s_i T = T  (row reading, but not adjacent boxes)
      Wait, standard YOR: the eigenvalue depends on d.

      Correct YOR:
        If T' = T with i and i+1 swapped is NOT a standard tableau:
          - If i+1 is directly right of i (d = 1): eigenvalue = +1
          - If i+1 is directly below i (d = -1): eigenvalue = -1
        If T' IS standard:
          Matrix entry: rho(s_i)_{T,T} = 1/d
                        rho(s_i)_{T',T} = sqrt(1 - 1/d^2)  [off-diagonal]
                        rho(s_i)_{T,T'} = sqrt(1 - 1/d^2)
                        rho(s_i)_{T',T'} = -1/d
    """
    d_SYT = len(SYT)
    T_index = {id(T): idx for idx, T in enumerate(SYT)}
    # Use list index instead
    T_index = {idx: SYT[idx] for idx in range(d_SYT)}

    # Build index by frozenset content for lookup
    def tableau_key(T):
        return tuple(sorted(T.items()))

    tab_to_idx = {tableau_key(T): idx for idx, T in enumerate(SYT)}

    M = zeros(d_SYT, d_SYT)

    for idx, T in enumerate(SYT):
        d = axial_distance(T, i, i+1)

        # Try to form T' by swapping i and i+1
        T_prime = dict(T)
        T_prime[i], T_prime[i+1] = T[i+1], T[i]

        # Check if T' is standard (rows and columns increasing)
        tp_key = tableau_key(T_prime)

        if tp_key in tab_to_idx:
            # T' is also a standard tableau
            idx_p = tab_to_idx[tp_key]
            d_rat = Rational(1, d)
            # off-diagonal: sqrt(1 - 1/d^2)
            off = sqrt(1 - Rational(1, d*d))
            M[idx, idx] += d_rat          # diagonal block: 1/d
            M[idx_p, idx] += off          # off-diagonal: T' row, T col
            # Note: sympy Matrix is M[row, col]
            # rho(s_i)|_{T,T} = 1/d
            # rho(s_i)|_{T',T} = sqrt(1-1/d^2)  (image of T has component in T')
            # These will be set when we process T as well as T'
            # To avoid double-setting, we set the full 2x2 block when idx < idx_p
            # Reset and do it properly:
        else:
            # T' is not standard; i and i+1 are in the same row or column
            # d = 1 => same row (i+1 immediately right): eigenvalue +1
            # d = -1 => same column (i+1 immediately below): eigenvalue -1
            M[idx, idx] += Rational(d, abs(d))  # +1 or -1

    # Redo more carefully to handle the 2x2 blocks
    M = zeros(d_SYT, d_SYT)
    processed = set()

    for idx, T in enumerate(SYT):
        if idx in processed:
            continue
        d = axial_distance(T, i, i+1)

        T_prime = dict(T)
        T_prime[i], T_prime[i+1] = T[i+1], T[i]
        tp_key = tableau_key(T_prime)

        if tp_key in tab_to_idx:
            idx_p = tab_to_idx[tp_key]
            d_rat = Rational(1, d)
            off = sqrt(1 - Rational(1, d*d))
            # 2x2 block in rows/cols {idx, idx_p}
            # rho(s_i) in this 2x2 basis {T, T'}: [[1/d, off], [off, -1/d]]
            M[idx, idx] = d_rat
            M[idx_p, idx_p] = -d_rat
            M[idx, idx_p] = off
            M[idx_p, idx] = off
            processed.add(idx)
            processed.add(idx_p)
        else:
            # eigenvalue = sign(d) = ±1
            M[idx, idx] = Rational(1, 1) if d > 0 else Rational(-1, 1)
            processed.add(idx)

    return M

# ============================================================
# 3. Main computation
# ============================================================

def compute_rank_for_partition(lam, word):
    """
    Compute rank of prod_{i in word} (I + rho(s_i)) in V_lam.
    """
    SYT = standard_young_tableaux(lam)
    d = len(SYT)

    if d == 0:
        return 0, 0

    # Build matrices for generators s_1, s_2, s_3, s_4
    n = sum(lam)
    rho = {}
    for i in range(1, n):
        rho[i] = yor_matrix(lam, i, SYT)

    # Compute product of (I + rho(s_i)) over word
    I = eye(d)
    product = eye(d)
    for i in word:
        product = product * (I + rho[i])

    # Compute rank over rationals (simplify first if needed)
    # Since entries may involve sqrt, use numerical rank
    # Convert to float for rank computation
    product_float = [[complex(product[r, c]) for c in range(d)] for r in range(d)]

    import numpy as np
    M_np = np.array(product_float, dtype=complex)
    rank = np.linalg.matrix_rank(M_np, tol=1e-8)

    return d, rank

def main():
    n = 5
    parts = partitions_of(n)
    word = [1, 2, 1, 3, 2, 1, 4, 3, 2, 1]  # reduced word for w_0 in S_5

    print("rank_decomposition_s5")
    print("=" * 60)
    print(f"Group: S_{n}")
    print(f"Element: prod_{{i in word}} (I + s_i)")
    print(f"Word (reduced word for w_0): {word}")
    print(f"At q=1 (group algebra level)")
    print()
    print(f"{'Lambda':>15}  {'dim':>5}  {'rank':>6}  {'rank/dim':>10}")
    print("-" * 50)

    total_weighted_rank = 0
    total_dim_sq = 0
    results = []

    for lam in parts:
        dim, rank = compute_rank_for_partition(lam, word)
        results.append((lam, dim, rank))
        total_weighted_rank += rank * dim
        total_dim_sq += dim * dim
        ratio = f"{rank}/{dim}" if dim > 0 else "—"
        print(f"{str(lam):>15}  {dim:>5}  {rank:>6}  {ratio:>10}")

    print("-" * 50)
    print(f"\nSum of dim(V_lambda)^2 = {total_dim_sq}  (should be |S_5| = 120)")
    print(f"Sum of rank(V_lambda) * dim(V_lambda) = {total_weighted_rank}")
    print()
    print("Expected: total weighted rank = ell(w_0) = 10")
    print("(This would mean the product has rank 10 in the regular rep,")
    print(" matching the conjecture from rank_hierarchy_s5.py)")
    print()

    # Also show the regular rep rank (should match rank_hierarchy_s5.py at q=1)
    # By Maschke / Artin-Wedderburn: reg rep = sum_lambda dim(lambda) * V_lambda
    # rank in reg rep = sum_lambda dim(lambda) * rank(lambda)
    reg_rank = sum(dim * rank for (lam, dim, rank) in results)
    print(f"Rank in regular representation = sum_lambda dim * rank = {reg_rank}")
    verified = "VERIFIED" if reg_rank == 30 else "MISMATCH (expected 30)"
    print(f"  -> {verified} (cross-check against rank_hierarchy_s5.py at q=1)")
    print()
    print("Partition details:")
    for lam, dim, rank in results:
        print(f"  {lam}: dim={dim}, rank={rank}")

    # Dimension check
    dim_sq_sum = sum(dim*dim for (_, dim, _) in results)
    print()
    print(f"Sum dim^2 = {dim_sq_sum} {'= 5! = 120 OK' if dim_sq_sum == 120 else '!= 120 ERROR'}")

if __name__ == "__main__":
    main()
