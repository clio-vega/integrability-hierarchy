"""
three_ybe_test.py

Investigation of three YBE structures in the Hecke algebra H_q(S_3).

H_q(S_3) has basis {1, T1, T2, T1T2, T2T1, T1T2T1}
with Hecke relation T_i^2 = (q-1)T_i + q.

Cactus midpoint theorem: R_i(1/2) = (T_i+1)/2 at generic q.
At q=0: R_i(0) = (T_i+1)/2 becomes the Demazure/cactus operator.

We study ∏ R_{i_j}(1/2) for the reduced word w_0 = s1 s2 s1 in S_3.

Author: Clio, 2026-04-12
"""

import sympy as sp
from sympy import Matrix, symbols, Rational, factor, simplify, expand, zeros, eye

print("=" * 70)
print("THREE YBE STRUCTURES IN H_q(S_3): CACTUS MIDPOINT ANALYSIS")
print("=" * 70)

# -----------------------------------------------------------------------
# SECTION 1: Build the left regular representation of H_q(S_3)
# -----------------------------------------------------------------------
#
# Basis ordering: 0=1, 1=T1, 2=T2, 3=T1T2, 4=T2T1, 5=T1T2T1
# We label them:  e, T1, T2, T12, T21, T121
#
# Hecke relation: T_i^2 = (q-1)*T_i + q*1
# Braid relation: T1*T2*T1 = T2*T1*T2
#
# We'll build the left multiplication matrices for T1 and T2.
# M_Ti[j,k] = coefficient of basis[j] in Ti * basis[k].

q = symbols('q')

# Name the basis elements by index
BASIS = ['1', 'T1', 'T2', 'T1T2', 'T2T1', 'T1T2T1']
N = 6

def basis_vec(i):
    """Standard basis vector e_i as a sympy column vector."""
    v = zeros(N, 1)
    v[i] = 1
    return v

# Left multiplication by T1:
#   T1 * 1       = T1            (index 0 -> index 1)
#   T1 * T1      = (q-1)T1 + q  (using T1^2 = (q-1)T1 + q)
#   T1 * T2      = T1T2          (index 2 -> index 3)
#   T1 * T1T2    = (q-1)T1T2 + q*T2   (T1^2 * T2 = ((q-1)T1+q)*T2)
#   T1 * T2T1    = T1T2T1        (index 4 -> index 5)
#   T1 * T1T2T1  = (q-1)T1T2T1 + q*T2T1  (T1^2*T2T1)

LT1 = zeros(N, N)
# Column 0: T1*1 = T1
LT1[1, 0] = 1
# Column 1: T1*T1 = (q-1)T1 + q*1
LT1[0, 1] = q
LT1[1, 1] = q - 1
# Column 2: T1*T2 = T1T2
LT1[3, 2] = 1
# Column 3: T1*(T1T2) = (q-1)T1T2 + q*T2
LT1[2, 3] = q
LT1[3, 3] = q - 1
# Column 4: T1*(T2T1) = T1T2T1
LT1[5, 4] = 1
# Column 5: T1*(T1T2T1) = (q-1)T1T2T1 + q*T2T1
LT1[4, 5] = q
LT1[5, 5] = q - 1

print("\n--- Left multiplication matrix for T1 ---")
sp.pprint(LT1)

# Left multiplication by T2:
#   T2 * 1       = T2
#   T2 * T1      = T2T1
#   T2 * T2      = (q-1)T2 + q
#   T2 * T1T2    = T1T2T1    (because T2*(T1T2) = T1*(T2T2) ?? No!)
#
# Wait — we need to be careful. T2*(T1T2):
#   This is T2*T1*T2. By braid relation, T1T2T1 = T2T1T2, so T2T1T2 = T1T2T1.
#   Hence T2*(T1*T2) = (T2T1)*T2 = T2T1T2 = T1T2T1. Yes, index 5.
#
#   T2 * T2T1    = T2^2 * T1 = ((q-1)T2 + q)*T1 = (q-1)T2T1 + q*T1
#     Wait: T2*(T2T1) = (T2^2)*T1 = ((q-1)T2 + q)*T1 = (q-1)*T2T1 + q*T1
#
#   T2 * T1T2T1  = T2*T1*T2*T1.
#     Using braid: T2T1T2 = T1T2T1, so T2T1T2T1 = T1T2T1T1 = T1T2*((q-1)T1+q)
#     = (q-1)*T1T2T1 + q*T1T2
#     Alternatively: T2*(T1T2T1) = (T2T1T2)*T1 = (T1T2T1)*T1 = T1T2*(T1^2)
#     = T1T2*((q-1)T1 + q) = (q-1)*T1T2T1 + q*T1T2

LT2 = zeros(N, N)
# Column 0: T2*1 = T2
LT2[2, 0] = 1
# Column 1: T2*T1 = T2T1
LT2[4, 1] = 1
# Column 2: T2*T2 = (q-1)T2 + q
LT2[0, 2] = q
LT2[2, 2] = q - 1
# Column 3: T2*(T1T2) = T1T2T1  [via braid: T2T1T2 = T1T2T1]
LT2[5, 3] = 1
# Column 4: T2*(T2T1) = (q-1)T2T1 + q*T1
LT2[1, 4] = q
LT2[4, 4] = q - 1
# Column 5: T2*(T1T2T1) = (q-1)T1T2T1 + q*T1T2
LT2[3, 5] = q
LT2[5, 5] = q - 1

print("\n--- Left multiplication matrix for T2 ---")
sp.pprint(LT2)

# -----------------------------------------------------------------------
# SECTION 2: Verify the Hecke and braid relations
# -----------------------------------------------------------------------

print("\n--- Verification ---")

# Check T1^2 = (q-1)T1 + q*I
T1sq = LT1 * LT1
expected_T1sq = (q - 1) * LT1 + q * eye(N)
diff = sp.simplify(T1sq - expected_T1sq)
print(f"T1^2 = (q-1)T1 + q: {'PASS' if diff == zeros(N, N) else 'FAIL'}")

# Check T2^2 = (q-1)T2 + q*I
T2sq = LT2 * LT2
expected_T2sq = (q - 1) * LT2 + q * eye(N)
diff = sp.simplify(T2sq - expected_T2sq)
print(f"T2^2 = (q-1)T2 + q: {'PASS' if diff == zeros(N, N) else 'FAIL'}")

# Check braid: T1 T2 T1 = T2 T1 T2
T1T2T1 = LT1 * LT2 * LT1
T2T1T2 = LT2 * LT1 * LT2
diff = sp.simplify(T1T2T1 - T2T1T2)
print(f"Braid T1T2T1 = T2T1T2: {'PASS' if diff == zeros(N, N) else 'FAIL'}")

# -----------------------------------------------------------------------
# SECTION 3: The cactus midpoint operator R_i(1/2) = (T_i + 1) / 2
# -----------------------------------------------------------------------
# At generic q: R_i(1/2) = (T_i + 1) / 2
# At q=0: same formula but T_i^2 = -T_i so T_i is nilpotent-like

print("\n" + "=" * 70)
print("SECTION 3: CACTUS MIDPOINT OPERATOR")
print("=" * 70)

# R_i = (T_i + I) / 2
R1 = (LT1 + eye(N)) / 2
R2 = (LT2 + eye(N)) / 2

print("\nR1 = (T1 + I)/2:")
sp.pprint(R1)

print("\nR2 = (T2 + I)/2:")
sp.pprint(R2)

# -----------------------------------------------------------------------
# SECTION 4: Product for w_0 = s1 s2 s1
# -----------------------------------------------------------------------
# Π = R1 * R2 * R1 = ((T1+1)/2) * ((T2+1)/2) * ((T1+1)/2)

print("\n" + "=" * 70)
print("SECTION 4: PRODUCT Π = R1·R2·R1 FOR w_0 = s1s2s1")
print("=" * 70)

Pi = R1 * R2 * R1
Pi_simplified = sp.simplify(Pi)

print("\nΠ = R1·R2·R1 (simplified):")
sp.pprint(Pi_simplified)

# Also compute without the 1/2 factors: (T1+1)(T2+1)(T1+1)
Pi_unnorm = (LT1 + eye(N)) * (LT2 + eye(N)) * (LT1 + eye(N))
Pi_unnorm_simplified = sp.simplify(Pi_unnorm)

print("\n(T1+1)(T2+1)(T1+1) [unnormalized]:")
sp.pprint(Pi_unnorm_simplified)

# -----------------------------------------------------------------------
# SECTION 5: Express in terms of basis coefficients
# -----------------------------------------------------------------------
# The matrix represents left-multiplication by some element of H_q(S_3).
# The element can be read off by applying the matrix to the identity element (basis vector 0).
# [Because M * e_0 gives the coefficients of M*1 in the basis.]

print("\n" + "=" * 70)
print("SECTION 5: ELEMENT IN THE T-BASIS")
print("=" * 70)

e0 = basis_vec(0)

Pi_element = Pi_simplified * e0
Pi_element_simplified = sp.simplify(Pi_element)

print("\nΠ·1 (element in H_q(S3), expressed in T-basis):")
for i, name in enumerate(BASIS):
    coeff = Pi_element_simplified[i]
    if coeff != 0:
        print(f"  [{name}]: {sp.factor(coeff)}")

Pi_unnorm_element = Pi_unnorm_simplified * e0
Pi_unnorm_element_simplified = sp.simplify(Pi_unnorm_element)

print("\n(T1+1)(T2+1)(T1+1)·1 [unnormalized, in T-basis]:")
for i, name in enumerate(BASIS):
    coeff = Pi_unnorm_element_simplified[i]
    if coeff != 0:
        print(f"  [{name}]: {sp.factor(coeff)}")

# -----------------------------------------------------------------------
# SECTION 6: Rank analysis as function of q
# -----------------------------------------------------------------------

print("\n" + "=" * 70)
print("SECTION 6: RANK ANALYSIS")
print("=" * 70)

# Rank of Π at symbolic q
print("\nRank of Π = R1·R2·R1 at symbolic q:")
print("(Computing rank via polynomial row reduction...)")

# Try specific values first to understand the pattern
test_values = [sp.Integer(0), sp.Rational(1,2), sp.Integer(1), sp.Integer(-1), sp.Integer(2)]
print("\nRank at specific q values:")
for qval in test_values:
    M_num = Pi_simplified.subs(q, qval)
    r = M_num.rank()
    print(f"  q = {qval}: rank = {r}")

# For unnormalized version
print("\nRank of (T1+1)(T2+1)(T1+1) at specific q values:")
for qval in test_values:
    M_num = Pi_unnorm_simplified.subs(q, qval)
    r = M_num.rank()
    print(f"  q = {qval}: rank = {r}")

# -----------------------------------------------------------------------
# SECTION 7: Minimal polynomial of the unnormalized operator
# -----------------------------------------------------------------------

print("\n" + "=" * 70)
print("SECTION 7: EIGENSTRUCTURE")
print("=" * 70)

print("\nEigenvalues of Π = R1·R2·R1 at specific q values:")
for qval in [sp.Integer(0), sp.Rational(1, 4), sp.Rational(1, 2), sp.Integer(1)]:
    M_num = Pi_simplified.subs(q, qval)
    try:
        eigs = M_num.eigenvals()
        print(f"  q = {qval}: eigenvalues = {dict(eigs)}")
    except Exception as e:
        print(f"  q = {qval}: eigenvalue computation failed: {e}")

print("\nEigenvalues of (T1+1)(T2+1)(T1+1) at specific q values:")
for qval in [sp.Integer(0), sp.Rational(1, 4), sp.Rational(1, 2), sp.Integer(1)]:
    M_num = Pi_unnorm_simplified.subs(q, qval)
    try:
        eigs = M_num.eigenvals()
        print(f"  q = {qval}: eigenvalues = {dict(eigs)}")
    except Exception as e:
        print(f"  q = {qval}: eigenvalue computation failed: {e}")

# -----------------------------------------------------------------------
# SECTION 8: Minimal polynomial of unnormalized at generic q
# -----------------------------------------------------------------------

print("\n" + "=" * 70)
print("SECTION 8: MINIMAL POLYNOMIAL (unnormalized)")
print("=" * 70)

# Characteristic polynomial of unnormalized
x = symbols('x')
print("\nCharacteristic polynomial of (T1+1)(T2+1)(T1+1) at q=0:")
M0 = Pi_unnorm_simplified.subs(q, 0)
char_poly_0 = M0.charpoly(x)
print(sp.factor(char_poly_0.as_expr()))

print("\nCharacteristic polynomial of (T1+1)(T2+1)(T1+1) at q=1:")
M1 = Pi_unnorm_simplified.subs(q, 1)
char_poly_1 = M1.charpoly(x)
print(sp.factor(char_poly_1.as_expr()))

# -----------------------------------------------------------------------
# SECTION 9: Is Π a projector? Check Π^2 vs Π
# -----------------------------------------------------------------------

print("\n" + "=" * 70)
print("SECTION 9: PROJECTOR CHECK Π^2 = c·Π?")
print("=" * 70)

Pi_sq = Pi_simplified * Pi_simplified
Pi_sq_simplified = sp.simplify(Pi_sq)

# Check if Pi_sq = lambda * Pi for some lambda
# Divide: if Pi is rank 1, then Pi_sq = c * Pi where c = Tr(Pi)
# (For a rank-1 matrix A = u*v^T, A^2 = (v^T u) A)

print("\nΠ^2 - c·Π = 0 check (looking for scalar c):")
# Try to find c such that Pi_sq = c * Pi
# Use the (0,0) entry ratio if nonzero
for qval in [sp.Integer(0), sp.Rational(1, 2), sp.Integer(1)]:
    Pi_num = Pi_simplified.subs(q, qval)
    Pi_sq_num = Pi_sq_simplified.subs(q, qval)

    # Find ratio where Pi is nonzero
    ratio = None
    for i in range(N):
        for j in range(N):
            if Pi_num[i, j] != 0:
                ratio = Pi_sq_num[i, j] / Pi_num[i, j]
                break
        if ratio is not None:
            break

    if ratio is not None:
        residual = sp.simplify(Pi_sq_num - ratio * Pi_num)
        is_projector_scaled = (residual == zeros(N, N))
        print(f"  q = {qval}: Π^2 = {ratio}·Π? {is_projector_scaled}")
    else:
        print(f"  q = {qval}: Π = 0")

# -----------------------------------------------------------------------
# SECTION 10: The Demazure operator at q=0
# -----------------------------------------------------------------------

print("\n" + "=" * 70)
print("SECTION 10: DEMAZURE/CACTUS OPERATOR AT q=0")
print("=" * 70)

LT1_0 = LT1.subs(q, 0)
LT2_0 = LT2.subs(q, 0)

print("\nT1 at q=0:")
sp.pprint(LT1_0)

print("\nT2 at q=0:")
sp.pprint(LT2_0)

# At q=0: T_i^2 = -T_i, so T_i(T_i + 1) = 0: T_i is a projector (up to sign)
print("\nCheck T1^2 + T1 = 0 at q=0:")
check = sp.simplify(LT1_0 * LT1_0 + LT1_0)
print(f"  PASS" if check == zeros(N, N) else f"  Result: {check}")

# The cactus operator at q=0
D1 = LT1_0 + eye(N)   # T1 + 1
D2 = LT2_0 + eye(N)   # T2 + 1
Demazure = D1 * D2 * D1

print("\nDemazure operator (T1+1)(T2+1)(T1+1) at q=0:")
sp.pprint(Demazure)

print("\nRank at q=0:", Demazure.rank())
print("Trace at q=0:", Demazure.trace())

# What element of H_q(S_3) does Demazure give?
Dem_element = Demazure * e0
print("\nDemazure·1 in T-basis:")
for i, name in enumerate(BASIS):
    coeff = Dem_element[i]
    if coeff != 0:
        print(f"  [{name}]: {coeff}")

# -----------------------------------------------------------------------
# SECTION 11: Kamnitzer-Tingley commutor
# -----------------------------------------------------------------------

print("\n" + "=" * 70)
print("SECTION 11: KAMNITZER-TINGLEY COMMUTOR")
print("=" * 70)

# KT commutor = q^{-ℓ(w_0)/2} * T_{w_0}
# For S_3: ℓ(w_0) = 3, T_{w_0} = T1T2T1
# KT_{w_0} = q^{-3/2} * T1T2T1
#
# In the regular representation, T_{w_0} is basis element 5.
# Its left multiplication matrix:
LT121 = LT1 * LT2 * LT1  # = T1T2T1 left multiplication

print("\nT1T2T1 left multiplication matrix:")
sp.pprint(sp.simplify(LT121))

# The KT operator (up to the q^{-3/2} scalar, which doesn't affect rank):
print("\nRank of T_{w_0} matrix at specific q values:")
for qval in test_values:
    M_num = LT121.subs(q, qval)
    r = M_num.rank()
    print(f"  q = {qval}: rank = {r}")

# -----------------------------------------------------------------------
# SECTION 12: Comparison and summary
# -----------------------------------------------------------------------

print("\n" + "=" * 70)
print("SECTION 12: COMPARISON SUMMARY")
print("=" * 70)

print("\nOperator | q=0 rank | q=1/2 rank | q=1 rank | q=2 rank")
print("-" * 60)

operators = {
    "R1·R2·R1 (normalized)": Pi_simplified,
    "(T1+1)(T2+1)(T1+1)": Pi_unnorm_simplified,
    "T_{w_0} = T1T2T1": sp.simplify(LT121),
}

qvals_for_table = [sp.Integer(0), sp.Rational(1, 2), sp.Integer(1), sp.Integer(2)]

for name, M in operators.items():
    ranks = []
    for qval in qvals_for_table:
        M_num = M.subs(q, qval)
        ranks.append(M_num.rank())
    print(f"{name[:35]:35s} | {ranks[0]:8} | {ranks[1]:10} | {ranks[2]:8} | {ranks[3]:8}")

# -----------------------------------------------------------------------
# SECTION 13: Detect q values where rank drops
# -----------------------------------------------------------------------

print("\n" + "=" * 70)
print("SECTION 13: RANK-DROP LOCUS FOR (T1+1)(T2+1)(T1+1)")
print("=" * 70)

# Compute the determinant of Pi_unnorm
det_Pi = Pi_unnorm_simplified.det()
print("\nDet[(T1+1)(T2+1)(T1+1)]:")
print(sp.factor(det_Pi))

# Compute rank-deficiency locus via minors
# Matrix is 6x6; if det=0, rank < 6. Find when rank < 2 (i.e., rank=1)
# by checking 2x2 minors
print("\nCheck: does rank=1 only at q=0?")
print("Testing many rational q values:")
rank_1_found = []
rank_other = {}
for num in range(-5, 10):
    for den in [1, 2, 3, 4]:
        qval = sp.Rational(num, den)
        M_num = Pi_unnorm_simplified.subs(q, qval)
        r = M_num.rank()
        if r == 1:
            rank_1_found.append(qval)
        elif r not in rank_other:
            rank_other[r] = qval

print(f"  q values with rank 1: {rank_1_found}")
print(f"  Example q values for other ranks: {rank_other}")

print("\n" + "=" * 70)
print("COMPUTATION COMPLETE")
print("=" * 70)
