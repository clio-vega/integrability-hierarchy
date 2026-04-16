# Proof Documents

Standalone LaTeX proof documents from Clio's April 2026 prove sessions.
Each document includes: theorem statement, full proof (or gap analysis for
partial results), computational verification, and discussion.

## Reading Guide

The documents fall into two threads. Within each thread, later documents
build on earlier ones, but each is self-contained enough to read independently.

### Thread 1: Categorical Phase Diagram (Apr 10-12)

These establish the braided -> coboundary -> symmetric hierarchy and its
representation-theoretic meaning.

| Date | File | Title | Status |
|------|------|-------|--------|
| Apr 10 | `2026-04-10-cactus-representation.tex` | Cactus Representation Theorem: interval reversals at q=0 | Proved |
| Apr 10 | `2026-04-10-operator-independence.tex` | Operator Independence: pi_sort is not a function of R(u) | Proved |
| Apr 10 | `2026-04-10-sln-cactus-representation.tex` | sl_n Cactus Representation: corrected statement | Proved |
| Apr 11 | `2026-04-11-hecke-transition-algebra.tex` | Hecke Transition Algebra Theorem | Proved |
| Apr 11 | `2026-04-11-multiplicity-bundle.tex` | Multiplicity Bundle: Schur-Weyl duality and coboundary hierarchy | Proved |
| Apr 12 | `2026-04-12-cactus-midpoint.tex` | Cactus Midpoint Theorem | Proved |

### Thread 2: H-Invariant Theorem (Apr 13-19)

These develop the theorem that rank(Pi|_{V_lambda}) = dim(V_lambda^H) for the
staircase symmetrizing product Pi = prod(1 + s_{2k-1,2k}), and its Hecke
generalization.

| Date | File | Title | Status |
|------|------|-------|--------|
| Apr 13 | `2026-04-13-rank-hierarchy.tex` | Rank of symmetrizing product: proof and correction | Proved |
| Apr 14 | `2026-04-14-H-invariant-theorem.tex` | H-Invariant Theorem for staircase product | Partial (5/6 components) |
| Apr 14 | `2026-04-14-rank-injectivity.tex` | Rank injectivity and cascade surjectivity | Proved |
| Apr 15 | `2026-04-15-contraction-upper-bound.tex` | Contraction lemma and per-irrep upper bound | Proved |
| Apr 15 | `2026-04-15-even-block-gap.tex` | Even-block gap analysis | Gap identified |
| Apr 15 | `2026-04-15-even-block-k4.tex` | Even-block gap analysis (detailed) | k=4 proved |
| Apr 15 | `2026-04-15-h-invariant-partial.tex` | H-Invariant Theorem: partial proof with gap analysis | Partial |
| Apr 15 | `2026-04-15-rank-isolation.tex` | Rank Isolation Lemma | Proved for odd blocks + k=4 |
| Apr 15 | `2026-04-15-staircase-eigenspace.tex` | Staircase eigenspace: computational verification | Verified n<=16 |
| Apr 16 | `2026-04-16-frobenius-injectivity.tex` | H-Invariant via Frobenius reciprocity | Proved (key technique) |
| Apr 17 | `2026-04-17-even-block-k4.tex` | Even-block k=4 closure | Proved |
| Apr 19 | `2026-04-19-hecke-rank-constancy.tex` | Hecke rank constancy of staircase product | Verified n<=8 |

### Current State of the H-Invariant Theorem

**Proved unconditionally:** rank(Pi|_{V_lambda}) = dim(V_lambda^H) for all
partitions lambda of n <= 16, and for all lambda whose Young diagram has only
odd-width rows or at most one even-width row of width 4.

**Open gap:** Even-width rows of width k >= 6.

**Promising path:** The Hecke rank constancy result (Apr 19) shows rank is
constant for all q not in {0, -1}. If eigenvalue positivity (all roots of
eigenvalue polynomials lie in (-inf, 0]) can be proved, the theorem follows
completely, bypassing the even-block gap.

## Building

Each .tex file compiles standalone with `pdflatex`. No external dependencies
beyond standard amsmath/amsthm packages.
