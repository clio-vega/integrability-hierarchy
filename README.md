# Integrability Hierarchy

Research on the categorical structure underlying Littlewood-Richardson coefficients,
integrable lattice models, and the representation theory of symmetric groups.

**Author:** Clio (mathematical research AI, collaborating with Robin Langer)

## Contents

### Expository Paper (top level)

A 38pp expository paper tracing LR coefficients from classical combinatorics
through integrable lattice models to the categorical phase diagram:

- `main.tex` / `main.pdf` -- full compiled paper
- `section1-introduction.tex` through `section7-open.tex` -- individual sections
- `phase-diagram.tex` -- standalone categorical phase diagram figure
- `computations/` -- supporting computation scripts

### Proof Documents (`proofs/`)

18 standalone proof documents produced during April 2026 prove sessions.
These develop two main threads:

1. **The categorical phase diagram:** braided -> coboundary -> symmetric monoidal
   categories, cactus group representations, and the multiplicity bundle theorem.

2. **The H-invariant theorem:** rank of the staircase symmetrizing product on
   Specht modules equals dim(V^H), with the Hecke rank constancy generalization.

See `proofs/README.md` for a guide to each document and the logical dependencies.
