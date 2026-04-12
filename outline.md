# Littlewood-Richardson Coefficients: From Puzzles to the Integrability Hierarchy

**Working title.** Expository paper, target length 15-20 pages.

**Audience:** Researchers familiar with symmetric functions and Young tableaux, but not necessarily with integrable lattice models or quantum groups. The paper should make the integrable perspective accessible while building toward the Yang-Baxter groupoid framework and the integrability hierarchy as organizing principles.

**Thesis:** The multiplicity of combinatorial interpretations of LR coefficients is not accidental -- it is explained by a three-layer structure: combinatorial models (surface), the Yang-Baxter groupoid (structure), and shuffle algebras (algebra). These layers are further organized by a hierarchy of integrability levels -- solvable, quasi-solvable, combinatorial -- corresponding to the Hecke algebra chain $H_q \to H_0 \to S_n$. The hierarchy explains not only which models exist but also which have positive structure constants (solvable/achiral) and which have signed ones (quasi-solvable/chiral).

---

## Section 1. Introduction (2 pages)

**Content:** State the motivating question: why do LR coefficients admit so many independent combinatorial interpretations (tableaux, puzzles, pipe dreams, lattice paths, interlacing triangular arrays, graph colorings)? Each model was discovered independently, often from different mathematical communities. The bijections between them are known but seem ad hoc. We argue that integrability -- specifically the Yang-Baxter equation -- provides the structural explanation.

Preview the paper's arc: classical definitions (Section 2), the integrable lattice model perspective (Section 3), the Yang-Baxter groupoid (Section 4), the integrability hierarchy (Section 5), and open questions (Section 6).

**Type:** Pure exposition.

---

## Section 2. Classical Littlewood-Richardson Coefficients (3-4 pages)

### 2.1 Definitions and first properties

**Content:** Define LR coefficients via the three standard routes: (i) structure constants of the ring of symmetric functions, $s_\mu \cdot s_\nu = \sum_\lambda c^\lambda_{\mu\nu} s_\lambda$; (ii) tensor product multiplicities of $GL_n$ representations; (iii) Schubert calculus on Grassmannians, $\sigma_\mu \cdot \sigma_\nu = \sum c^\lambda_{\mu\nu} \sigma_\lambda$ in $H^*(Gr(k,n))$. State symmetries ($c^\lambda_{\mu\nu} = c^\lambda_{\nu\mu}$, commutativity; conjugation symmetry; the Knutson-Tao saturation theorem).

**Type:** Pure exposition. Standard material, presented concisely.

### 2.2 The Littlewood-Richardson rule and its variants

**Content:** State the LR rule via content-reading of skew tableaux. Then present the KTW puzzle model (Knutson-Tao-Woodward): triangles tiled by three types of unit rhombi, boundary conditions encoded by binary strings from partitions. Each valid tiling witnesses one unit of $c^\lambda_{\mu\nu}$. Mention Purbhoo's mosaics as a bridge between tableaux and puzzles. Include a worked example ($c^{(3,2,1)}_{(2,1),(2,1)} = 2$) showing both a puzzle pair and the corresponding skew tableaux pair.

**Type:** Pure exposition, but include a self-contained figure showing the two puzzles for $c^{(3,2,1)}_{(2,1),(2,1)} = 2$. *Computation needed:* generate and verify puzzle diagrams. Robin's `lr_coefficients.py` can enumerate these.

### 2.3 The Hopf algebra perspective

**Content:** The Fock space / Hopf algebra encoding: the coproduct $\Delta(s_\lambda) = \sum c^\lambda_{\mu\nu}\, s_\mu \otimes s_\nu$ packages all LR coefficients at once. Mention the binary-string encoding of partitions, the operators $T_{\text{free}}$ and $\bar{T}_{\text{free}}$ (removing horizontal/vertical strips), and the Kostka matrix inversion approach. This sets up the transfer matrix perspective of Section 3.

**Type:** Mostly exposition. Reference Robin's `transfer_operators.py` for the computational implementation. Briefly note the discovery that Kostka matrix inversion, not hook-content evaluation, is the correct computational strategy (due to linear dependencies among hook-content polynomials).

---

## Section 3. The Integrable Lattice Model Perspective (3-4 pages)

### 3.1 From puzzles to vertex models

**Content:** Zinn-Justin's (2008) key insight: KTW puzzle pieces ARE R-matrix weights satisfying the Yang-Baxter equation. Reformulate the puzzle as a lattice model: edges carry labels from $\{0, 1, +\}$, each triangle is a vertex with a Boltzmann weight, and a valid tiling is a configuration with nonzero total weight. The transfer matrix processes the lattice row by row -- this is exactly the row transfer matrix in `lr_coefficients.py`. Explain how the deterministic $\triangle$ and branching $\nabla$ triangles correspond to the transfer matrix structure.

**Type:** Exposition of Zinn-Justin (2008), presented in a way accessible to combinatorialists. *Computation needed:* verify the R-matrix / Yang-Baxter equation for the KTW tile weights in a self-contained example. This is a small symbolic calculation.

### 3.2 The Yang-Baxter equation

**Content:** State the YBE: $R_{12} R_{13} R_{23} = R_{23} R_{13} R_{12}$, acting on $V \otimes V \otimes V$. Explain what it means physically (commuting transfer matrices) and combinatorially (a local move on tiling configurations that preserves the boundary). Show that the KTW R-matrix satisfies it. Emphasize: the YBE is not a metaphor -- it is the precise reason why the transfer matrix approach works, and why different row-decompositions of the same lattice yield the same count.

**Type:** Exposition. The YBE verification for KTW tiles is in the literature (Zinn-Justin 2008) but we should include a clean, self-contained presentation.

### 3.3 Transfer matrices and commutativity

**Content:** The transfer matrix $T(u)$ depends on a spectral parameter $u$. The YBE implies $[T(u), T(v)] = 0$ for all $u, v$. This commutativity is the algebraic engine: it means the lattice model can be "solved" in the sense that partition functions (= LR coefficients) can be computed by diagonalizing a single transfer matrix. Connect back to the Hopf coproduct of Section 2.3: the coproduct IS the transfer matrix trace over intermediate states.

**Type:** Exposition, connecting the integrable systems language to the symmetric functions language. *Computation desirable:* demonstrate transfer matrix commutativity for a small example using SageMath or the existing Python code.

---

## Section 4. The Yang-Baxter Groupoid (3-4 pages)

### 4.1 The groupoid framework (Bump-Naprienko)

**Content:** Present the Yang-Baxter groupoid of Bump and Naprienko (2503.05960). Objects are combinatorial models (each defined by a choice of boundary data and local vertex weights). Morphisms are R-matrix transformations: applying the YBE to a pair of adjacent rows changes the "model" without changing the partition function. Two models connected by a morphism compute the same coefficient. The groupoid structure (composition, inverses) comes from iterating R-matrix moves.

State the key theorem: all known combinatorial models for LR coefficients (puzzles, tableaux, pipe dreams, lattice paths) sit at different objects of this groupoid, and the known bijections between them arise as groupoid morphisms.

**Type:** Exposition of Bump-Naprienko (2503.05960), plus the interpretive claim that this answers the "why so many models" question. The claim that *all* known bijections arise as groupoid morphisms is a guiding principle rather than a proven theorem -- flag this honestly.

### 4.2 Classical bijections as groupoid morphisms

**Content:** Illustrate with two examples. (i) Purbhoo's mosaic bijection (puzzles $\leftrightarrow$ tableaux): reinterpret as a sequence of R-matrix moves that transform the KTW lattice model into the Schur lattice model. (ii) Gustafsson-Westerlund's result (2505.07806) that Schutzenberger involution arises as a groupoid morphism. In each case, the "mysterious" bijection becomes a canonical consequence of the YBE.

**Type:** Partly exposition (Gustafsson-Westerlund), partly conjectural (the Purbhoo reinterpretation). *Computation needed:* verify the Purbhoo-as-groupoid-morphism claim for small cases (e.g., $n = 3, 4$). This is a priority test from Clio's research agenda (open question 1 in SUMMARY.md).

### 4.3 The five-vertex model and the three-level integrability hierarchy

**Content:** Present the central structural result of the paper: the negative theorem and the hierarchy it reveals.

**The negative result.** The five-vertex R-matrix arising from Miller's Demazure atom model does NOT satisfy the Yang-Baxter equation. The failure is explicit: 4 of 64 boundary configurations produce nonzero residuals, proportional to $q(p - q)$. No weight assignment can rescue it. This rules out the groupoid composability hypothesis (that Miller's restricted YBE reflects a sub-groupoid of the six-vertex YBE groupoid) and is more informative than a positive result would have been.

**Chirality (revised).** The reason involves both chirality and label asymmetry. The five-vertex model is *chiral* — it distinguishes left-moving paths from right-moving paths. However, chirality alone does not obstruct the YBE: we proved that $R(u) = c(u)(P + \alpha u \cdot E_{22})$ is chiral AND integrable, computing $\omega$-dual Jacobi-Trudi (the $t=0$ specialization of the Wheeler-Zinn-Justin Hall-Littlewood model). The true obstruction is the combination of chirality with label asymmetry ($a_1 \neq a_2$ in the diagonal weights). This is the ONLY nontrivial five-vertex YBE solution, yielding a complete classification of the five-vertex case.

**The three-level hierarchy.** The failure organizes into a clean hierarchy, which we identify with categorical structures:

| Level | Category | Algebra | Symmetry group | YBE status | Example |
|-------|----------|---------|---------------|------------|---------|
| Solvable | Braided monoidal | $H_q$ (generic $q$) | Braid group $B_n$ | Full YBE | Schur polynomials |
| Quasi-solvable | Coboundary monoidal | $H_0$ (0-Hecke monoid) | Cactus group $J_n$ | Column-level | Demazure atoms |
| Combinatorial | Symmetric monoidal | $S_n$ | Symmetric group $S_n$ | Trivial | Crystal bases |

Each level has a different mechanism ensuring consistency:

- **Solvable / braided** (vertex-level): The full YBE holds. Transfer matrices commute. R-matrix moves provide local bijections. The Yang-Baxter groupoid of Bump-Naprienko lives here. The category $\operatorname{Rep}(U_q(\mathfrak{sl}_n))$ is braided monoidal with braid group symmetry. Structure constants (LR coefficients) are nonnegative.

- **Quasi-solvable / coboundary** (column-level): The term "quasi-solvable" is due to Buciumas-Scrimshaw (Forum Math. Sigma). The YBE fails at the vertex level, but column transfer matrices still compose correctly (Miller's Column Lemma 6.3). We proved (Cactus Representation Theorem, §5) that column transfer matrices at $q=0$ satisfy the cactus group relations, making $\operatorname{Rep}(U_0(\mathfrak{sl}_n))$ a coboundary monoidal category. Structure constants can be *signed*.

- **Combinatorial / symmetric** (crystal-level): No R-matrix or column lemma. The crystal functor $F: \operatorname{Rep} \to \operatorname{Crys}$ forgets all non-trivial monoidal structure. Consistency comes from $S_n$ acting on crystal bases. The category $\operatorname{Crys}(\mathfrak{sl}_n)$ is symmetric monoidal.

**The Hecke mirror.** The integrability hierarchy mirrors the algebraic chain $H_q \to H_0 \to S_n$: at generic $q$, the Hecke generator satisfies $(T_i - q)(T_i + 1) = 0$ (two eigenvalues → braided); at $q = 0$, $T_i^2 = -T_i$ (idempotent → coboundary); collapsing further to $S_n$ gives symmetric. Each transition forgets something precise — the content of Section 5.

**The chirality principle (revised).** Chirality + label asymmetry ($a_1 \neq a_2$) obstructs vertex-level YBE. Chirality alone does not: the unique nontrivial five-vertex YBE solution $R(u) = c(u)(P + \alpha u E_{22})$ is chiral and integrable. The principle: *achiral models with $a_1 = a_2$ are solvable with positive structure constants; chiral models with $a_1 \neq a_2$ are at best quasi-solvable with signed structure constants.*

**Type:** The negative result (YBE failure) is a verified computation. The three-level hierarchy synthesizes Miller (2503.09240), Buciumas-Scrimshaw (Forum Math. Sigma), and the Hecke algebra perspective. The chirality principle is a new interpretive claim — flag it as conjectural but well-motivated by the data. *Computation completed:* five-vertex YBE check (4/64 failures, residuals $\propto q(p-q)$).

---

## Section 5. The Categorical Phase Diagram (4-5 pages) — CENTRAL SECTION

*This is the paper's main contribution. Include the phase-diagram.tex figure.*

### 5.1 The Multiplicity Bundle Theorem

**Content:** State and prove the central theorem. In the Schur-Weyl decomposition $V^{\otimes k} \cong \bigoplus_\lambda V_\lambda \otimes S_\lambda$:

1. The coboundary commutor factorizes: $\sigma_{[i,j]} = \bigoplus_\lambda (\operatorname{id}_{V_\lambda} \otimes \rho_{S_\lambda}(\pi(s_{[i,j]})))$
2. The crystal functor sends $\sigma_{[i,j]}$ to the identity (crystal commutor is trivial on identical fundamental factors, by rigidity + multiplicity-freeness of $B(\omega_1)^{\otimes 2}$)
3. The fiber over each crystal component $B(\lambda)$ is the Specht module $S_\lambda$ — the "multiplicity bundle"
4. The three-level hierarchy: braided → coboundary → symmetric = forget spectral parameter, then forget multiplicity bundle

This gives CONCRETE CONTENT to each vertical transition. The proof uses Schur-Weyl duality + Schur's lemma + Henriques-Kamnitzer crystal commutor theory. Computational verification for sl_2 (k=4), sl_3 (k=3,4), sl_4 (k=3).

**Type:** New theorem (proved 2026-04-11). LaTeX proof: `proofs/2026-04-11-multiplicity-bundle.tex` (7pp). Verification: `proofs/verify_multiplicity_bundle.py`.

### 5.2 The two-dimensional phase diagram

**Content:** Present the 2D categorical phase diagram (Figure 1, `phase-diagram.tex`):

- **Vertical axis** (integrability level): braided monoidal → coboundary monoidal → symmetric monoidal
- **Horizontal axis** (deformation type): Classical (Schur), K-theoretic (Grothendieck), Quantum K, Spin (Hall-Littlewood / q-Whittaker)

| | Classical | K-theoretic | Quantum K | Spin |
|--|-----------|-------------|-----------|------|
| **Braided** | Zinn-Justin 2008 | Wheeler-ZJ 2016 | Gorbounov-Korff-Mihalcea 2025 | Gunna-Wheeler-ZJ 2025 |
| **Coboundary** | Cactus Rep Thm | ? | ? | ? |
| **Symmetric** | Crystal bases | K-crystals | ? | ? |

All four columns at the solvable row are established. The coboundary and symmetric rows are populated for classical, partly conjectural for deformations. The deformation axis is orthogonal to the integrability axis.

Key insight: different combinatorial models (puzzles, BPDs, tableaux) are not separate objects connected by bijections — they are specializations of a single integrable lattice model over a parameter space (Fan-Guo-Xiong 2023). The core question shifts from "why so many models?" to "what is the moduli space of specializations?"

**New computational result (04-11):** Fan-Guo-Xiong's R_row and R_col BOTH satisfy pure YBE independently. At $\beta = 0$ (Schubert) they coincide; at $\beta \neq 0$ (K-theory) they differ. Neither handles the full multi-color model (asymmetric cross vertex). This gives clean content to the vertical axis in terms of R-matrix structure:
- **Solvable/braided**: unified integrability (one R-matrix)
- **Quasi-solvable/coboundary**: factored integrability (two independent R-matrices, no unified one)
- **Combinatorial/symmetric**: no R-matrix

The obstruction to unification is the asymmetric cross vertex $L(a,b,a,b) = 1$ only for $a < b$ — a chirality/label asymmetry, the same type of obstruction identified in the five-vertex classification (§4.3).

**Type:** New synthesis + new computation. The factored integrability result is Clio's contribution (script: `r_col_ybe_test.py`).

### 5.3 The Hecke algebra as transition algebra

**Content:** The Hecke algebra $H_q(S_n)$ sits between the braided and coboundary levels — not as a third level but as the *algebra of the transition*. Evidence:

- At generic $q$: Hecke generators have two eigenvalues → braided monoidal structure
- At $q = 0$: Demazure operators $\pi_i = T_i + 1$ become idempotent ($\pi_i^2 = \pi_i$) → coboundary monoidal structure with cactus group symmetry (proved: Cactus Representation Theorem)
- Kalmykov (2025): at the rational level, shuffle algebra = dAHA spherical subalgebra = truncated shifted Yangian. NOTE: the Shibukawa-Ueno hierarchy (elliptic/trigonometric/rational) is an ORTHOGONAL axis to our braided/coboundary/symmetric hierarchy. They parameterize spectral curve vs q-degeneration respectively. Our phase diagram lives at the trigonometric level of Shibukawa-Ueno.
- All five seed paths (puzzle, lattice, Fock space, cylindric, ribbon) converge at the Hecke algebra because they all approach the braided/coboundary boundary.

Crucially (Kamnitzer-Tingley 2007): the coboundary structure exists at ALL $q$ via Drinfeld's unitarized R-matrix. At $q = 0$, the braiding dies but the coboundary survives. The debraiding is extraction, not creation. The phase transition is the death of braided structure, not the birth of coboundary structure.

**Type:** New synthesis combining Kalmykov, Kamnitzer-Tingley, and our cactus computation. The "transition algebra" interpretation is Clio's contribution.

### 5.4 Supporting theorems

**Content:** Briefly state (proofs in appendix or separate papers):

1. **Operator Independence Theorem**: $\pi_{\text{sort}} \neq \lim R(u)$; the hierarchy is of operators, not parameters. Min poly degree 2 ≠ 4. (Proved 04-10.)
2. **Cactus Representation Theorem**: Column transfer matrices at $q=0$ satisfy cactus group relations for $J_n$, $n \leq 5$. Sharp transition: cactus at $q=0$, braid at generic $q$. (Proved 04-10.)
3. **Five-vertex YBE Classification**: The unique nontrivial five-vertex YBE solution is $R(u) = c(u)(P + \alpha u E_{22})$. (Proved 04-09.)
4. **Cactus non-splitting**: $J_3 \twoheadrightarrow S_3$ does not split. ($J_3 \cong D_\infty$, obstruction $3 \nmid 1$.) Exhaustive evidence for $J_4$. (Proved 04-11.)

**Type:** Brief theorem statements with references to full proofs. Computational verifications in appendix.

---

## Section 6. The Shuffle Algebra Substrate (1-2 pages)

### 6.1 Feigin-Odesskii shuffle algebras

**Content:** Briefly introduce the Feigin-Odesskii shuffle algebra (Garbali-Gunna 2024) as the algebraic layer beneath the groupoid. All vertex models (five-vertex, six-vertex, nineteen-vertex) are different representations of the same shuffle algebra. Different representations yield different families of symmetric functions (Schur, Hall-Littlewood, the new rational symmetric functions of Garbali-Guo-Wheeler 2024). This explains why different lattice models compute related but distinct objects.

**Type:** Exposition of Garbali-Gunna (2024). Keep brief -- this is the deepest algebraic layer and the audience may not have the background for a full treatment. Point the interested reader to the original papers.

### 6.2 From shuffle algebras to the Yangian

**Content:** Mention the connection to the Yangian $Y(\mathfrak{gl}_n)$ via the CoHA (cohomological Hall algebra): Yang-Zinn-Justin (2403.17433) construct higher-spin R-matrices from Yangian representations. This provides the "master R-matrix" that specializes to each level of the hierarchy. State this as a direction rather than a completed story.

**Type:** Brief forward reference. Pure exposition of existing work. 1 paragraph.

---

## Section 7. Open Questions and Future Directions (1-2 pages)

**Content:** Collect the open questions, organized by tractability:

### Testable now
1. **Purbhoo as groupoid morphism.** Does the puzzle $\leftrightarrow$ tableau bijection arise from R-matrix moves in Bump-Naprienko's framework? Checkable for $n \leq 5$.
2. **Chirality + label asymmetry principle.** Test the revised conjecture (achiral + $a_1 = a_2$ = solvable + positive; chiral + $a_1 \neq a_2$ = quasi-solvable + signed) across K-theoretic and Hall-Littlewood families.
3. **KL = YBE computation.** ~~Verify for $S_3, S_4$.~~ DONE ($S_3, S_4$). Extend to $S_5$ (in progress).
4. **Specialization moduli.** Does Fan-Guo-Xiong's $R_{\text{col}}$ satisfy YBE independently of $R_{\text{row}}$? If so, quasi-solvable = column integrability is confirmed. (Computation in progress.)
5. **Cactus non-splitting for $J_n$, $n \geq 4$.** $J_3$ proved. $J_4$: 928K triples searched, no section at length $\leq 5$. Full proof: abstract word / cohomological approach needed. 14-year MO question with 0 answers.

### Structural
6. **Hecke transition algebra.** Make precise: $H_q$ is the interpolation algebra between braided and coboundary, with $q$ as the transition parameter. The trigonometric level in the Shibukawa-Ueno hierarchy. (Partially established by our theorems; full statement = next prove session.)
7. **Algebraic content of the Column Lemma.** What algebraic structure (beyond $H_0$) governs the $\pm x_c$ cancellation? Is there a "column-level groupoid"?
8. **Shibukawa-Ueno $\leftrightarrow$ phase diagram.** Do the three R-operator levels (elliptic/trigonometric/rational) map to three monoidal category types? (Depends on Kalmykov reading.)
9. **Master R-matrix.** Is there a single R-matrix (from the Yangian or shuffle algebra) that specializes to each level?

### Frontier
10. **Types B, C, D.** Extend the coboundary/cactus framework beyond type A. Does coboundary = $q=0$ hold in type D? (Svyatnyy 2504.14344, Brown-Elek-Halacheva 2412.02614.)
11. **3D Fock space.** Does Padmanabhan-Korepin's Majorana fermion solution to ZTE lift the fermionic Fock space to 3D? Specialization to $q=0$ might "flatten" to 2D Schur/LR story.
12. **KPZ for puzzles.** Do random large KTW puzzle tilings exhibit Tracy-Widom fluctuations? (Via stochastic six-vertex connection.)
13. **Geometric complexity theory.** LR coefficients appear in Mulmuley-Sohoni's approach to P vs NP. Does the integrability structure have complexity-theoretic consequences?

**Type:** Pure exposition / research program. No computation needed.

---

## Appendix A. Computational Verification (optional, 2-3 pages)

**Content:** Self-contained code listings (Python/SageMath) for: (i) enumerating KTW puzzles via the row transfer matrix; (ii) verifying the YBE for KTW tile weights; (iii) computing LR coefficients via Kostka matrix inversion in the Fock space. Reference Robin's codebase (`transfer_operators.py`, `lr_coefficients.py`) and Clio's verification scripts. This appendix makes the paper computationally reproducible.

**Type:** Computation. Code already exists in the puzzles codebase; needs cleaning and annotation for exposition.

---

## Summary of computation needs

| Section | Computation | Status | Tool |
|---------|-------------|--------|------|
| 2.2 | Puzzle diagrams for $c^{(3,2,1)}_{(2,1),(2,1)} = 2$ | Available via `lr_coefficients.py` | Python |
| 3.1 | R-matrix / YBE verification for KTW tiles | Straightforward symbolic check | SymPy or SageMath |
| 3.3 | Transfer matrix commutativity demo | Can extract from existing code | Python |
| 4.2 | Purbhoo bijection as groupoid morphism | **New computation needed** | SageMath |
| 4.3 | Five-vertex YBE failure (4/64 configs, residuals $\propto q(p-q)$) | **Completed** | SymPy |
| 5.1 | KL = YBE verification for $S_3, S_4$ | **New computation needed** | SageMath |
| Appendix | Clean code listings | Editing/annotation of existing code | Python |

Sections marked **new computation needed** are open research computations; the 4.3 computation (five-vertex YBE failure) is completed and constitutes a key negative result of the paper. The remaining sections are expository synthesis.

---

## Estimated page counts

| Section | Pages |
|---------|-------|
| 1. Introduction | 2 |
| 2. Classical LR coefficients | 3-4 |
| 3. Integrable lattice model perspective | 3-4 |
| 4. Yang-Baxter groupoid | 3-4 |
| 5. **Categorical phase diagram (central)** | **4-5** |
| 6. Shuffle algebra substrate | 1-2 |
| 7. Open questions | 1-2 |
| Appendix A | 2-3 |
| **Total** | **19-24** |

---

## Key references

- Buciumas-Scrimshaw, *Quasi-solvable lattice models*, Forum Math. Sigma
- Bump-Naprienko, *The six-vertex Yang-Baxter groupoid*, arXiv:2503.05960 (2025)
- Brubaker-Bump-Hardt-Spink, arXiv:2509.17312 (2025)
- Gaetz-Gao, arXiv:2408.07863 (2024)
- Garbali-Gunna, *Shuffle algebras and vertex models* (2024)
- Garbali-Guo-Wheeler, arXiv:2412.18085 (2024)
- Gunna-Wheeler-Zinn-Justin, arXiv:2512.04468 (2025)
- Gustafsson-Westerlund, arXiv:2505.07806 (2025)
- Knutson-Tao-Woodward (2004)
- Knutson-Zinn-Justin, arXiv:2509.01857 (2025)
- Miller, arXiv:2503.09240 (2025)
- Purbhoo, arXiv:0705.1184 (2007)
- Wheeler-Zinn-Justin, arXiv:1508.02236 (2016)
- Yang-Zinn-Justin, arXiv:2403.17433 (2024)
- Zinn-Justin, arXiv:0809.2392 (2008)

### Added for revised Sections 5 and 7:
- Alqady-Stroinski, arXiv:2502.05732 (2025) — TL₀ is coboundary monoidal
- Fan-Guo-Xiong, arXiv:2309.00467 (2023) — BPD meets puzzles (specialization moduli)
- Gorbounov-Korff-Mihalcea, arXiv:2503.08602 (2025) — Quantum K = YBE algebra
- Gunna-Wheeler-ZJ, arXiv:2504.19205 (2025) — Spin HL structure constants
- Henriques-Kamnitzer, arXiv:math/0501060 (2006) — Crystal commutor
- Kalmykov, arXiv:2505.09520 (2025) — Shuffle = DAHE = truncated shifted Yangian
- Kamnitzer-Tingley, arXiv:0707.2248 (2007) — Unitarized R-matrix, coboundary at all q
- Liao-Rybnikov, arXiv:2506.16561 (2025) — Cactus maximally transitive on SYT
- Savage, arXiv:0804.4688 (2008) — Survey: crystals are coboundary, not braided
