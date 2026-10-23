# Scope and gaps in the weighted determinant connection

The product connection is complete for both current worked examples. The notebook gives an analytic proof of each step and executes exact and numerical checks of the actual polynomial models.

The crucial new observation is that a double-weight point can belong once to the main fixed-point list and once to the residual list. Both lists then remain internally distinct. The earlier composite-factor proof extends to these lists without treating a squared mode factor as a new dynamical mode or inserting two identical determinant rows.

The identity `q^2-F=4r` supplies the contact order of each fixed point, hence the weighted Abel divisor relation. The Vandermonde inverse recovers the actual restriction `q` from its fixed ordinates. After moving the ordinate column to the end, the Schur complement is `(rho'-q)/2`, with the original determinant recovered by its column-order sign. The main determinant therefore represents the actual Hamiltonian restriction `q`; its sheet-involution counterpart uses opposite fixed ordinates with the moving point unchanged. The auxiliary Vandermonde determinant resolves the summed sigma factor. An analytic comparison of logarithmic derivatives proves the weighted gauge cancellation. The initial product condition, rather than an empirical scale, gives the final normalized `U,V` determinants.

Three convention issues are resolved explicitly:

- The determinant's ordinate is `rho'/2`, so the curve is monic in Onishi's sense and the Abel differential remains `d rho/rho'`.
- The ordinate column precedes the highest polynomial column in the present examples. Its column permutation gives `V=D/Delta` and `U=-D*/Delta`.
- Fixed Abel lifts are balanced only after their weighted sum has been verified to lie in the period lattice. Sigma normalizers and second-kind coefficient formulas are then evaluated on those same lifts. The raw weighted gauge sum need not vanish on arbitrary principal lifts.

The proof covers the two standard interpolation families, `degree(q)<=g` and `degree(q)=g+1`, assuming the complete weighted point multiset can be split into a main list and a residual list with distinct coordinates within each. The residual list has `g` points. Compatible balanced weighted lifts and regular fixed sigma normalizers are also required. Both actual examples satisfy these assumptions.

The following are extensions, rather than gaps in either current example:

1. **Higher scalar degree.** This notebook tests only the two short column lists above. The longer pole-ordered lists, their nonzero interpolation minors and sigma factorization are now treated in the companion [higher_degree](../higher_degree/README.md) package.
2. **Forced repeated rows.** Weights or degeneracies that cannot be shared between two internally distinct lists require confluent Onishi limits and Hermite interpolation. Simply duplicating a row gives a zero determinant and is invalid. Those limits are not proved in this notebook.
3. **Obstructed balanced mode lifts.** If all weights have a nontrivial common divisor, it may be impossible to absorb a given period vector using identical lifts for each mode's copies. One must retain quasi-period exponential multipliers or justify a lift choice for individual divisor copies. The examples each have unit-weight points, so this obstruction does not arise.
4. **Reconstructing the entire Hamiltonian from the determinant.** Interpolation gives `q` on the invariant power line. It does not specify the off-line polynomial `Q` or its independent partial derivatives, which control the individual phase terms. Many different `Q` share the same `q` and determinant. The paper supplies `Q` before forming these products.
5. **General global continuation.** The determinant products are meromorphic curve functions. This does not make every individual sigma mode globally single-valued in the propagation variable. The existing three-cycle experiment and its mode multipliers remain consistent with this distinction.

Primary reference: [Onishi, Theorem 7.2 and Lemma 7.1](https://arxiv.org/pdf/math/0105189), with the author's [published errata](https://math-meijo-u.heteml.net/research/pub/correct_pem.pdf). The distinct-residual proof retains the mode weights. The numerical checks instantiate the corrected sign directly and do not fit it from determinant values.
