# Onishi determinants for the complete weighted mode products

This is a new proof notebook connecting the paper's sigma mode solutions to its normalized Onishi determinants. It derives the interpolation and composite-factor construction and checks the weighted version for both worked examples.

Run from the repository root:

```sh
.venv/bin/marimo edit --watch code/hyperelliptic_scalar_first/onishi_products/notebook.py
```

`notebook.html` is the executed walkthrough. To regenerate it:

```sh
.venv/bin/python -m marimo check code/hyperelliptic_scalar_first/onishi_products/notebook.py
.venv/bin/python -m marimo export html code/hyperelliptic_scalar_first/onishi_products/notebook.py -o code/hyperelliptic_scalar_first/onishi_products/notebook.html --force
```

The notebook imports the exact polynomial models from `code/hyperelliptic_scalar_first/model.py`; all determinant algebra, sigma evaluation and product checks are implemented in the notebook. It does not import the old determinant notebooks or their verification reports. It does not rerun RK4 or alter either trajectory notebook.

## How the matrix size is chosen

The known inputs are the genus, the scalar polynomial `q` and the Hamiltonian's mode count. An Onishi determinant must use an initial segment of the monomials ordered by their poles at infinity. The first list containing the ordinate has `g+2` columns and accommodates `degree(q)<=g`. The next polynomial column permits `degree(q)=g+1`, giving size `g+3`. Both lists have only one ordinate column, so its moving-row cofactor is constant. This notebook explains these short lists before point counting or composite sigma factors. The companion [higher_degree](../higher_degree/README.md) package treats longer lists and reduces their mixed ordinate columns using the fixed-point interpolation conditions.

## The proof

1. The scalar identity `q^2-F=4r` identifies the full zero divisor of `q-rho'`, including the multiplicities supplied by the weights.
2. The pole-ordered Onishi matrices use the ordinate `rho'/2`, making their convention exactly the author's monic curve convention. Interpolation recovers the actual `q` and gives the correct column sign and constant Vandermonde minor.
3. An auxiliary polynomial-column Onishi identity splits the composite sigma factor into the residual elementary factors.
4. Comparing the logarithmic derivative of the resulting complete raw product with the scalar-chain-rule identity proves `sum_j beta_j(c_jk-h_jk)=0` on balanced fixed Abel lifts. This is an analytic coefficient identity, not numerical parameter fitting.
5. The paper's initial product condition fixes the proportionality constant. In the present higher-degree interpolation family, `V=D/Delta` and `U=-D*/Delta`, where only the fixed points are moved to the opposite hyperelliptic sheet in `D*`. The star denotes the sheet involution, not complex conjugation.

The new Schur-complement subsection first moves the ordinate column to the last position. The resulting top-left block is Vandermonde. In that order, the determinant divided by the Vandermonde determinant is exactly `(rho'-q)/2`; the original Onishi column order supplies the stated sign. The fixed-ordinate vector is unhalved, so `p(rho)^T W^{-1} rho'_A=q(rho)`.

For genus two, the main points are modes 1–4 and the residual points are modes 1 and 5. For genus three, the main points are modes 1–5 and the residual points are modes 1, 6 and 7. Each list is internally distinct. Sharing mode 1 between the two lists accounts for its double weight without repeated rows in either determinant.

## Executed checks

- **40 exact checks** for the actual two examples: scalar contact multiplicities, Vandermonde-inverse recovery of q, Schur complements and block determinants, column-order signs, normalizing minors and auxiliary Vandermonde identities.
- **4 additional exact checks** of the lower-degree interpolation column orientation, in genera 2–5.
- **180 numerical checks** at 30 dps, including both sheets, exterior coordinates, and the bounded genus-three branch.
- Maximum scaled errors: **9.19e-30** in genus two and **3.33e-28** in genus three.
- The weighted gauge-coefficient residual is at most **1.49e-29**. The sigma-solution product versus normalized determinant error is at most **3.66e-29**.

The main Vandermonde minors are exactly `21/50000` (genus two) and `81/4000000000` (genus three). Signs and composite-factor constants come from the corrected identities, not calibration. The integer Abel-lattice path class is verified before balancing a unit-weight lift; this changes its integration path, not the curve point.

See `proof_scope.md` for the scope of this notebook. Companion checks of
the paper's branch-point and longer-column arguments are in
[branch_points](../branch_points/README.md) and
[higher_degree](../higher_degree/README.md). The notebook and paper use the
fixed-point notation `rho(mu_j), rho'(mu_j), w(mu_j)`.
