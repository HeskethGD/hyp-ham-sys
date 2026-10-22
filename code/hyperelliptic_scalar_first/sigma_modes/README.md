# Sigma modes on normalised odd-degree curves

This notebook evaluates the corrected sigma factors for the
three models established in `../special_q/`. Both u and v are reconstructed
from independent sigma ratios. First and second derivatives come from
analytic sigma jets, rather than the target differential equations.

From the repository root:

```sh
.venv/bin/python -m marimo edit code/hyperelliptic_scalar_first/sigma_modes/notebook.py
PYTHONPATH=code .venv/bin/python -m hyperelliptic_scalar_first.sigma_modes.verification
.venv/bin/python -m marimo check code/hyperelliptic_scalar_first/sigma_modes/notebook.py
.venv/bin/python -m marimo export html code/hyperelliptic_scalar_first/sigma_modes/notebook.py -o code/hyperelliptic_scalar_first/sigma_modes/notebook.html --force
```

The CLI and notebook recompute checks and write `results.json`; `notebook.html`
is the executed walkthrough. To isolate a model without overwriting the main
report, use e.g. `--example constant_weighted_g2 --output /tmp/sigma_g2.json`.
Default precision is 35 digits. Direct commands accept `--dps` for investigation.

## Normalisation

For original leading coefficient L>0, set s=sqrt(L/4), t=s*z and Y=y/s.
Only the curve and Abel differentials are normalised. Original modes,
powers and Hamiltonian parameters remain unchanged. Then dw_k/dz=s*rho^(k-1).
First sigma derivatives acquire s, second derivatives s^2. All canonical,
Wronskian, and potential comparisons use the original propagation variable z.

## Checks

- Exact leading-four normalisation, curve norm, and potential scaling.
- Fixed-point sigma expressions for h against regularised second-kind Abel
  integrals transported along the same prescribed contour as the first-kind
  integrals; fixed points are also checked on the sigma stratum.
- Moving Abel lifts continued by quadrature versus direct unreduced Abel maps.
- Corrected u and v sigma formulas against independent elementary quadratures.
- Raw and reconstructed paired powers, full weighted products and Hamiltonian.
- First derivatives versus the full canonical field in independent variables.
- Wronskians and second derivatives computed from sigma jets.

The selected curves include constant Q in genera two and three and a
nonconstant special Q with weights (5,5) in genus two. Two noninitial points
on both sheets are tested per curve. The tolerance is 10^(-(dps-12)).
Errors are scaled by max(1, abs(actual), abs(expected)); precision is not
asserted as an achieved accuracy. The general-mode examples are checked separately.

See `normalization.md` for the scaling argument. Global continuation,
branch-point limits, and high-degree determinant identities are not established
by these local checks.

The fixed-point paths start at a common regular exterior anchor and have
three segments above all finite branch points. Root-factor square roots
transport the sheet continuously; quadrature is split near close roots.
This resolves a failed comparison at one genus-three point when separate
automatic endpoint calls were used. Tolerances and sigma formulas were not
changed to accommodate it. See `normalization.md` for details.

## Executed results

All 680 numerical checks and 30 exact normalization checks pass:

| Example | Numerical checks | Maximum scaled discrepancy |
|---|---:|---:|
| Constant Q, weighted genus two | 220 | 1.25e-33 |
| Constant Q, weighted genus three | 334 | 1.75e-31 |
| Nonconstant special Q, weights (5,5), genus two | 126 | 9.27e-33 |

The curve scaling factors are 2, 2 and 25*sqrt(5), respectively. The source
passes `marimo check`, and the executed HTML contains all successful results.
These are observed local discrepancies at 35-digit working precision, not
certified accuracy bounds. See `results.json` for every individual check.
