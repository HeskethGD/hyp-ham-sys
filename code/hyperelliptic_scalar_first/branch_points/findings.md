# Branch fixed points: computational findings

Historical record of the initial checks. The branch-point results and
arbitrary-genus justification have since been incorporated into the main
manuscript; the code now lives in the main computational package.

The executed notebook evaluates the existing shifted sigma ratios at branch
fixed points directly. It does not substitute square roots for those ratios.
The square-root expression is only an independent local reference.

## What passed

At 35 working digits, 1,003 symbolic/numerical checks passed with numerical
threshold 1e-23. Maximum scaled errors were approximately 1.21e-34 in genus
two and 3.01e-32 in genus three. The curves are

    g=2: F = 4 rho^5 - 20 rho^3 + 16 rho
    g=3: F = 4 rho^7 - 56 rho^5 + 196 rho^3 - 144 rho

Both have Q=0, unit weights, distinct fixed coordinates a=-g,...,g, and
squarefree F. Every fixed mode point is a branch point.

1. All fixed sigma-sharp normalizers were nonzero. Their smallest moduli
   were approximately 0.618 in genus two and 0.233 in genus three.
2. The fixed-point sigma stratum identities and short formulas for h agree
   with independent, compatible-path first/second-kind integration.
3. Independently evaluated u and v sigma ratios agree with local square-root
   ratios and satisfy the canonical equations, paired products, zero
   Wronskians and polynomial second-order equations at regular points on
   both sheets.
4. The same sigma factors, with the extra correction for Q=I1-I2+1, satisfy
   the different canonical field on the same scalar curve. Its restricted
   q is zero but B=(-1,1,0,...,0). Thus the checks are not limited to B=0.
5. Signed local-parameter continuation through the largest branch point,
   including its exact value, gives finite sigma-mode values and derivatives
   satisfying the canonical and second-order equations. The attached mode
   passes through zero with nonzero derivative. With t0=0.15 and u(t0)=1,
   its values at t=(0.15,0.075,0,-0.075,-0.15) are (1,1/2,0,-1/2,-1),
   up to numerical error. The code obtains these from sigma, not square roots.

## The exponential is still needed

The branch condition b=0 does not imply h=0. For the largest branch point,
the h vectors on the tested lifts are approximately

    g=2: (-0.341788328067284, 0.727020016739567)
    g=3: (0.361295378615655, -1.436528402584083, 0.668629666350438).

The raw shifted sigma ratios and natural exponentials must be combined.
Setting h=0 would discard the compensation for sigma quasi-periodicity.
The geometric equality of a branch point and its involution does not justify
identifying their unreduced Abel and second-kind lifts.

## Polynomial potentials and independence

For the mode attached to the largest branch point, the directly tested
potential is

    g=2, a=2: G = 4 rho^3 + 6 rho^2 - 2 rho - 2
    g=3, a=3: G = 6 rho^5 + 15 rho^4 - 20 rho^3 - 45 rho^2 + 8 rho + 12.

These equal (1/4) d/drho [F(rho)/(rho-a)]. Both canonical modes satisfy
the equation, but their Wronskian is zero when B=0. They are proportional;
a second independent solution has not been constructed in this exploration.

## Follow-up: arbitrary genus and determinants

The arbitrary-genus kernel extension and the two existing determinant
families are now treated in [proof.md](proof.md). The proof uses Onishi's
nonvanishing result for sigma-sharp at every finite curve point, local
holomorphic continuation of both kinds of Abel integral, and the existing
nonbranch kernel identity. The branch-point restriction is not needed for
the resulting shifted sigma mode representation.

The determinant divisor count needs a precise replacement: at a simple
branch point, both q-y and q+y have a simple zero in the local curve
parameter. The point is counted once in each divisor. The main column
operation and auxiliary residual factorization then extend unchanged,
provided branch opposite lifts remain minus their original vectors and
the same distinct-list assumptions are retained.

The additional executable [determinant notebook](determinant_notebook.py)
tests the full sigma identity, both fixed-sheet choices, the residual
auxiliary determinant, composite-factor splitting, complete products,
compatible period shifts, and corrected-mode invariance. It includes
numerical mixed branch/nonbranch rows and exact interpolation checks of
both column families in genera two through six. See
[determinant_results.json](determinant_results.json) for the executed report.
All 292 additional checks passed. The numerical maximum scaled errors are
approximately 1.56e-34 in genus two and 2.97e-34 in genus three, at the
same 1e-23 acceptance threshold. Of these checks, 40 are exact symbolic
checks and 252 are numerical sigma/period/determinant checks.

Before balancing, raw products in the tested intervals have defects of
approximately 0.0141 (genus two) and 0.0233 (genus three), measured against
the normalized algebraic product. This is the expected lift-dependent
exponential, not a failure of the corrected modes. With balanced lifts,
the complete natural exponential cancels and the raw products recover
the determinant up to a constant.

## Changes that would be needed in the main paper

- Integrate the arbitrary-genus branch-point specialization and its citations
  from proof.md. Numerical sigma evidence still covers genera two and three.
- Allow b=0 in the mode theorem while retaining b!=0 specifically for the
  assertion that the two canonical modes form a fundamental pair.
- Retain smoothness: F'(a)=-4r'(a) at a branch fixed point. With distinct
  coordinates, a weight greater than one therefore forces a singular curve.
  The exact code confirms this in a double-weight example.
- Replace the determinant proof's nonbranch-only zero-count sentence by the
  local branch calculation in proof.md; retain explicit opposite Abel lifts.
- Global continuation and all-cycle multipliers remain outside these tests.

No main manuscript, production implementation or Generapy source was changed.
