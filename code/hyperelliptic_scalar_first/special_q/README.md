# Special Q and polynomial mode potentials

This notebook supports the paper's special-polynomial classification and
second-order mode equations. It includes shared exact models, executable
verification, a marimo walkthrough, an executed HTML export and saved results.

From the repository root, using the existing environment:

```sh
.venv/bin/python -m marimo edit code/hyperelliptic_scalar_first/special_q/notebook.py
PYTHONPATH=code .venv/bin/python -m hyperelliptic_scalar_first.special_q.verification
.venv/bin/python -m marimo check code/hyperelliptic_scalar_first/special_q/notebook.py
.venv/bin/python -m marimo export html code/hyperelliptic_scalar_first/special_q/notebook.py -o code/hyperelliptic_scalar_first/special_q/notebook.html --force
```

`notebook.html` is the executed walkthrough; `results.json` records checks,
example polynomials, coefficient-search dimensions and numerical errors.
Both the CLI and notebook recompute the checks and write the report. For
exact algebra alone use `--exact-only --output /tmp/special_q_exact.json`.

- `model.py`: weighted special-polynomial construction, independent coefficient
  equations, scalar restrictions and full canonical fields.
- `verification.py`: exact checks and short numerical comparisons.
- `proof_scope.md`: arbitrary-degree weighted proof, Wronskian equivalence,
  general second derivatives, polynomial-potential uniqueness and limitations.

The first examples are constant-Q weighted genus-two and genus-three curves,
and a nonconstant special Q with weights (5,5) giving a smooth quintic. The
last example demonstrates that odd scalar degree does not require odd total
weight or constant Q. Its high-degree q needs an extension of the paper's
short determinant construction.

Numerics use 50-digit Abel quadratures versus full canonical RK4, with step
refinement, and a separate second-order RK4 solve. The following notebook in
[`../sigma_modes/`](../sigma_modes/README.md) now carries out leading-four
normalization and independent sigma evaluation on these same curves. Finite
coefficient searches support rather than replace the proof. See
`proof_scope.md` for exact scope.

## Executed results

- 154 exact checks passed, including 60 independent homogeneous coefficient
  searches, the original examples' genera, and constant-Q determinant signs.
- Canonical RK4 refinement factors are 16.00–16.02. At 64 steps, scaled mode
  discrepancies are 4.89e-20, 6.34e-20 and 1.01e-16 for the three examples.
- Separate second-order RK4 mode discrepancies are at most 8.28e-17;
  scaled Wronskian discrepancies are at most 3.40e-22.
- The marimo source passes `marimo check` and has been executed to HTML.

These are local finite-step discrepancies at 50-digit working precision,
not certified bounds or claims of 50-digit solution accuracy.
