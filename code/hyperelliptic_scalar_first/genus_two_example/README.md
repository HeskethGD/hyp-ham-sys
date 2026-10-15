# Genus-two worked example: extended propagation

This folder isolates the first worked example for review. Its marimo notebook contains the complete calculation, with code and commentary explaining the exact algebra, analytic identities, numerical evaluation, scalar inverse quadrature, and independent RK4 propagation. The original combined notebook and genus-three example are unchanged.

Open the interactive notebook from the repository root:

```sh
.venv/bin/marimo edit --watch code/hyperelliptic_scalar_first/genus_two_example/notebook.py
```

`notebook.html` is an executed export including source code and outputs. Regenerate it with:

```sh
.venv/bin/python -m marimo check code/hyperelliptic_scalar_first/genus_two_example/notebook.py
.venv/bin/python -m marimo export html code/hyperelliptic_scalar_first/genus_two_example/notebook.py -o code/hyperelliptic_scalar_first/genus_two_example/notebook.html --force
```

The fine RK4 runs take a few minutes. The experiment settings are ordinary constants near the top of the notebook, so expensive propagation is not tied to a frequently changing UI widget.

## What is prescribed, and what is numerical?

- The Hamiltonian, weights, offsets, restricted curve, partial derivatives, and polynomial corrections are defined or derived with exact rational algebra. The fast numerical field is checked symbolically against direct differentiation of the full Hamiltonian.
- The sigma mode solutions and fixed-point formulas for the exponential coefficients are analytic expressions. Generapy supplies numerical period data, characteristics, forward Abel maps, second-kind primitives, and sigma jets. No coefficients are fitted to a trajectory.
- Scalar root finding inverts a specified quadrature for a prescribed propagation value. It does not estimate Hamiltonian or gauge parameters. The quadrature is checked at two Gauss–Legendre orders and against the first Abel increment.
- RK4 integrates the full polynomial field independently. It receives initial modes, duration, and step count; it never calls sigma, the Abel map, or the scalar inverse during its integration loop.

The notebook uses 30-digit arithmetic and propagates on the positive sheet from `rho=4.5` to `rho=10`, staying outside the largest real branch point. All five initial seed amplitudes receive the same sixth-root scale to satisfy the initial product condition. This is a declared choice of initial data, not a fit; it differs from the original plot's mode-amplitude choices while preserving the same Hamiltonian and invariant offsets.

## Sampling and error measurements

The analytic solution is evaluated at 33 equally spaced propagation values. RK4 is run on 32,768-, 65,536-, and 131,072-step grids. These grids pass exactly through all analytic observation times, so the comparison uses no interpolation. The finest run saves 513 states for plotting; it does not store every integration step.

The reported mode error is the maximum over all ten mode components and all 33 analytic observations of `abs(RK4 - sigma) / max(1, abs(sigma))`. Hamiltonian drift is scaled by the largest Hamiltonian term; weighted-invariant drift is scaled by its conserved offset. The notebook checks step-halving ratios and requires the finest mode error to be below `1e-20`.

Working precision and RK4 truncation accuracy are distinct. Smaller integration steps reduce fourth-order truncation error without requiring sigma evaluation on a dense grid. The convergence plot exposes the remaining error instead of hiding it beneath overlapping trajectory curves.

## Artifacts

- `notebook.py`: complete, commented marimo source.
- `notebook.html`: executed walkthrough with embedded code and figures.
- `results.json`: precision, versions, numerical settings, analytic checks, convergence, timings, and curvature diagnostics.
- `u_trajectories.png`, `v_trajectories.png`: real parts of every mode, with thin dashed black RK4 lines and black analytic markers.
- `rk4_refinement.png`: scaled discrepancies on all three integration grids.

The mode plots use actual amplitudes, without independently normalizing the curves. The endpoint-chord statistic is a curvature diagnostic, not a fitted line or part of the solution construction. These computations cover this regular genus-two path; they do not test turning-point passage or global period continuation.

## Combined paper layouts

`modes_combined.png` groups all real mode components on a single page (three columns, four rows in genus two; three columns, five rows in genus three). In genus three, `scalar_powers_combined.png` places the scalar coordinate above the seven modal powers. The refinement figures remain separate.

`plot_data.json` saves the plotted samples after numerical evaluation. Layouts can be regenerated without recomputing trajectories using `render_paper_plots` in `code/hyperelliptic_scalar_first/paper_plots.py`. The notebook also generates these layouts on execution.
