# Genus-three worked example: three oscillations on a real oval

This folder isolates the paper's second worked example. Its Hamiltonian, weights and invariant offsets are unchanged. We choose initial data on a bounded real branch instead of the original exterior branch. The notebook contains the exact algebra, generapy calls, scalar inverse, continued Abel lift, explicit sigma solutions and independent full-system RK4 integration.

From the repository root:

```sh
.venv/bin/marimo edit --watch code/hyperelliptic_scalar_first/genus_three_example/notebook.py
```

Read `notebook.html` for the executed walkthrough without rerunning the calculation. To validate the notebook and regenerate that export:

```sh
.venv/bin/python -m marimo check code/hyperelliptic_scalar_first/genus_three_example/notebook.py
.venv/bin/python -m marimo export html code/hyperelliptic_scalar_first/genus_three_example/notebook.py -o code/hyperelliptic_scalar_first/genus_three_example/notebook.html --force
```

## Why this branch oscillates

All seven roots of the existing scalar polynomial are real. Between the fifth and sixth roots, approximately `0.3505953510803235` and `0.5739427200213239`, the polynomial is nonnegative. The scalar coordinate travels back and forth between them. Its period is approximately `6.624556082295094`; we propagate for three periods, approximately `19.87366824688528`.

Write `F(rho) = (rho-L)(R-rho) G(rho)` on this interval. The substitution `rho = c-d cos(theta)` gives `dz/dtheta = 1/sqrt(G(rho))`, with `c=(L+R)/2` and `d=(R-L)/2`. Increasing theta passes smoothly through turning points and changes the sign of the curve ordinate. This is a parametrization of the exact scalar curve, not a sinusoidal approximation in propagation time.

The modal powers `I_j=u_j v_j=gamma_j-beta_j rho` are periodic and signed. The sixth and seventh powers cross zero. We evaluate both sigma mode formulas directly and integrate the full polynomial field without dividing by a mode. We do not use logarithmic mode equations across their zeros.

A single generapy forward Abel map initializes the lift. Subsequent Abel coordinates are obtained by quadrature along the actual angular path, preserving the accumulated periods through sheet changes and completed cycles. Calling a fresh principal Abel map at each point would not preserve this continuation automatically.

## What is analytical, numerical, or fitted?

- SymPy derives the restricted curve, phase coefficients, polynomial corrections and canonical equations with exact rational arithmetic. Exact real-root isolation establishes the seven simple real branch points.
- The sigma formulas and fixed-point second-kind coefficient expressions are prescribed analytic formulas. Generapy evaluates periods, characteristics, forward Abel integrals, regularized second-kind integrals and sigma jets numerically.
- The notebook supplies angular quadrature, the scalar inverse and the continued Abel lift. Quadrature is checked at 64 and 96 Gauss–Legendre nodes per quarter-turn panel. Safeguarded Newton iteration inverts a known integral; it does not fit coefficients.
- All seven declared real seed amplitudes receive the same eighth-root scale to enforce the weighted initial product. This is an algebraic choice of initial conditions. No individual mode is singled out.
- RK4 independently integrates all fourteen canonical mode components. Its loop receives only the initial state, polynomial field, duration and step count; it never calls sigma or uses the scalar trajectory.
- Nothing is fitted to an RK4 trajectory. No envelope is removed and no individual plot is rescaled.

These initial data do not impose `v_j=conjugate(u_j)`. Signed paired powers should not be interpreted as positive optical intensities. Periodic paired powers need not give periodic individual amplitudes: reciprocal cycle multipliers can produce an envelope or sign reversal. The notebook reports the actual mode multipliers.

## Sampling, accuracy and artifacts

The default run uses 30-digit arithmetic, 97 equally spaced analytic observations and additional observations exactly at turning points and modal-power zeros. The three RK4 grids have 24,576, 49,152 and 98,304 steps. All 97 comparison times belong to every grid, so there is no interpolation. The finest grid saves 1,537 states for the thin dashed black curves; analytic evaluations appear as black markers.

The scaled mode error is the maximum over all fourteen components and all 97 comparison nodes of `abs(RK4-sigma)/max(1,abs(sigma))`. Hamiltonian drift uses `max(1,abs(U),abs(V),abs(Q))`; weighted-invariant drift uses `max(1,abs(gamma_j-beta_j gamma_7))`. The notebook requires step-halving factors between 12 and 20 and finest mode discrepancy below `1e-12`. The measured errors in `results.json` distinguish arithmetic precision from RK4 truncation accuracy over this much longer path.

The executed run gives:

| RK4 steps | Maximum scaled mode error |
|---:|---:|
| 24,576 | 2.5450e-13 |
| 49,152 | 1.5925e-14 |
| 98,304 | 9.9593e-16 |

The step-halving factors are **15.9808** and **15.9904**. The finest scaled Hamiltonian and invariant drifts are **1.5018e-17** and **3.6849e-18**. Analytic identity errors, including turning and zero observations, remain below **6.0e-29**. The separate final-to-initial paired-power repeat error is **9.16e-29**.

All six turning points and twelve modal-power zero crossings are explicitly sampled by the sigma checks. The 97 uniform sigma samples and extra event evaluations take about 63 seconds; the three RK4 grids take about 13, 26 and 53 seconds. A full run takes roughly three minutes on the recorded environment. The zero-crossing checks supplement the uniform-grid comparisons; RK4 is compared with sigma at the uniform observation nodes.

The saved artifacts are:

- `notebook.py` and `notebook.html`: commented source and executed walkthrough.
- `results.json`: settings, initial data, versions, geometric checks, convergence and cycle multipliers.
- `scalar_oscillations.png`: scalar coordinate over three cycles.
- `powers_trajectories.png`: all seven signed modal powers.
- `u_trajectories.png` and `v_trajectories.png`: all actual real mode amplitudes.
- `rk4_refinement.png`: discrepancies on the three grids.

This tests continuation on one real oval through its turning points and mode zeros. It does not establish an arbitrary global continuation algorithm for all curves. 

## Combined paper layouts

`modes_combined.png` groups all real mode components on a single page (three columns, four rows in genus two; three columns, five rows in genus three). In genus three, `scalar_powers_combined.png` places the scalar coordinate above the seven modal powers. The refinement figures remain separate.

`plot_data.json` saves the plotted samples after numerical evaluation. Layouts can be regenerated without recomputing trajectories using `render_paper_plots` in `code/hyperelliptic_scalar_first/paper_plots.py`. The notebook also generates these layouts on execution.
