# Hyperelliptic sigma solutions of polynomial Hamiltonian systems

Paper and computational examples by Graham Hesketh. This directory is self-contained and can be used as the root of a separate repository.

The paper derives scalar hyperelliptic motion from polynomial canonical Hamiltonians, reconstructs individual modes with sigma functions and Abelian integrals, and relates their weighted products to Ōnishi determinants. It also classifies special Hamiltonian polynomials and derives self-generating second-order polynomial potentials. Worked examples cover genera two and three, including an oscillatory trajectory.

## Contents

- `paper/main.tex`, `paper/main.pdf`: manuscript source and compiled paper.
- `paper/verify_paper.py`: fast checks of printed algebra, saved results, figure assets and references.
- `paper/check_scalar_inversion.py`: additional genus-three and genus-four coordinate-ratio checks.
- `code/hyperelliptic_scalar_first/genus_two_example/`: the five-pair genus-two example.
- `code/hyperelliptic_scalar_first/genus_three_example/`: the seven-pair oscillatory genus-three example.
- `code/hyperelliptic_scalar_first/onishi_products/`: determinant construction and product-factorization notebook.
- `code/hyperelliptic_scalar_first/special_q/`: weighted special-polynomial classification, exact checks and independent Hamiltonian/second-order integration.
- `code/hyperelliptic_scalar_first/sigma_modes/`: sigma reconstruction and analytic derivative checks for the three special-polynomial examples.
- `code/hyperelliptic_scalar_first/notebook.py`: supporting analytic derivations, including alternative coefficient formulas and short checks on both sheets.
- `code/hyperelliptic_scalar_first/model.py`, `sigma_jets.py`, `verification.py`: shared algebra, sigma coefficient calculations and verification.
- `code/hyperelliptic_scalar_first/paper_plots.py`: publication layouts, regenerated from saved plot data.
- `requirements.txt`, `environment.json`: direct dependencies and the recorded computation environment.

Each notebook directory includes an executed `notebook.html`, its `results.json`, and relevant figures. No execution is needed to read the HTML or PDF. The source notebooks are authoritative; the HTML files record the executed walkthroughs. The paper uses the combined figure layouts.

## License

The Python code, notebook source files and accompanying software documentation are available under the [MIT license](LICENSE). See [LICENSE_SCOPE.md](LICENSE_SCOPE.md) for the exact scope. The manuscript, figures, generated HTML exports and numerical reports are excluded from this software license; no manuscript license is selected here. Third-party dependencies retain their own licenses.

## Python setup

Python **3.14** was used for the supplied results. Git is needed to install the two pinned source dependencies. From this repository root:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Generapy and SymPy are pinned to the exact Git commits used for these calculations. The remaining direct dependencies are pinned to their recorded versions. `environment.json` records those versions; this is not a platform-wide lockfile.

## Read and run the examples

After setup, from the repository root:

```sh
python -m marimo edit code/hyperelliptic_scalar_first/genus_two_example/notebook.py
python -m marimo edit code/hyperelliptic_scalar_first/genus_three_example/notebook.py
python -m marimo edit code/hyperelliptic_scalar_first/onishi_products/notebook.py
python -m marimo edit code/hyperelliptic_scalar_first/special_q/notebook.py
python -m marimo edit code/hyperelliptic_scalar_first/sigma_modes/notebook.py
```

Launch one notebook at a time. The supporting derivation notebook can be opened in the same way at `code/hyperelliptic_scalar_first/notebook.py`.

To execute and export a notebook, for example:

```sh
python -m marimo export html code/hyperelliptic_scalar_first/genus_two_example/notebook.py -o code/hyperelliptic_scalar_first/genus_two_example/notebook.html --force
```

Full runs recompute periods, sigma samples and RK4 comparisons and may take several minutes. They overwrite generated artifacts in their own directories. The general-mode examples use 30-digit working precision; the special-polynomial quadrature checks use 50 digits and their sigma reconstruction uses 35 digits. Observed discrepancies are reported separately from working precision.

The scalar coordinate is obtained by numerical inversion of its Abelian integral. Generapy evaluates the Abel vectors and sigma mode formulas; an independent full-system RK4 calculation supplies the comparison. No mode or exponential coefficients are fitted to trajectory data. The scalar sigma ratio is a coordinate-representation check.

## Verification

Fast paper checks, using the supplied numerical reports:

```sh
python paper/verify_paper.py
```

Recompute the additional scalar coordinate checks:

```sh
python paper/check_scalar_inversion.py
```

Recompute the supporting short-path verification:

```sh
PYTHONPATH=code python -m hyperelliptic_scalar_first.verification --dps 30
```

The worked-example and determinant notebooks execute their own checks. Saved numerical reports describe the scope of each calculation; the short-path verification complements the longer trajectories.

Recompute the special-polynomial and second-order checks:

```sh
PYTHONPATH=code python -m hyperelliptic_scalar_first.special_q.verification
PYTHONPATH=code python -m hyperelliptic_scalar_first.sigma_modes.verification
```

Each package has a README describing its checks and command options. Both
commands write `results.json` in their own package directory.

## Rebuild the paper

Install a LaTeX distribution with `latexmk` and `pdflatex` (for example TeX Live or MacTeX). Python requirements do not install LaTeX. From the repository root:

```sh
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

All figure paths point into `../code/hyperelliptic_scalar_first/`; no original workspace files are needed.

## Change plot layouts without recomputing solutions

The two example folders contain `plot_data.json` with the plotted real mode components, powers, scalar coordinate and sample times. To regenerate the combined paper layouts from these saved data:

```sh
PYTHONPATH=code MPLBACKEND=Agg python - <<'PY'
from hyperelliptic_scalar_first.paper_plots import render_paper_plots
for genus in ('two', 'three'):
    render_paper_plots(f'code/hyperelliptic_scalar_first/genus_{genus}_example')
PY
```

The refinement figures are separate. Layout changes do not alter the Hamiltonians, trajectory samples or numerical comparisons.

## Scope

The paper treats smooth odd-degree hyperelliptic curves, with the fixed-point and polynomial-correction conditions stated in the manuscript. Local formulas and continuation on the genus-three example do not assert global single-valued meromorphicity in propagation for arbitrary curves.
