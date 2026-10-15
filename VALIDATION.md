# Repository validation

## Audit: 15 October 2026

- Python 3.14.7 on macOS; installed dependencies match every recorded version and both pinned Git commits. `python -m pip check` passes. A fresh network installation was not performed during this audit.
- `python paper/verify_paper.py` passes, including in isolated Python mode from outside the repository.
- `PYTHONPATH=code python -m hyperelliptic_scalar_first.verification --dps 30` passes: 46 exact checks, 441 genus-two checks and 637 genus-three checks.
- `python paper/check_scalar_inversion.py` passes all eight coordinate checks.
- `python -m marimo check` passes for all four notebooks. All Python sources compile.
- Fully executed all four notebooks with `python -m marimo export html` in a separate copy under `/tmp/hyp-ham-audit`, preserving supplied HTML and figures. All exports exited successfully.
- Genus two passes all three RK4 grids through 131,072 steps; finest maximum scaled mode error: `8.519223872e-23`.
- Genus three passes all three RK4 grids through 98,304 steps, including 18 turning/zero observations; finest maximum scaled mode error: `9.959331834e-16`.
- The determinant notebook passes 40 exact checks, four lower-family checks and 180 numerical checks. Its recomputed JSON report matches the supplied report exactly.
- Forced a rebuild of the 33-page paper in the separate copy with `latexmk -g -pdf -interaction=nonstopmode -halt-on-error main.tex`. No LaTeX warnings, undefined references or overfull/underfull boxes were reported. This is a build check, not a visual or mathematical peer review.
- No absolute user-directory paths were found in Python sources or supplied HTML exports. The source/document scan found no obvious credential markers.
- Added ignore rules for macOS metadata and LaTeX build intermediates; refreshed the stale checksum manifest.

There is no separate pytest suite; verification uses the scripts and notebook assertions above. The short verification and scalar-inversion reports were refreshed by their documented commands. Supplied notebook exports and trajectory reports were retained.

Before publication, choose a license (none is supplied) and make the initial Git commit (all repository content is currently untracked). These decisions remain with the author.

## Paper software acknowledgements: 15 October 2026

Added the project repository URL near the end of the paper and Generapy/mpmath references in the introduction, including the literature-reproduction checks performed during Generapy development. The mpmath citation follows the project website's format with the version used here (1.4.1). Rebuilt the updated 34-page paper, checked the affected rendered pages, and reran the paper verification and checksum checks.

The introduction's literature review was subsequently revised to explain the roles of hyperelliptic Kleinian functions, inversion on strata, prime forms and interpolation determinants. Added the Eilbeck--Enolskii--Previato (2003) addition-formula reference. The updated PDF was rebuilt and the affected pages visually checked; paper verification and manifest checks pass.

The ``Our contribution'' subsection was also rewritten to retain the Hamiltonian with explicit definitions of its complex modes, paired powers, integer weights and polynomial, while describing the reduction and reconstruction in words. The rebuilt introduction was visually checked, and paper verification and manifest checks pass.

## Original packaging record

- Rebuilt `paper/main.pdf` from the copied LaTeX and five figure assets within this bundle.
- Ran `paper/verify_paper.py` in isolated Python mode from outside the repository; all printed-algebra, saved-result, figure and reference checks passed.
- Ran `marimo check` on all four copied notebooks; all passed.
- Confirmed copied source/docs have no absolute workspace paths or required references to exploratory files outside the bundle.
- Python checks used the existing environment recorded in `environment.json`. A fresh network installation of `requirements.txt` was not performed.
- Long trajectory computations were not rerun for packaging. Executed HTML and numerical reports are included from the verified example runs; supporting notebook prose may retain historical presentation in those snapshots.
- Original workspace files were copied, not moved or deleted. Paths and selected documentation were adjusted only in this bundle.

`MANIFEST.sha256` records the packaged file hashes, excluding itself. Regenerate it after editing the bundle.
