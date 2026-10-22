# Repository validation

## Readable coefficient table: 22 October 2026

Transposed the genus-three correction table so modes are columns and
polynomial coefficients are rows. Enlarged the entries and used display
fractions with increased row spacing. Updated the paper verifier to read
the transposed table and display fractions; all 21 coefficients still
reconstruct the exact symbolic corrections. Rebuilt and visually checked
the 39-page manuscript; paper verification and checksums pass.

## Mode-factor exposition: 22 October 2026

Moved the Baker--Akhiezer interpretation beside the corrected point-factor
definition in Section 4.4 and the polynomial-correction integral to Section
4.5. Section 4.6 now links back briefly to the normalized Baker--Akhiezer
factors. Removed the redundant substitution sentence and made the barred
uncorrected-factor terminology consistent in the determinant argument and
worked examples. Corrected the paragraph boundary around the genus-three
coefficient table so it no longer splits a sentence.

The rebuilt manuscript remains 39 pages. Paper verification, reference
resolution, visual checks of the affected sections and checksum checks pass.
No mathematical formulas or computational algorithms were changed.

## Simple branch fixed points: 22 October 2026

Extended the main manuscript's local sigma mode construction to simple
branch fixed points. Added the arbitrary-genus continuation argument using
Onishi's nonvanishing theorem, clarified weights as multiplicities after a
compatible degenerate-mode reduction, and revised the determinant divisor
count and opposite-lift convention. Smoothness forces unit weight at a
branch fixed point; nonzero Wronskians remain necessary for independence.

The exploratory work had passed 1,003 mode/turning-point checks and 292
additional determinant checks, including exact interpolation in genera
two through six and full numerical sigma identities in genera two and
three. Those calculations remain local exploratory material; the manuscript
extension is justified analytically and does not depend on their files.

Rebuilt the 39-page manuscript, visually checked the affected assumptions,
kernel proof, mode theorem, determinant argument and discussion, and reran
`paper/verify_paper.py`. All checks pass; the LaTeX log has no undefined
references or overfull/underfull box warnings. Refreshed the manuscript and
validation entries in the checksum manifest.

## Consolidated computational material: 22 October 2026

Moved the special-polynomial and sigma-mode packages into
`code/hyperelliptic_scalar_first/special_q/` and
`code/hyperelliptic_scalar_first/sigma_modes/`, alongside the general
worked examples and determinant notebook. Updated imports, notebook
bootstrap paths, CLI commands, README navigation and paper-verification
paths. The packages retain their models, proof/normalisation notes,
executed HTML walkthroughs and saved numerical reports.

Both marimo notebooks pass static checks and were executed to regenerate
their HTML and reports from the new locations. The runs pass 154 exact
algebra checks, the independent Hamiltonian/second-order comparisons,
30 exact normalisation checks and 680 numerical sigma checks. The paper
verifier passes with the consolidated paths, including isolated Python
execution; both CLI entry points and Python compilation checks pass.
No numerical algorithms were changed.

Added `/part_2/` to `.gitignore` and removed its previously staged files
from Git's index. Historical notes and the separate draft remain locally;
the paper and public code have no dependency on that folder.

## Clarified worked examples: 22 October 2026

Rewrote the opening sentences of Sections 7.1--7.3 to identify the
theorems being demonstrated. Section 7.3 now displays the second-order
mode equations, explicit cubic and quintic potentials, and a polynomial
equation in the paired power for the first example. Numerical curve
normalisation and verification summaries are in Section 8.3; example
prose no longer uses unexplained `LC` notation or implementation paths.
Extended `paper/verify_paper.py` to check all newly printed potentials
and the paired-power equation against exact symbolic expressions.
Verification passes, the 41-page PDF builds without LaTeX warnings, and
the revised examples and validation subsection were visually checked.

## Integration of special-polynomial results: 22 October 2026

Integrated the weighted special-polynomial classification into Section 4 of
`paper/main.tex`, using the common-translation uniqueness proof and a direct
generating-function existence proof. Added the Chan--Chyan--Srivastava
identification and its references. A new Section 6 combines constant
Wronskians, prescribed-product polynomiality and explicit self-generating
potentials. Section 7.3 adds the three verified special-polynomial examples
from Part 2 without repeating the sigma construction.

Extended `paper/verify_paper.py` to compare the new printed curves and
potential with their exact models, check squarefreeness and validate the
reported counts and sigma error bounds against the saved Part 2 results.
The verification passes. Rebuilt the paper from 35 to 40 pages with no
LaTeX warnings or overfull/underfull boxes, and visually checked the added
material, introduction and references. Existing numerical algorithms and
saved computations were not changed or rerun.

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

At the time of this initial audit, no license or initial Git commit was present. Subsequent licensing changes are recorded below.

## Paper software acknowledgements: 15 October 2026

Added the project repository URL near the end of the paper and Generapy/mpmath references in the introduction, including the literature-reproduction checks performed during Generapy development. The mpmath citation follows the project website's format with the version used here (1.4.1). Rebuilt the updated 34-page paper, checked the affected rendered pages, and reran the paper verification and checksum checks.

The introduction's literature review was subsequently revised to explain the roles of hyperelliptic Kleinian functions, inversion on strata, prime forms and interpolation determinants. Added the Eilbeck--Enolskii--Previato (2003) addition-formula reference. The updated PDF was rebuilt and the affected pages visually checked; paper verification and manifest checks pass.

The ``Our contribution'' subsection was also rewritten to retain the Hamiltonian with explicit definitions of its complex modes, paired powers, integer weights and polynomial, while describing the reduction and reconstruction in words. The rebuilt introduction was visually checked, and paper verification and manifest checks pass.

## Response to pre-submission review: 15 October 2026

Added explicit odd-degree scope to the abstract, a local ODE uniqueness sentence, a distinction from algebraic complete integrability, and a table of the existing scalar-coordinate, paired-power and complete-product checks. Added the original Frobenius--Stickelberger citation and its 2026 English translation, crediting both translators. The paper verification script now compares every entry of the new table with the saved numerical reports. No numerical algorithm or trajectory data was changed.

Applied MIT to the software with exact coverage and exclusions in `LICENSE_SCOPE.md`; the manuscript and figures remain outside that license. Recorded a point-by-point response in `reviews/response_to_critical_review.md`. The existing `.gitignore` excludes `reviews/`, so the review and response remain local and are omitted from the public-bundle checksum manifest.

Rebuilt the 35-page PDF, checked the affected rendered pages, and reran paper verification and the checksum checks. The earlier full notebook execution results remain applicable to the unchanged numerical code.

The abstract and introduction were subsequently revised to emphasize the association of determinant identities with conservation laws as a central contribution extending the author's genus-one nonlinear-optics work. Repeated explanations were consolidated into the contribution subsection, and the review response was aligned with that framing. Rebuilt and visually checked the revised introduction; paper verification and checksum checks pass.

## Original packaging record

- Rebuilt `paper/main.pdf` from the copied LaTeX and five figure assets within this bundle.
- Ran `paper/verify_paper.py` in isolated Python mode from outside the repository; all printed-algebra, saved-result, figure and reference checks passed.
- Ran `marimo check` on all four copied notebooks; all passed.
- Confirmed copied source/docs have no absolute workspace paths or required references to exploratory files outside the bundle.
- Python checks used the existing environment recorded in `environment.json`. A fresh network installation of `requirements.txt` was not performed.
- Long trajectory computations were not rerun for packaging. Executed HTML and numerical reports are included from the verified example runs; supporting notebook prose may retain historical presentation in those snapshots.
- Original workspace files were copied, not moved or deleted. Paths and selected documentation were adjusted only in this bundle.

`MANIFEST.sha256` records the packaged file hashes, excluding itself. Regenerate it after editing the bundle.
