import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import sys as _sys
    from pathlib import Path
    example_dir = Path(__file__).resolve().parent
    if str(example_dir.parent.parent) not in _sys.path:
        _sys.path.insert(0, str(example_dir.parent.parent))
    import json
    import marimo as mo
    from hyperelliptic_scalar_first.sigma_modes.verification import run_verification
    return example_dir, json, mo, run_verification


@app.cell
def _(mo):
    mo.md(r"""
    # Sigma reconstruction of the special-Q modes

    The special-polynomial notebook establishes the weighted classification and
    compares elementary Abel-integral solutions with the canonical and
    second-order equations. We now independently evaluate the paper's
    corrected sigma point factors for the same three models: constant $Q$
    in genera two and three, and nonconstant special $Q$ with weights $(5,5)$
    on a smooth genus-two curve.

    Both $u_j$ and $v_j$ are evaluated from their own sigma formulas. We do
    not recover either mode by division by the other. Analytic sigma jets
    supply first and second derivatives, so the Wronskian and potential
    checks test the sigma formulas themselves. No coefficients are fitted.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 1. Normalise the curve without changing the original Hamiltonian

    Let $(\rho'_z)^2=F(\rho)$ with positive leading coefficient $L$.
    Define
    \[
    s=\sqrt{L/4},\qquad t=sz,\qquad
    Y=\rho'_t=y/s,\qquad \widehat F=F/s^2.
    \]
    The curve $Y^2=\widehat F(\rho)$ has leading coefficient four.
    Its fixed points are $(a_j,b_j/s)$, where $b_j=q(a_j)$.
    The Abel coordinates satisfy
    \[
    dw_k=\frac{\rho^{k-1}\,d\rho}{Y},\qquad
    \frac{dw_k}{dt}=\rho^{k-1},\qquad
    \frac{dw_k}{dz}=s\rho^{k-1}.
    \]
    Thus $\Delta w_1=s\Delta z$. First sigma derivatives along the curve
    must be multiplied by $s$, and second derivatives by $s^2$, before
    comparison with the original canonical field and potential.

    The amplitudes and powers are unchanged. We do **not** reinterpret
    $q/s$ as a new canonical Hamiltonian with the same unit product
    coefficients: its norm is
    $(q/s)^2-\widehat F=4r/s^2$. All full-field and product checks below
    use the original $Q$, $q$, $r$, and propagation variable $z$.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 2. Fixed-point lifts and the natural exponential

    Generapy computes first- and second-kind periods on the normalised
    curve. We use the paper's sigma convention, with its algebraic
    second-kind basis and hyperelliptic normalisation. Starting at one regular
    exterior anchor, **both kinds of differential follow the same explicit
    contour** to each fixed point. Its three segments go above all finite
    branch points; continuous root-factor square roots transport the sheet.
    Integration is split near close branch points. Fixed points have
    unreduced Abel lifts $\eta_j$, checked on the sigma stratum. Put
    \[
    \bar u_j=\sqrt{\beta_j}
      \frac{\sigma_\flat(w+\eta_j)}
      {\sigma_\sharp(w)\sigma_\sharp(\eta_j)},\quad
    \bar v_j=(-1)^{g+1}\sqrt{\beta_j}
      \frac{\sigma_\flat(w-\eta_j)}
      {\sigma_\sharp(w)\sigma_\sharp(\eta_j)}.
    \]
    Here $(\sharp,\flat)=(\partial_2,1)$ for genus two and
    $(\partial_2,\partial_3)$ for genus three. The paper's short
    fixed-point formulas determine $h_{jk}$ from sigma derivatives at
    $\eta_j$. We check these independently against
    $-\operatorname{Reg}\int_\infty^{(a_j,b_j/s)}dr_k$ on that prescribed path.

    An initial test using separate automatic endpoint calls disagreed at one
    genus-three fixed point. A common contour, with resolved quadrature near
    close roots, is necessary for this comparison. No coefficient is fitted
    and the acceptance tolerance is unchanged.

    Since $B_j=0$, the full initial-value solutions are
    \[
    u_j=u_{j0}\frac{\bar u_j(w)}{\bar u_j(w_0)}
                  e^{-h_j\cdot(w-w_0)},\qquad
    v_j=v_{j0}\frac{\bar v_j(w)}{\bar v_j(w_0)}
                  e^{+h_j\cdot(w-w_0)}.
    \]
    Initial modes obey the paired powers and complete-product
    constraint of the original Hamiltonian. The moving lift is continued
    by integrating $d\rho/Y$ along each short path; its endpoints are
    checked against independent unreduced Abel-map calls.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 3. Derivatives from sigma jets

    Write $A=\sum_{k=0}^{g-1}\rho^k\partial_{w_{k+1}}$.
    For any nonvanishing sigma derivative $f(w)$,
    \[
    \frac{d\log f}{dt}=A\log f,\qquad
    \frac{d^2\log f}{dt^2}
    =\sum_{i,j=0}^{g-1}\rho^{i+j}\partial_{w_{i+1}}\partial_{w_{j+1}}\log f
     +Y\sum_{k=1}^{g-1}k\rho^{k-1}\partial_{w_{k+1}}\log f.
    \]
    We apply these to each numerator and denominator and differentiate the
    $h_j$ correction explicitly. Order-three sigma jets suffice. The test
    of $u_j''/u_j$ uses $(\log u_j)''+((\log u_j)')^2$ in physical units,
    independently of the prescribed polynomial
    \[
    G_j(\rho)=\frac14\partial_\rho
               \frac{F(\rho)-q(a_j)^2}{\rho-a_j}.
    \]
    The physical Wronskian is $u_jv_j[(\log v_j)'-(\log u_j)']$;
    its target is $-\beta_jq(a_j)$.
    """)
    return


@app.cell
def _(run_verification):
    sigma_report = run_verification(dps=35)
    if not sigma_report["all_passed"]:
        _failures = [(entry["name"], row) for entry in sigma_report["examples"]
                     for row in entry["checks"] if not row["passed"]]
        raise AssertionError(f"Sigma reconstruction failed: {_failures}")
    return (sigma_report,)


@app.cell
def _(mo, sigma_report):
    mo.vstack([
        mo.md(f"## 4. Executed results\n\n**{sigma_report['check_count']} numerical checks passed**, with exact curve-normalisation checks for each model. Working precision is 35 digits; the acceptance threshold is $10^{{-23}}$, not a claim of 35-digit accuracy."),
        mo.ui.table([
            dict(example=entry["name"], genus=entry["genus"], weights=entry["beta"],
                 scale=entry["scale"], checks=entry["check_count"],
                 maximum_scaled_error=entry["maximum_scaled_error"], seconds=entry["seconds"])
            for entry in sigma_report["examples"]
        ], selection=None),
    ])
    return


@app.cell
def _(mo, sigma_report):
    _sections = []
    for _entry in sigma_report["examples"]:
        _sections.append(mo.md(f"### `{_entry['name']}`\n\nOriginal $F$: `{_entry['original_F']}`\n\nNormalised $\\widehat F$: `{_entry['normalized_F']}`"))
        _sections.append(mo.ui.table(_entry["checks"], selection=None, page_size=12))
    mo.vstack(_sections)
    return


@app.cell
def _(example_dir, json, sigma_report):
    results_path = example_dir / "results.json"
    results_path.write_text(json.dumps(sigma_report, indent=2)+"\n")
    return (results_path,)


@app.cell
def _(mo, results_path):
    mo.md(r"""
    ## 5. Coverage and remaining questions

    Each curve is checked at two noninitial coordinates on both sheets.
    Tests cover the scalar sigma inversion, one-point stratum, Abel-lift
    continuation, fixed-point second-kind coefficients, individual modes
    versus Abel quadrature, independently reconstructed paired powers,
    canonical equations, Wronskians, analytic second derivatives, weighted
    products, and the Hamiltonian.

    The first notebook's full canonical and second-order RK4 refinement
    comparisons remain independent evidence for these same quadratures.
    This notebook adds sigma evaluation, not another fitted reconstruction.

    These are short regular paths, not turning-point or global continuation
    tests. The nonconstant example has $\deg q=5>g+1$; its sigma mode checks
    do not prove a new high-degree Onishi determinant identity. Branch-point
    limits, global multipliers, and a common spectral operator remain separate
    work. No LaTeX manuscript has been started.
    """+f"\n\nResults saved to `{results_path.name}`.")
    return


if __name__ == "__main__":
    app.run()
