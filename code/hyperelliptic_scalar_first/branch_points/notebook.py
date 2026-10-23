import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import sys
    from pathlib import Path
    directory = Path(__file__).resolve().parent
    if str(directory.parents[1]) not in sys.path:
        sys.path.insert(0, str(directory.parents[1]))
    import json
    import marimo as mo
    from hyperelliptic_scalar_first.branch_points.verification import run_verification
    return directory, json, mo, run_verification


@app.cell
def _(mo):
    mo.md(r"""
    # Simple branch fixed points: shifted sigma modes

    Both modes are evaluated independently from their own shifted sigma
    numerators and natural exponentials. Square roots supply a reference;
    they do not replace the sigma representation.

    We test unit weights and $Q=0$ on two smooth curves:
    \[
    F(\rho)=4\prod_{a=-g}^{g}(\rho-a),\qquad g=2,3.
    \]
    All fixed mode points are simple branch points. On a common contour,
    the endpoint substitution $\rho=a-\delta t^2$ removes the square-root
    singularity from both kinds of Abel integrals. Independently integrated
    second-kind values are compared with the paper's sigma formulas for $h$.
    The fixed sigma-sharp normalizers are reported below.

    At regular moving points on both sheets we check sigma modes against
    square-root ratios, paired products, canonical equations, zero Wronskians
    and polynomial potentials. A second Hamiltonian, $Q=I_1-I_2+1$, vanishes
    on this same power line but gives $B=(-1,1,0,\ldots,0)$; it tests the
    canonical equations with nonzero polynomial corrections.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Through the turning point

    At the largest branch point $a=g$, use the signed local parameter
    \[
    \rho=a+t^2,\quad Y=t\sqrt{F(a+t^2)/t^2},\quad
    \frac{dw_k}{dt}=\frac{2(a+t^2)^{k-1}}{\sqrt{F(a+t^2)/t^2}}.
    \]
    The Abel lift is continued through $t=0$ onto the opposite sheet.
    Sigma jets give the modes and their first two propagation derivatives
    directly, without dividing by the vanishing numerator. We check
    the canonical equations, $U=Y/2$, $V=-Y/2$, paired products and
    \[
    \psi''=G_a(\rho)\psi,\qquad
    G_a(\rho)=\frac14\partial_\rho\frac{F(\rho)}{\rho-a}.
    \]
    The attached mode must continue proportionally to signed $t$, not $|t|$.
    All these checks include the turning point itself.
    """)
    return


@app.cell
def _(run_verification):
    report = run_verification(dps=35)
    assert report["all_passed"], "Inspect the detailed verification rows."
    return (report,)


@app.cell
def _(directory, json, report):
    results_path = directory / "results.json"
    results_path.write_text(json.dumps(report, indent=2)+"\n")
    return


@app.cell
def _(mo, report):
    mo.vstack([
        mo.md(f"## Executed results\n\n**{report['check_count']} checks passed**, at 35 working digits with numerical acceptance threshold $10^{{-23}}$."),
        mo.ui.table([dict(genus=e["genus"], checks=e["check_count"], maximum_scaled_error=e["maximum_scaled_error"], seconds=e["seconds"]) for e in report["examples"]], selection=None),
        mo.ui.table([dict(check=k, passed=v) for k, v in report["exact_checks"].items()], selection=None),
    ])
    return


@app.cell
def _(mo, report):
    sections = []
    for entry in report["examples"]:
        sections.extend([
            mo.md(f"### Genus {entry['genus']}\n\nCurve: `{entry['F']}`\n\nPotential for the largest branch point: `{entry['potential_last']}`\n\nFixed sigma-sharp moduli: `{entry['fixed_sigma_sharp_moduli']}`"),
            mo.ui.table(entry["turning_samples"], selection=None),
            mo.ui.table(entry["checks"], selection=None, page_size=12),
        ])
    mo.vstack(sections)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Scope

    These computations cover every fixed branch point on these two curves,
    regular moving points on both sheets, and continuation through the largest
    branch point on each curve. The arbitrary-genus
    justification is in `proof.md`. These numerical checks do not test
    global multipliers; `determinant_notebook.py` checks the determinants.

    For special $Q$, the two canonical modes are proportional at a branch
    fixed point. The polynomial potential remains, but this pair is not
    fundamental. The exact checks also confirm the singularity caused by
    repeated-weight branch fixed points. We have not constructed a second
    independent solution or treated singular curves here.
    """)
    return


if __name__ == "__main__":
    app.run()
