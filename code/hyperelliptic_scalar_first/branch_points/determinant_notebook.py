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
    from hyperelliptic_scalar_first.branch_points.determinants import run_verification
    return directory, json, mo, run_verification


@app.cell
def _(mo):
    mo.md(r"""
    # Branch fixed points: arbitrary genus and Onishi determinants

    The arbitrary-genus argument is recorded in `proof.md`, with precise
    source references. Its main points are:

    1. Onishi's Proposition 6.6(1) gives nonzero sigma-sharp at every finite
       curve point. The Abel lift and second-kind primitive are holomorphic
       in a local branch parameter. Thus the sigma kernel and its logarithmic
       identities extend from nonbranch fixed points by analytic continuation.
    2. At a branch fixed point, write $x-a=t^2$. Then both $q-y$ and $q+y$
       have simple zeros in $t$. This replaces the nonbranch step that assumes
       the other factor is nonzero; the weighted zero divisor is unchanged.
    3. The determinant column operation and auxiliary Vandermonde splitting
       still apply. Although a branch point equals its geometric involution,
       its opposite lift must be **minus** its Abel vector. Half-period
       quasi-period factors must be retained.

    The computational tests below check these consequences rather than
    numerically claiming a proof for every genus.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Checks performed

    Exact symbolic determinants test both main-paper column families, with
    mixed branch and nonbranch rows, in genera two through six. Smoothness,
    the branch row, interpolation, and opposite-sheet interpolation are checked.

    Numerical tests on the two earlier genus-two/genus-three curves evaluate
    the **full** Onishi sigma determinant formula, using its corrected sign.
    The main and opposite-sheet identities, residual auxiliary identity,
    composite-factor splitting, and complete raw sigma products are checked
    on both sheets. A further matrix mixes branch and regular fixed rows.

    We compute the integer coordinates of the total Abel period and move
    one unit-weight point by minus this period. Its second-kind coefficient
    is shifted using independently computed second-kind periods and checked
    against the sigma derivatives at the shifted point. The corrected modes
    must remain unchanged. With balanced lifts, the sum of $h_j$ vanishes
    and the raw products reproduce the determinants up to a fixed constant.
    """)
    return


@app.cell
def _(run_verification):
    report = run_verification(dps=35)
    assert report["all_passed"], "Inspect the verification rows."
    return (report,)


@app.cell
def _(directory, json, report):
    path = directory / "determinant_results.json"
    path.write_text(json.dumps(report, indent=2)+"\n")
    return


@app.cell
def _(mo, report):
    mo.vstack([
        mo.md(f"## Executed results\n\n**{report['check_count']} checks passed**, at 35 working digits and numerical threshold $10^{{-23}}$."),
        mo.ui.table([dict(genus=e["genus"], checks=e["check_count"], maximum_scaled_error=e["maximum_scaled_error"], period_coordinates=e["integer_period_coordinates"], unbalanced_raw_product_defect=e["unbalanced_raw_product_defect"]) for e in report["examples"]], selection=None),
        mo.ui.table([dict(check=k, passed=v) for k,v in report["exact_checks"].items()], selection=None, page_size=12),
    ])
    return


@app.cell
def _(mo, report):
    sections = []
    for entry in report["examples"]:
        sections.extend([
            mo.md(f"### Genus {entry['genus']}\n\nUnbalanced sum of h: `{entry['unbalanced_sum_h']}`"),
            mo.ui.table(entry["checks"], selection=None, page_size=12),
        ])
    mo.vstack(sections)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Interpretation and scope

    Branch rows do not force a confluent determinant: their coordinates
    remain distinct. At zero energy all fixed ordinates vanish, so the two
    algebraic main determinants coincide. This is consistent with
    $U=y/2$, $V=-y/2$; it does not force either determinant to vanish.

    Unbalanced raw products generally retain an exponential multiplier.
    The compatible period shift removes it from the complete raw product;
    individual corrected modes are invariant under that change of lift.

    These checks cover the two short column families with internally distinct
    main/residual lists. Longer column lists are tested in `../higher_degree/`.
    This notebook does not address repeated rows, singular or even-degree
    curves, or independence of the two canonical modes at a branch point.
    """)
    return


if __name__ == "__main__":
    app.run()
