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
    from hyperelliptic_scalar_first.higher_degree.verification import run_verification
    return directory, json, mo, run_verification


@app.cell
def _(directory, mo):
    mo.md((directory / "proof.md").read_text())
    return


@app.cell
def _(run_verification):
    report = run_verification(dps=40)
    assert report["all_passed"], "Inspect the detailed verification rows."
    return (report,)


@app.cell
def _(directory, json, report):
    (directory / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    return


@app.cell
def _(mo, report):
    mo.vstack([
        mo.md(f"## Executed results\n\n**{report['check_count']} checks passed.** Numerical checks use 40 working digits and a scaled-error threshold of $10^{{-26}}$."),
        mo.ui.table([dict(example=e["name"], genus=e["genus"], modes=e["N"], degree_q=e["degree_q"], degrees_B=e["degrees_B"], maximum_scaled_error=e["maximum_scaled_error"], seconds=e["seconds"]) for e in report["examples"]], selection=None),
        mo.ui.table([dict(check=k, passed=v) for k, v in report["exact_checks"].items()], selection=None, page_size=12),
    ])
    return


@app.cell
def _(mo, report):
    sections = []
    for entry in report["examples"]:
        sections.extend([
            mo.md(f"### {entry['name']}\n\nRestricted polynomial: `{entry['q']}`\n\nScalar polynomial: `{entry['F']}`\n\nFull multivariate polynomial: `{entry['Q']}`"),
            mo.ui.table([dict(mode=j + 1, B=b, **r) for j, (b, r) in enumerate(zip(entry["B"], entry["reductions"]))], selection=None),
            mo.ui.table([dict(check=k, passed=v) for k, v in entry["exact_checks"].items()], selection=None),
            mo.md(f"Longer determinant metadata: `{entry['determinant']}`"),
            mo.ui.table(entry["checks"], selection=None, page_size=15),
        ])
    mo.vstack(sections)
    return


if __name__ == "__main__":
    app.run()
