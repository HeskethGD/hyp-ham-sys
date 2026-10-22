"""Sigma versus Abel quadrature in original Hamiltonian units.

PYTHONPATH=code .venv/bin/python -m hyperelliptic_scalar_first.sigma_modes.verification
"""

import argparse
import importlib.metadata
import json
from pathlib import Path
from time import perf_counter

import sympy as sp

from ..special_q.model import canonical_field, example_models
from .model import SigmaModes


def verify_model(data, dps=35, progress=None):
    started = perf_counter()
    model = SigmaModes(data, dps=dps, progress=progress)
    c, n = model.ctx, model.n
    tolerance = c.mpf(10)**(-(dps-12))
    rows = []
    samples = []

    def record(name, actual, expected):
        error = abs(actual-expected)/max(c.one, abs(actual), abs(expected))
        rows.append(dict(check=name, scaled_error=c.nstr(error, 12), passed=bool(error < tolerance)))

    for j, (a, b, h) in enumerate(zip(model.a, model.b, model.h)):
        fixed_jet = model.jet(model.eta[j])
        record(f"mode {j+1} fixed-point sigma divisor", fixed_jet[(0,)*model.genus], 0)
        if model.genus == 3:
            record(f"mode {j+1} fixed-point sigma_3 stratum", fixed_jet[(0, 0, 1)], 0)
        integral = model.second_at_fixed[j]
        for k in range(model.genus):
            record(f"mode {j+1} h_{k+1}: sigma versus prescribed-path second-kind integral", h[k], -integral[k])

    arguments, derivatives, Hamiltonian = canonical_field(data)
    field = model.compile(arguments, derivatives, cse=True)
    energy = model.compile(arguments, Hamiltonian, cse=True)
    rho = data["rho"]
    potential_functions = []
    for a, b in zip(data["a"], data["ordinate"]):
        quotient, remainder = sp.div(data["F"]-b*b, rho-a, rho)
        if remainder != 0:
            raise AssertionError("Fixed ordinate must satisfy b^2=F(a)")
        potential_functions.append(model.compile(rho, sp.diff(quotient, rho)/4))

    for sheet in (1, -1):
        if progress:
            progress(f"checking sigma modes on sheet {sheet}")
        start = model.evaluate(model.rho_start, sheet)
        initial = model.initial(start)
        for fraction in (c.mpf("0.5"), c.one):
            coordinate = model.rho_start+fraction*(model.rho_end-model.rho_start)
            increment = model.abel_increment(model.rho_start, coordinate, sheet)
            # Continue the moving lift by quadrature; do not independently
            # reduce each endpoint modulo periods or fit any exponent.
            end = model.evaluate(coordinate, sheet, w=start["w"]+increment)
            direct_lift = model.curve.abel_map_kind_1((coordinate, end["Y"]), reduce=False).value
            tag = f"sheet {sheet}, rho={c.nstr(coordinate, 8)}"
            for k in range(model.genus):
                record(tag+f" continued versus direct Abel lift {k+1}", end["w"][k], direct_lift[k])
            duration = c.quad(lambda x: 1/(sheet*c.sqrt(model.F_original(x))), [model.rho_start, coordinate])
            record(tag+" physical propagation normalization", increment[0], model.scale*duration)
            record(tag+" sigma coordinate inversion", end["inversion"], coordinate)
            record(tag+" sigma one-point divisor", end["sigma"], 0)
            if model.genus == 3:
                record(tag+" sigma_3 one-point stratum", end["sigma3"], 0)
            state = model.solution(start, end, initial)
            reference = model.quadrature_solution(model.rho_start, coordinate, initial, sheet)
            canonical_rates = field(*state)
            for j in range(n):
                u, v = state[j], state[n+j]
                du, dv = end["rates_u"][j], end["rates_v"][j]
                power = data["beta"][j]*(model.a[j]-coordinate)
                bare_u = (model.b_original[j]-end["y"])/(2*(model.a[j]-coordinate))
                bare_v = (-model.b_original[j]-end["y"])/(2*(model.a[j]-coordinate))
                G = potential_functions[j](coordinate)
                record(tag+f" u_{j+1} sigma versus quadrature", u, reference[j])
                record(tag+f" v_{j+1} sigma versus quadrature", v, reference[n+j])
                record(tag+f" mode {j+1} raw sigma product", end["raw_u"][j]*end["raw_v"][j], power)
                record(tag+f" mode {j+1} independently reconstructed product", u*v, power)
                record(tag+f" u_{j+1} physical log derivative", du, bare_u)
                record(tag+f" v_{j+1} physical log derivative", dv, bare_v)
                record(tag+f" u_{j+1} canonical equation", u*du, canonical_rates[j])
                record(tag+f" v_{j+1} canonical equation", v*dv, canonical_rates[n+j])
                record(tag+f" mode {j+1} physical Wronskian", u*v*(dv-du), -data["beta"][j]*model.b_original[j])
                record(tag+f" u_{j+1} sigma second derivative potential", end["potentials_u"][j], G)
                record(tag+f" v_{j+1} sigma second derivative potential", end["potentials_v"][j], G)
            U = c.fprod(value**weight for value, weight in zip(state[:n], data["beta"]))
            V = c.fprod(value**weight for value, weight in zip(state[n:], data["beta"]))
            record(tag+" weighted U", U, (model.q(coordinate)+end["y"])/2)
            record(tag+" weighted V", V, (model.q(coordinate)-end["y"])/2)
            record(tag+" Hamiltonian", energy(*state), 0)
            samples.append(dict(sheet=sheet, rho=c.nstr(coordinate, 16), physical_z=c.nstr(duration, 16),
                                modes=[dict(real=c.nstr(value.real, 18), imag=c.nstr(value.imag, 18)) for value in state]))
    maximum = max(c.mpf(row["scaled_error"]) for row in rows)
    return dict(name=data["name"], dps=dps, genus=model.genus, beta=list(data["beta"]),
                original_F=str(data["F"]), normalized_F=str(model.normalization["F"]),
                scale=str(model.normalization["scale"]),
                convention="t=s*z, Y=y/s, dw_k/dz=s*rho**(k-1); original canonical amplitudes unchanged",
                exact_normalization_checks={name: bool(value) for name, value in model.normalization["checks"].items()},
                fixed_paths=model.fixed_paths,
                tolerance=c.nstr(tolerance, 6), check_count=len(rows), maximum_scaled_error=c.nstr(maximum, 12),
                all_passed=all(row["passed"] for row in rows), checks=rows, samples=samples,
                seconds=round(perf_counter()-started, 3))


def run_verification(dps=35, selected=None, progress=None):
    if dps < 30:
        raise ValueError("Use at least 30 working digits")
    reports = []
    for data in example_models():
        if selected is not None and data["name"] not in selected:
            continue
        callback = (lambda message, name=data["name"]: progress(f"{name}: {message}")) if progress else None
        report = verify_model(data, dps=dps, progress=callback)
        reports.append(report)
        if progress:
            progress(f"{data['name']}: {report['check_count']} checks; max error {report['maximum_scaled_error']}; passed={report['all_passed']}")
    if not reports:
        raise ValueError("No example names matched")
    return dict(environment={name: importlib.metadata.version(name) for name in ("generapy", "sympy", "mpmath", "marimo")},
                scope="Short paths on both sheets; independently reconstructed u and v; analytic sigma-jet first and second derivatives. No new determinant or global monodromy claims.",
                all_passed=all(r["all_passed"] for r in reports),
                check_count=sum(r["check_count"] for r in reports), examples=reports)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=35)
    parser.add_argument("--example", action="append", dest="selected")
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("results.json"))
    args = parser.parse_args()
    report = run_verification(args.dps, args.selected, progress=lambda message: print(message, flush=True))
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(f"Wrote {args.output}")
    if not report["all_passed"]:
        failed = [dict(example=r["name"], **row) for r in report["examples"] for row in r["checks"] if not row["passed"]]
        raise AssertionError(f"Sigma checks failed: {failed}")


if __name__ == "__main__":
    main()
