"""Branch fixed-point checks for the main paper's sigma-mode formulas.

Run as a package or directly with the project's Python interpreter. Modes
are evaluated from their own sigma ratios, including at a turning point.
"""

import argparse
import json
from pathlib import Path
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "code"))

import sympy as sp
from hyperelliptic_scalar_first.special_q.model import scalar_data, canonical_field
from hyperelliptic_scalar_first.branch_points.model import BranchSigmaModes, example, quotient_jet


def verify(genus, dps=35, progress=print):
    started = perf_counter()
    data = example(genus)
    model = BranchSigmaModes(data, dps=dps, progress=progress)
    c, g, n = model.ctx, genus, model.n
    rows = []
    tolerance = c.mpf(10)**(-(dps-12))
    def record(name, actual, expected):
        error = abs(actual-expected)/max(c.one, abs(actual), abs(expected))
        rows.append(dict(check=name, scaled_error=c.nstr(error, 12), passed=bool(error < tolerance)))

    rho = data["rho"]
    potentials = []
    for j, (a, eta, h, integral) in enumerate(zip(model.a, model.eta, model.h, model.second_at_fixed)):
        jet = model.jet(eta)
        record(f"fixed {j+1}: sigma stratum", jet[(0,)*g], 0)
        if g == 3:
            record(f"fixed {j+1}: sigma_3 stratum", jet[(0, 0, 1)], 0)
        for k in range(g):
            record(f"fixed {j+1}: h_{k+1} sigma versus second-kind quadrature", h[k], -integral[k])
        quotient = sp.cancel(data["F"]/(rho-data["a"][j]))
        assert quotient.is_polynomial(rho)
        potentials.append(model.compile(rho, sp.diff(quotient, rho)/4))

    arguments, field, H = canonical_field(data)
    field_fn, H_fn = model.compile(arguments, field), model.compile(arguments, H)
    # General-Q branch example: Q=I1-I2-(gamma1-gamma2) vanishes on
    # this leaf but has partials (1,-1,0,...); B=(-1,1,0,...).
    general_Q = data["variables"][0]-data["variables"][1]+1
    general = scalar_data(general_Q, data["variables"], data["beta"], data["gamma"], rho)
    assert general["q"] == 0 and general["F"] == data["F"]
    general_args, general_field, general_H = canonical_field(general)
    general_field_fn = model.compile(general_args, general_field)
    general_H_fn = model.compile(general_args, general_H)

    for sheet in (1, -1):
        progress(f"genus {g}: regular points, sheet {sheet}")
        start = model.evaluate(model.rho_start, sheet)
        initial = model.initial(start)
        for fraction in (c.mpf("0.5"), c.one):
            x = model.rho_start+fraction*(model.rho_end-model.rho_start)
            w = start["w"]+model.abel_increment(model.rho_start, x, sheet)
            end = model.evaluate(x, sheet, w=w)
            state = model.solution(start, end, initial)
            reference = model.quadrature_solution(model.rho_start, x, initial, sheet)
            field_values = field_fn(*state)
            tag = f"sheet {sheet}, fraction {fraction}"
            record(tag+": Hamiltonian", H_fn(*state), 0)
            for j in range(n):
                u, v = state[j], state[n+j]
                du, dv = end["rates_u"][j], end["rates_v"][j]
                for label, actual, expected in (
                    ("u sigma versus square root", u, reference[j]),
                    ("v sigma versus square root", v, reference[n+j]),
                    ("product", u*v, model.a[j]-x),
                    ("raw sigma product", end["raw_u"][j]*end["raw_v"][j], model.a[j]-x),
                    ("u canonical", u*du, field_values[j]),
                    ("v canonical", v*dv, field_values[n+j]),
                    ("Wronskian", u*v*(dv-du), 0),
                    ("u polynomial potential", end["potentials_u"][j], potentials[j](x)),
                    ("v polynomial potential", end["potentials_v"][j], potentials[j](x)),
                ):
                    record(tag+f", mode {j+1}: "+label, actual, expected)
            duration = (w[0]-start["w"][0])/model.scale
            B = [model.number(value) for value in general["B"]]
            general_state = [state[j]*c.exp(B[j]*duration) for j in range(n)]
            general_state += [state[n+j]*c.exp(-B[j]*duration) for j in range(n)]
            values = general_field_fn(*general_state)
            record(tag+": general-Q Hamiltonian", general_H_fn(*general_state), 0)
            for j in range(n):
                record(tag+f": general-Q u_{j+1} canonical", general_state[j]*(end["rates_u"][j]+B[j]), values[j])
                record(tag+f": general-Q v_{j+1} canonical", general_state[n+j]*(end["rates_v"][j]-B[j]), values[n+j])

    progress(f"genus {g}: analytic continuation through the largest branch point")
    # rho=a+t^2, Y=t*sqrt(F(a+t^2)/t^2). Signed t continues the sheet.
    a, eta = model.a[-1], model.eta[-1]
    symbolic_a = data["a"][-1]
    divided = sp.cancel(model.normalization["F"]/(rho-symbolic_a))
    R = model.compile(rho, divided)
    def local_point(t):
        x = a+t*t
        Y = t*c.sqrt(R(x))
        increment = c.matrix([c.quad(lambda q: 2*(a+q*q)**k/c.sqrt(R(a+q*q)), [0, t]) for k in range(g)])
        return x, Y, eta+increment
    t0 = c.mpf("0.15")
    x0, Y0, w0 = local_point(t0)
    base = model.evaluate(x0, w=w0)
    initial = model.initial(base)
    base_u = [quotient_jet(model, w0, x0, Y0, e, norm, h, 1) for e, norm, h in zip(model.eta, model.normals, model.h)]
    base_v = [quotient_jet(model, w0, x0, Y0, e, norm, h, -1) for e, norm, h in zip(model.eta, model.normals, model.h)]
    turning_samples = []
    for t in (t0, t0/2, c.zero, -t0/2, -t0):
        x, Y, w = local_point(t)
        uj = [quotient_jet(model, w, x, Y, e, norm, h, 1) for e, norm, h in zip(model.eta, model.normals, model.h)]
        vj = [quotient_jet(model, w, x, Y, e, norm, h, -1) for e, norm, h in zip(model.eta, model.normals, model.h)]
        u = [tuple(initial[j]*v/base_u[j][0] for v in uj[j]) for j in range(n)]
        v = [tuple(initial[n+j]*val/base_v[j][0] for val in vj[j]) for j in range(n)]
        state = [vals[0] for vals in u]+[vals[0] for vals in v]
        values = field_fn(*state)
        record(f"turn t={t}: Hamiltonian", H_fn(*state), 0)
        record(f"turn t={t}: U", c.fprod(state[:n]), Y/2)
        record(f"turn t={t}: V", c.fprod(state[n:]), -Y/2)
        for j in range(n):
            tag = f"turn t={t}, mode {j+1}"
            record(tag+": product", u[j][0]*v[j][0], model.a[j]-x)
            record(tag+": u canonical derivative", u[j][1], values[j])
            record(tag+": v canonical derivative", v[j][1], values[n+j])
            record(tag+": Wronskian", u[j][0]*v[j][1]-v[j][0]*u[j][1], 0)
            record(tag+": u polynomial ODE", u[j][2], potentials[j](x)*u[j][0])
            record(tag+": v polynomial ODE", v[j][2], potentials[j](x)*v[j][0])
        record(f"turn t={t}: vanishing mode signed continuation", u[-1][0], initial[n-1]*t/t0)
        turning_samples.append(dict(t=c.nstr(t, 12), rho=c.nstr(x, 18), Y=c.nstr(Y, 18),
                                    u_last=c.nstr(u[-1][0], 18), du_last=c.nstr(u[-1][1], 18)))

    maximum = max(c.mpf(row["scaled_error"]) for row in rows)
    return dict(genus=g, dps=dps, F=str(data["F"]),
                fixed_sigma_sharp_moduli=[c.nstr(abs(v), 12) for v in model.normals],
                potential_last=str(sp.diff(sp.cancel(data["F"]/(rho-symbolic_a)), rho)/4),
                h_last=[c.nstr(v, 18) for v in model.h[-1]],
                tolerance=c.nstr(tolerance, 8), maximum_scaled_error=c.nstr(maximum, 12),
                check_count=len(rows), all_passed=all(row["passed"] for row in rows),
                checks=rows, turning_samples=turning_samples,
                seconds=round(perf_counter()-started, 3))


def exact_checks():
    x, a = sp.symbols("x a")
    checks = {}
    for g in (2, 3):
        data = example(g)
        checks[f"genus {g} squarefree"] = data["squarefree"]
        checks[f"genus {g} all fixed points branch"] = all(b == 0 for b in data["ordinate"])
        checks[f"genus {g} all corrections vanish"] = all(b == 0 for b in data["B"])
    rho = sp.Symbol("rho")
    variables = sp.symbols("I1:5")
    repeated = scalar_data(0, variables, (2, 1, 1, 1), (-3, -1, 1, 3), rho)
    checks["repeated-weight zero energy curve singular"] = not repeated["squarefree"]
    checks["double-weight fixed point has F'=0"] = sp.diff(repeated["F"], rho).subs(rho, repeated["a"][0]) == 0
    # For any simple branch point, the square-root solution obeys the
    # divided-difference polynomial potential, using rho'^2=F, rho''=F'/2.
    P = sum(sp.Symbol(f"p{k}")*(x-a)**k for k in range(5))
    F = (x-a)*P
    derived = sp.diff(F, x)/(4*(x-a))-F/(4*(x-a)**2)
    checks["branch square-root ODE identity"] = sp.cancel(derived-sp.diff(P, x)/4) == 0
    return checks


def run_verification(dps=35, genera=(2, 3), progress=print):
    exact = exact_checks()
    reports = [verify(g, dps, progress) for g in genera]
    return dict(scope="Simple branch fixed points on smooth odd-degree curves; explicit genera 2 and 3; no global monodromy or determinant extension asserted.",
                exact_checks=exact, examples=reports,
                check_count=len(exact)+sum(r["check_count"] for r in reports),
                all_passed=all(exact.values()) and all(r["all_passed"] for r in reports))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=35)
    parser.add_argument("--genus", type=int, choices=(2, 3), action="append")
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("results.json"))
    args = parser.parse_args()
    report = run_verification(args.dps, args.genus or (2, 3), lambda s: print(s, flush=True))
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(f"{report['check_count']} checks; all passed={report['all_passed']}")
    for example_report in report["examples"]:
        print(f"genus {example_report['genus']}: max error {example_report['maximum_scaled_error']}")
        for row in example_report["checks"]:
            if not row["passed"]:
                print(row)
    if not report["all_passed"]:
        raise SystemExit(1)
