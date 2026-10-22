"""Exact identities, independent coefficient searches and short ODE checks.

Run from the repository root:
    PYTHONPATH=code .venv/bin/python -m hyperelliptic_scalar_first.special_q.verification
"""

import argparse
import importlib.metadata
import json
from pathlib import Path

import sympy as sp
from sympy.printing.pycode import MpmathPrinter
from mpmath import mp

from .model import (cancellation_residuals, canonical_field, example_models,
                    homogeneous_kernel, scalar_data, special_polynomial)


def exact_checks():
    checks = {}
    classification = []
    patterns = [(1,), (2,), (3,), (4,), (1, 1), (1, 2), (1, 3),
                (2, 2), (1, 1, 1), (2, 1, 1), (1, 1, 1, 1), (3, 1, 1, 1)]
    for beta in patterns:
        for degree in range(5):
            variables, monomials, kernel = homogeneous_kernel(beta, degree)
            expected = int(degree == 0 or 2*degree == sum(beta))
            key = f"coefficient kernel {beta}, degree {degree}"
            checks[key] = len(kernel) == expected
            if kernel and degree:
                candidate = sum(v*m for v, m in zip(kernel[0], monomials))
                diagonal = candidate.subs(dict(zip(variables, beta)), simultaneous=True)
                H = special_polynomial(variables, beta)
                checks[key + " agrees with H"] = diagonal != 0 and sp.expand(
                    candidate/diagonal - H/(-1)**degree) == 0
            classification.append(dict(beta=list(beta), degree=degree,
                                       dimension=len(kernel), expected=expected))

    for beta in [(1,)*6, (1,)*8, (2, 1, 1, 1, 1), (2, 1, 1, 1, 1, 1, 1)]:
        variables = sp.symbols(f"I1:{len(beta)+1}")
        H = special_polynomial(variables, beta)
        checks[f"full weighted identities {beta}"] = all(
            residual == 0 for residual in cancellation_residuals(H, variables, beta))

    rho, a, b, c, y, beta_symbol = sp.symbols("rho a b c y beta", nonzero=True)
    F = sp.Function("F")(rho)
    B = sp.Function("B")(rho)
    d = a-rho
    L = (b-y)/(2*d)+B
    R = (-b-y)/(2*d)-B
    tangent = lambda expression: y*sp.diff(expression, rho) + sp.diff(F, rho)*sp.diff(expression, y)/2

    def curve_zero(expression):
        numerator = sp.together(expression).as_numer_denom()[0]
        return sp.Poly(numerator, y).rem(sp.Poly(y*y-F, y)).is_zero

    G_b = sp.diff((F-b*b)/(rho-a), rho)/4
    even_part = G_b+B*B+b*B/d
    odd_part = y*(sp.diff(B, rho)-B/d)
    checks["general u second derivative"] = curve_zero(tangent(L)+L*L-even_part-odd_part)
    checks["general v second derivative"] = curve_zero(tangent(R)+R*R-even_part+odd_part)
    W = beta_symbol*d*(R-L)
    checks["weighted Wronskian"] = sp.simplify(W+beta_symbol*(b+2*d*B)) == 0
    checks["Wronskian derivative"] = sp.simplify(
        tangent(W)+2*beta_symbol*y*(d*sp.diff(B, rho)-B)) == 0
    L_c = (-c-y)/(2*d)
    R_c = (c-y)/(2*d)
    G_c = sp.diff((F-c*c)/(rho-a), rho)/4
    checks["normalised Wronskian convention"] = sp.simplify(d*(R_c-L_c)-c) == 0
    checks["constant W u potential"] = curve_zero(tangent(L_c)+L_c*L_c-G_c)
    checks["constant W v potential"] = curve_zero(tangent(R_c)+R_c*R_c-G_c)

    g2, g3 = sp.symbols("g2 g3")
    cubic = 4*rho**3-g2*rho-g3
    lame = sp.diff(sp.cancel((cubic-cubic.subs(rho, a))/(rho-a)), rho)/4
    checks["Weierstrass Lame recovery"] = sp.expand(lame-2*rho-a) == 0
    for degree in range(2, 9):
        polynomial = sum((k+1)*rho**k for k in range(degree+1))
        at_a = polynomial.subs(rho, a)
        quotient, remainder = sp.div(polynomial-c*c, rho-a, rho)
        potential = sp.diff((polynomial-c*c)/(rho-a), rho)/4
        checks[f"pole obstruction degree {degree}"] = sp.cancel(
            potential-sp.diff(quotient, rho)/4-(c*c-at_a)/(4*(rho-a)**2)) == 0
        checks[f"potential degree {degree}"] = sp.degree(sp.diff(quotient, rho), rho) == degree-2
        checks[f"division remainder degree {degree}"] = sp.expand(remainder-at_a+c*c) == 0

    variables = sp.symbols("I1:4")
    ell = sp.symbols("ell1:4")
    linear = sum(k*x for k, x in zip(ell, variables))
    linear_data = scalar_data(linear, variables, (1, 1, 1), (-1, 0, 1), rho)
    checks["three-mode linear Q has constant shifts"] = all(
        sp.expand(value-sum(ell)/2+coefficient) == 0
        for value, coefficient in zip(linear_data["B"], ell))
    leaf = variables[0]-variables[1]+1
    invisible = scalar_data(1+leaf**2*(1+variables[2]), variables,
                            (1, 1, 1), (-1, 0, 1), rho)
    checks["fixed-leaf nonuniqueness"] = invisible["q"] == 1 and all(v == 0 for v in invisible["B"])
    checks["fixed-leaf example is not universal"] = any(
        v != 0 for v in cancellation_residuals(invisible["Q"], variables, (1, 1, 1)))
    branch_variables = sp.symbols("J1:6")
    branch_data = scalar_data(0, branch_variables, (1,)*5, (-2, -1, 0, 1, 2), rho)
    checks["zero-energy curve can remain smooth"] = branch_data["squarefree"] and branch_data["genus"] == 2
    checks["zero-energy fixed points are branch points"] = all(value == 0 for value in branch_data["ordinate"])

    examples = []
    for data in example_models():
        r = data["rho"]
        name = data["name"]
        checks[name+" smooth odd curve"] = data["squarefree"] and data["degree"] % 2 == 1
        checks[name+" expected genus"] = data["genus"] == (3 if name.endswith("g3") else 2)
        checks[name+" special corrections"] = all(v == 0 for v in data["B"])
        checks[name+" distinct nonbranch points"] = len(set(data["a"])) == len(data["a"]) and all(
            v != 0 for v in data["ordinate"])
        potentials = []
        for j, (point, ordinate) in enumerate(zip(data["a"], data["ordinate"])):
            quotient, remainder = sp.div(data["F"]-ordinate**2, r-point, r)
            potential = sp.diff(quotient, r)/4
            checks[f"{name} polynomial potential {j+1}"] = remainder == 0 and sp.degree(potential, r) == data["degree"]-2
            potentials.append(str(sp.expand(potential)))
        if sp.degree(data["q"], r) <= data["genus"]:
            # The actual lower-degree Onishi columns, with the paper's y/2
            # normalization. This checks interpolation, not sigma evaluation.
            g = data["genus"]
            indices = list(range(g+1))
            fixed = data["a"][:g+1]
            delta = sp.prod(fixed[j]-fixed[i] for j in range(g+1) for i in range(j))
            rows = [[point**k for k in range(g+1)] + [data["q"]/2] for point in fixed]
            moving = [r**k for k in range(g+1)] + [y/2]
            D = sp.Matrix(rows+[moving]).det()
            opposite_rows = [row[:-1]+[-row[-1]] for row in rows]
            Dstar = sp.Matrix(opposite_rows+[moving]).det()
            checks[name+" determinant V"] = sp.expand(D-delta*(y-data["q"])/2) == 0
            checks[name+" determinant U"] = sp.expand(Dstar-delta*(y+data["q"])/2) == 0
            checks[name+" determinant energy"] = sp.expand((Dstar-D)/delta-data["q"]) == 0
            occurrences = [j for j, weight in enumerate(data["beta"]) for _ in range(weight)]
            for j in indices:
                occurrences.remove(j)
            checks[name+" distinct residual list"] = len(occurrences) == g and len(set(occurrences)) == g
        examples.append(dict(name=name, beta=list(data["beta"]), gamma=list(data["gamma"]),
                             degree_Q=int(sp.Poly(data["Q"], *data["variables"]).total_degree()),
                             degree_q=int(sp.degree(data["q"], r)), degree_F=data["degree"],
                             genus=data["genus"], squarefree=data["squarefree"],
                             q=str(data["q"]), F=str(data["F"]),
                             fixed_ordinates=list(map(str, data["ordinate"])), potentials=potentials,
                             determinant_scope="covered" if sp.degree(data["q"], r) <= data["genus"]+1 else "requires longer-column argument"))

    # Check the even-total degree dichotomy using scalar polynomial parts,
    # independently of the multivariable coefficient-extraction implementation.
    degree_cases = []
    for beta, gamma in [((1,)*6, (-5, -3, -1, 1, 3, 5)),
                        ((1,)*8, (-7, -5, -3, -1, 1, 3, 5, 7)),
                        ((5, 5), (-1, 1))]:
        m = sum(beta)//2
        t = sp.Symbol("t")
        points = [sp.Rational(g, w) for g, w in zip(gamma, beta)]
        reversed_S = sp.Poly(sp.prod((1-p*t)**w for p, w in zip(points, beta)), t)
        # Independent formal square root via h(t)^2 = reversed_S(t).
        # This avoids expanding a product of half-integer series.
        h = [sp.Integer(1)]
        for k in range(1, m+1):
            h.append((reversed_S.nth(k)-sum(h[i]*h[k-i] for i in range(1, k)))/2)
        P = sum(h[k]*rho**(m-k) for k in range(m+1))
        S = sp.prod((rho-p)**w for p, w in zip(points, beta))
        C = sp.prod(sp.Integer(w)**w for w in beta)
        checks[f"polynomial square-root remainder {beta}"] = sp.degree(sp.expand(P*P-S), rho) <= m-1
        for label, kappa, epsilon in [("generic", 1, 1), ("cancelled", 2*sp.sqrt(C), 1),
                                      ("cancelled zero offset", 2*sp.sqrt(C), 0)]:
            curve = sp.Poly(sp.expand((kappa*P+epsilon)**2-4*C*S), rho)
            bound = 2*m if label == "generic" else m if epsilon else m-1
            checks[f"degree dichotomy {beta} {label}"] = curve.degree() <= bound
            degree_cases.append(dict(beta=list(beta), case=label, scalar_degree=int(curve.degree()), bound=bound))

    # Read-only regression: the original paper examples retain their genera
    # and their nonzero corrections. No original generated files are written.
    from ..model import example, genus_three_example
    original_examples = []
    for factory in (example, genus_three_example):
        data = factory()
        degree = int(sp.degree(data["F"], data["rho"]))
        checks[f"original genus {data['genus']} unchanged"] = degree == 2*data["genus"]+1 and data["squarefree"]
        checks[f"original genus {data['genus']} not special"] = any(value != 0 for value in data["B"])
        original_examples.append(dict(genus=int(data["genus"]), total_weight=sum(data["beta"]),
                                      degree_F=degree, special=False))
    failed = [name for name, value in checks.items() if not value]
    if failed:
        raise AssertionError(f"Exact checks failed: {failed}")
    return dict(check_count=len(checks), all_passed=True, checks={k: bool(v) for k, v in checks.items()},
                coefficient_search=classification, examples=examples,
                degree_cases=degree_cases, original_examples=original_examples)


def rk4(field, initial, duration, steps):
    state = list(initial)
    h = duration/steps
    for _ in range(steps):
        k1 = field(state)
        k2 = field([x+h*k/2 for x, k in zip(state, k1)])
        k3 = field([x+h*k/2 for x, k in zip(state, k2)])
        k4 = field([x+h*k for x, k in zip(state, k3)])
        state = [x+h*(a+2*b+2*c+d)/6 for x, a, b, c, d in zip(state, k1, k2, k3, k4)]
    return state


def numerical_checks(dps=50):
    """Compare curve quadratures with independently differentiated Hamiltonians.

    This checks elementary Abel-integral solutions, not sigma evaluation.
    A second-order IVP is also integrated independently with the same curve.
    """
    if dps < 35:
        raise ValueError("Use at least 35 digits for the refinement checks")
    ctx = mp.clone()
    ctx.dps = dps
    # A cloned context alone selects Python's float-producing printer.
    # Explicit mpf literals preserve the requested precision for rationals.
    printer = MpmathPrinter({"fully_qualified_modules": False})
    compile_numeric = lambda arguments, expression, **options: sp.lambdify(
        arguments, expression, modules=ctx, printer=printer, **options)
    reports = []
    for data in example_models():
        rho = data["rho"]
        numeric = lambda expression: ctx.mpf(str(sp.N(expression, dps+5)))
        q = compile_numeric(rho, data["q"])
        F = compile_numeric(rho, data["F"])
        Fprime = compile_numeric(rho, sp.diff(data["F"], rho))
        start, end = map(numeric, (data["rho_start"], data["rho_end"]))
        duration = ctx.quad(lambda x: 1/ctx.sqrt(F(x)), [start, end])
        y0 = ctx.sqrt(F(start))
        U0 = (q(start)+y0)/2
        n = len(data["beta"])
        u0 = [ctx.mpf(1)]*n
        u0[0] = ctx.root(U0, data["beta"][0])
        initial_powers = [numeric(g)-w*start for g, w in zip(data["gamma"], data["beta"])]
        v0 = [power/u for power, u in zip(initial_powers, u0)]
        expected_u, expected_v = [], []
        for point, ordinate, u, v in zip(data["a"], data["ordinate"], u0, v0):
            a, b = numeric(point), numeric(ordinate)
            amplitude = ctx.sqrt((a-end)/(a-start))
            integral = ctx.quad(lambda x: b/(2*(a-x)*ctx.sqrt(F(x))), [start, end])
            expected_u.append(u*amplitude*ctx.exp(integral))
            expected_v.append(v*amplitude*ctx.exp(-integral))
        expected = expected_u+expected_v
        variables, derivatives, H = canonical_field(data)
        evaluate = compile_numeric(variables, derivatives, cse=True)
        energy = compile_numeric(variables, H, cse=True)
        field = lambda state: evaluate(*state)
        error = lambda actual: max(abs(x-y)/max(ctx.mpf(1), abs(y)) for x, y in zip(actual, expected))
        errors = []
        final = None
        for steps in (16, 32, 64):
            final = rk4(field, u0+v0, duration, steps)
            errors.append(error(final))
        factors = [errors[j]/errors[j+1] for j in range(2)]
        if not all(ctx.mpf(12) < ratio < ctx.mpf(20) for ratio in factors):
            raise AssertionError(f"{data['name']}: no fourth-order convergence: {factors}")
        if errors[-1] >= ctx.mpf("1e-9"):
            raise AssertionError(f"{data['name']}: mode error too large: {errors[-1]}")

        potentials = []
        for point, ordinate in zip(data["a"], data["ordinate"]):
            quotient = sp.div(data["F"]-ordinate**2, rho-point, rho)[0]
            potentials.append(compile_numeric(rho, sp.diff(quotient, rho)/4))
        du0, dv0 = [], []
        for point, ordinate, u, v in zip(data["a"], data["ordinate"], u0, v0):
            a, b = numeric(point), numeric(ordinate)
            du0.append(u*(b-y0)/(2*(a-start)))
            dv0.append(v*(-b-y0)/(2*(a-start)))

        def second_order_field(state):
            coordinate, speed = state[:2]
            us, ps = state[2:2+n], state[2+n:2+2*n]
            vs, qs = state[2+2*n:2+3*n], state[2+3*n:]
            return [speed, Fprime(coordinate)/2] + list(ps) + [G(coordinate)*u for G, u in zip(potentials, us)] + list(qs) + [G(coordinate)*v for G, v in zip(potentials, vs)]

        second = rk4(second_order_field, [start, y0]+u0+du0+v0+dv0, duration, 64)
        second_error = error(second[2:2+n]+second[2+2*n:2+3*n])
        if second_error >= ctx.mpf("1e-9"):
            raise AssertionError(f"{data['name']}: second-order solution error {second_error}")
        wronskian_error = max(abs(u*s-v*p+w*numeric(b))/max(ctx.mpf(1), abs(w*numeric(b)))
                              for u, p, v, s, w, b in zip(second[2:2+n], second[2+n:2+2*n],
                                  second[2+2*n:2+3*n], second[2+3*n:], data["beta"], data["ordinate"]))
        if wronskian_error >= ctx.mpf("1e-9"):
            raise AssertionError(f"{data['name']}: Wronskian error {wronskian_error}")
        reports.append(dict(name=data["name"], dps=dps, sheet="positive real",
                            rho_start=str(start), rho_end=str(end), duration=ctx.nstr(duration, 18),
                            steps=[16, 32, 64], mode_errors=[ctx.nstr(e, 12) for e in errors],
                            refinement_factors=[ctx.nstr(v, 10) for v in factors],
                            second_order_mode_error=ctx.nstr(second_error, 12),
                            second_order_wronskian_error=ctx.nstr(wronskian_error, 12),
                            final_hamiltonian_residual=ctx.nstr(abs(energy(*final)), 12)))
    return dict(all_passed=True, scope="Short real paths; Abel quadratures versus full canonical and second-order RK4. No sigma evaluation.", examples=reports)


def run_verification(dps=50, exact_only=False):
    report = dict(environment={name: importlib.metadata.version(name) for name in ("sympy", "mpmath", "marimo")},
                  exact=exact_checks())
    if not exact_only:
        report["numerical"] = numerical_checks(dps)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--exact-only", action="store_true")
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("results.json"))
    args = parser.parse_args()
    report = run_verification(args.dps, args.exact_only)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Passed {report['exact']['check_count']} exact checks; wrote {args.output}")
    for example in report.get("numerical", {}).get("examples", []):
        print(example["name"], "RK4 errors:", example["mode_errors"], "refinement:", example["refinement_factors"])


if __name__ == "__main__":
    main()
