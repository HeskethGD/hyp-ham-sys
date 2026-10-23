"""Branch-point Onishi identities, composite factors and lift balancing.

This is separate from verification.py: it tests the full sigma determinant
formula and its residual factorization, not just the algebraic determinant.
"""

import argparse
import json
from pathlib import Path
import sys
from time import perf_counter

import sympy as sp

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "code"))
from hyperelliptic_scalar_first.branch_points.model import BranchSigmaModes, example
from hyperelliptic_scalar_first.abel_lifts import integer_lattice_coordinates


def verify_determinants(genus, dps=35, progress=print):
    started = perf_counter()
    model = BranchSigmaModes(example(genus), dps=dps, progress=progress)
    c, g, n = model.ctx, genus, model.n
    rows = []
    tolerance = c.mpf(10)**(-(dps-12))
    def record(name, actual, expected):
        error = abs(actual-expected)/max(c.one, abs(actual), abs(expected))
        rows.append(dict(check=name, scaled_error=c.nstr(error, 12), passed=bool(error < tolerance)))

    total = sum(model.eta, c.matrix(g, 1))
    integers, residual = integer_lattice_coordinates(model, total)
    for k in range(g):
        record(f"unbalanced total Abel sum is a period, coordinate {k+1}", residual[k], 0)
    second = model.curve.periods_kind_2()
    m, nn = c.matrix(integers[:g]), c.matrix(integers[g:])
    # Shift one unit-weight point by minus the total first-kind period.
    # h=-integral dr shifts by minus the matching 2*eta*m+2*eta_prime*n.
    h_shift = 2*second.eta*m+2*second.eta_prime*nn
    eta = [v.copy() for v in model.eta]
    h = [v[:] for v in model.h]
    period_shift = 2*model.periods.omega*m+2*model.periods.omega_prime*nn
    eta[-1] -= period_shift
    h[-1] = [h[-1][k]-h_shift[k] for k in range(g)]
    balanced_sum = sum(eta, c.matrix(g, 1))
    normals = [model.jet(v)[model.sharp] for v in eta]
    for k in range(g):
        record(f"balanced total Abel sum {k+1}", balanced_sum[k], 0)
        record(f"balanced total h {k+1}", sum(v[k] for v in h), 0)
    # Independently extract h at the shifted lift from its sigma derivatives.
    jet = model.jet(eta[-1])
    zeta = [jet[model.extend(model.sharp, k)]/normals[-1] for k in range(g)]
    h2 = zeta[1]/2
    shifted_h = [zeta[0]-jet[(1,)+(0,)*(g-1)]*h2/normals[-1], h2]
    if g == 3:
        shifted_h.append(zeta[2])
    for k in range(g):
        record(f"period-shifted h sigma versus second-kind periods {k+1}", shifted_h[k], h[-1][k])

    def sigma(w, index=None):
        return model.jet(w)[index or (0,)*g]

    def fs(vectors):
        """Full corrected Theorem 7.2(2); n>=g in every call here."""
        size = len(vectors)
        exponent = (size-g-1)*(size+g*g+2*g)//2
        sign = (-1)**exponent
        top = sigma(sum(vectors, c.matrix(g, 1)))
        for i in range(size):
            for j in range(i+1, size):
                top *= sigma(vectors[i]-vectors[j], model.flat)
        bottom = c.fprod(sigma(v, model.sharp)**size for v in vectors)
        return sign*top/bottom

    def raw(w, lifts, which):
        denominator = sigma(w, model.sharp)
        return [sigma(w+which*v, model.flat)/(denominator*sigma(v, model.sharp)) for v in lifts]

    A, B = list(range(g+1)), list(range(g+1, n))
    E = sum((eta[j] for j in A), c.matrix(g, 1))
    delta_A = c.fprod(model.a[j]-model.a[i] for p, i in enumerate(A) for j in A[p+1:])
    delta_B = c.fprod(model.a[j]-model.a[i] for p, i in enumerate(B) for j in B[p+1:])
    samples = []
    for sheet in (1, -1):
        progress(f"genus {g}: determinant sigma identities, sheet {sheet}")
        start = model.evaluate(model.rho_start, sheet)
        raw0u, raw0v = raw(start["w"], eta, 1), raw(start["w"], eta, -1)
        unbalanced0u = raw(start["w"], model.eta, 1)
        initial = model.initial(start)
        for fraction in (c.zero, c.mpf("0.5"), c.one):
            x = model.rho_start+fraction*(model.rho_end-model.rho_start)
            w = start["w"]+model.abel_increment(model.rho_start, x, sheet)
            Y = sheet*c.sqrt(model.F(x))
            tag = f"sheet {sheet}, fraction {fraction}"
            # Actual determinant, rather than its predicted closed form.
            matrix = c.matrix([[a**k for k in range(g+1)]+[c.zero] for a in model.a[:g+1]]
                              +[[x**k for k in range(g+1)]+[Y/2]])
            D = c.det(matrix)
            record(tag+": algebraic determinant", D, delta_A*Y/2)
            record(tag+": full main sigma determinant", fs([eta[j] for j in A]+[w]), D)
            record(tag+": full opposite-lift sigma determinant", fs([-eta[j] for j in A]+[w]), D)
            auxiliary = c.det(c.matrix([[model.a[j]**k for k in range(g+1)] for j in B]
                                       +[[x**k for k in range(g+1)]]))
            record(tag+": residual Vandermonde", auxiliary, delta_B*c.fprod(x-model.a[j] for j in B))
            record(tag+": full auxiliary sigma determinant", fs([-eta[j] for j in B]+[w]), auxiliary)
            ru, rv = raw(w, eta, 1), raw(w, eta, -1)
            unbalanced_ru = raw(w, model.eta, 1)
            record(tag+": balanced raw u product ratio", c.fprod(ru)/c.fprod(raw0u), Y/start["Y"])
            record(tag+": balanced raw v product ratio", c.fprod(rv)/c.fprod(raw0v), Y/start["Y"])
            composite = sigma(w+E)/sigma(w, model.sharp)
            samples.append(dict(tag=tag, composite_over_residual=composite/c.fprod(rv[j] for j in B),
                                main_over_raw_v=D/c.fprod(rv),
                                unbalanced_raw_product_ratio=c.fprod(unbalanced_ru)/c.fprod(unbalanced0u),
                                algebraic_product_ratio=Y/start["Y"]))
            dw = w-start["w"]
            u = [initial[j]*ru[j]/raw0u[j]*c.exp(-sum(h[j][k]*dw[k] for k in range(g))) for j in range(n)]
            v = [initial[n+j]*rv[j]/raw0v[j]*c.exp(sum(h[j][k]*dw[k] for k in range(g))) for j in range(n)]
            record(tag+": corrected U from sigma equals determinant", c.fprod(u), D/delta_A)
            record(tag+": corrected V from sigma equals determinant", c.fprod(v), -D/delta_A)
            for j in range(n):
                old_factor = c.exp(-sum(model.h[j][k]*dw[k] for k in range(g)))
                record(tag+f": mode {j+1} invariant under compatible lift shift",
                       u[j], initial[j]*unbalanced_ru[j]/unbalanced0u[j]*old_factor)
    for sample in samples[1:]:
        record(sample["tag"]+": composite factor constant", sample["composite_over_residual"], samples[0]["composite_over_residual"])
        record(sample["tag"]+": complete raw factorization constant", sample["main_over_raw_v"], samples[0]["main_over_raw_v"])
    unbalanced_defect = max(abs(s["unbalanced_raw_product_ratio"]-s["algebraic_product_ratio"]) for s in samples)
    # A genuinely mixed fixed-row matrix on the same curve. Interpolate q
    # from g branch rows and one regular row, rather than testing only q=0.
    fixed_x = model.rho_start+c.mpf("0.1")
    fixed_y = c.sqrt(model.F(fixed_x))
    fixed_w = model.anchor_w+model.abel_increment(model.rho_start, fixed_x, 1)
    mixed_x = model.a[:g]+[fixed_x]
    mixed_y = [c.zero]*g+[fixed_y]
    mixed_lifts = model.eta[:g]+[fixed_w]
    vandermonde = c.matrix([[a**k for k in range(g+1)] for a in mixed_x])
    coefficients = c.lu_solve(vandermonde, c.matrix(mixed_y))
    mixed_delta = c.det(vandermonde)
    for sheet in (1,-1):
        start = model.evaluate(model.rho_start, sheet)
        for fraction in (c.mpf("0.5"), c.one):
            x = model.rho_start+fraction*(model.rho_end-model.rho_start)
            w = start["w"]+model.abel_increment(model.rho_start,x,sheet)
            Y = sheet*c.sqrt(model.F(x))
            q = sum(coefficients[k]*x**k for k in range(g+1))
            D = c.det(c.matrix([[a**k for k in range(g+1)]+[b/2] for a,b in zip(mixed_x,mixed_y)]
                              +[[x**k for k in range(g+1)]+[Y/2]]))
            Dstar = c.det(c.matrix([[a**k for k in range(g+1)]+[-b/2] for a,b in zip(mixed_x,mixed_y)]
                                  +[[x**k for k in range(g+1)]+[Y/2]]))
            tag = f"mixed rows, sheet {sheet}, fraction {fraction}"
            record(tag+": interpolation", D, mixed_delta*(Y-q)/2)
            record(tag+": opposite interpolation", Dstar, mixed_delta*(Y+q)/2)
            record(tag+": full sigma determinant", fs(mixed_lifts+[w]), D)
            record(tag+": full opposite sigma determinant", fs([-v for v in mixed_lifts]+[w]), Dstar)
    return dict(genus=g, dps=dps, integer_period_coordinates=integers,
                unbalanced_sum_h=[c.nstr(sum(v[k] for v in model.h), 18) for k in range(g)],
                unbalanced_raw_product_defect=c.nstr(unbalanced_defect, 12),
                check_count=len(rows), all_passed=all(row["passed"] for row in rows),
                maximum_scaled_error=c.nstr(max(c.mpf(r["scaled_error"]) for r in rows), 12),
                tolerance=c.nstr(tolerance, 8), checks=rows, seconds=round(perf_counter()-started, 3))


def exact_determinants():
    """Check both column families with a mixture of branch/nonbranch rows.

    For each genus choose q of the largest permitted degree, with q(a0)=0.
    Define F=q^2+4 product(x-a_j), which is odd degree and squarefree for
    these test data. The branch row has zero ordinate; all rows stay distinct.
    """
    x, y = sp.symbols("x y")
    checks = {}
    for g in range(2, 7):
        for extended in (False, True):
            degree = g+1 if extended else g
            points = list(range(-g, g+1))
            q = sp.expand((x-points[0])*(x**(degree-1)+2))
            F = sp.expand(q*q+4*sp.prod(x-a for a in points))
            # For degree g+1 use 2g+2 intersections and leading q=2*x^(g+1)
            # so the even leading term cancels while fixed points stay on C.
            if extended:
                points = list(range(-g-1, g+1))
                q = sp.expand(2*(x-points[0])*(x**g+2))
                F = sp.expand(q*q-4*sp.prod(x-a for a in points))
            tag = f"g={g}, n=g+{3 if extended else 2}"
            checks[tag+": odd smooth curve"] = bool(sp.degree(F, x)==2*g+1 and sp.degree(sp.gcd(F, sp.diff(F,x)),x)==0)
            size = g+3 if extended else g+2
            chosen = points[:size-1]
            def row(a, b):
                base = [a**k for k in range(g+1)]+[b/2]
                return base+[a**(g+1)] if extended else base
            D = sp.det(sp.Matrix([row(sp.Integer(a), q.subs(x,a)) for a in chosen]+[row(x,y)]))
            Dstar = sp.det(sp.Matrix([row(sp.Integer(a), -q.subs(x,a)) for a in chosen]+[row(x,y)]))
            delta = sp.prod(chosen[j]-chosen[i] for i in range(len(chosen)) for j in range(i+1,len(chosen)))
            epsilon = -1 if extended else 1
            checks[tag+": branch row"] = q.subs(x,chosen[0])==0 and sp.diff(F,x).subs(x,chosen[0])!=0
            checks[tag+": main interpolation"] = sp.expand(D-epsilon*delta*(y-q)/2)==0
            checks[tag+": opposite-sheet interpolation"] = sp.expand(Dstar-epsilon*delta*(y+q)/2)==0
    return checks


def run_verification(dps=35, progress=print):
    exact = exact_determinants()
    reports = [verify_determinants(g, dps, progress) for g in (2,3)]
    return dict(scope="Arbitrary-genus proofs are in proof.md. Numerical full sigma determinant and residual checks cover genera 2 and 3; exact mixed-row tests cover genera 2 through 6.",
                exact_checks=exact, examples=reports,
                check_count=len(exact)+sum(r["check_count"] for r in reports),
                all_passed=all(exact.values()) and all(r["all_passed"] for r in reports))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=35)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("determinant_results.json"))
    args = parser.parse_args()
    report = run_verification(args.dps, lambda s: print(s, flush=True))
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(f"{report['check_count']} checks; all passed={report['all_passed']}")
    for k,v in report["exact_checks"].items():
        if not v: print(k, v)
    for r in report["examples"]:
        print(f"genus {r['genus']}: max error {r['maximum_scaled_error']}")
        for row in r["checks"]:
            if not row["passed"]: print(row)
    if not report["all_passed"]: raise SystemExit(1)
