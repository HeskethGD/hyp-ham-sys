"""Independent exact, sigma, quadrature and full-field trajectory checks.

Run: .venv/bin/python -m hyperelliptic_scalar_first.verification
No helpers from earlier research notebooks are imported.
"""
import argparse
import json
from pathlib import Path
import sympy as sp
from .model import example, genus_three_example, scalar_data, canonical_system
from .sigma_jets import extract_h, explicit_h_at_zero, fixed_point_h


def symbolic_checks():
    checks = {}
    rho, y = sp.symbols('rho y')
    for beta, gamma in [((1, 1, 1), (-1, 0, 1)),
                        ((2, 1, 1, 1), (-4, -1, 1, 4)),
                        ((3, 2, 1), (-3, 1, 2))]:
        x = sp.symbols(f'x1:{len(beta)+1}')
        mean = -sum(x)/sum(beta)
        Q = 1+mean+mean**2+x[0]*x[1]/7
        data = scalar_data(x, Q, beta, gamma, rho)
        q, r, F = (data[k] for k in ('q', 'r', 'F'))
        L = [b*(q-y)/(2*power)-partial for b, power, partial in
             zip(beta, data['powers'], data['partials'])]
        M = [-b*(q+y)/(2*power)+partial for b, power, partial in
             zip(beta, data['powers'], data['partials'])]
        tangent = lambda e: y*sp.diff(e, rho)+sp.diff(F, rho)*sp.diff(e, y)/2
        def zero(e):
            num = sp.together(e).as_numer_denom()[0]
            return sp.Poly(num, y).rem(sp.Poly(y*y-F, y)).is_zero
        tag = str(beta)
        checks[tag+' u product propagation'] = zero(tangent((q+y)/2)-(q+y)*sum(b*l for b, l in zip(beta, L))/2)
        checks[tag+' v product propagation'] = zero(tangent((q-y)/2)-(q-y)*sum(b*m for b, m in zip(beta, M))/2)
        checks[tag+' product norm'] = zero((q+y)*(q-y)/4-r)
        checks[tag+' scalar acceleration chain rule'] = sp.expand(
            sp.diff(F, rho)/2+q*sum(b*a for b, a in zip(beta, data['partials']))
            -2*sum(b*b*r/p for b, p in zip(beta, data['powers']))) == 0
        checks[tag+' point-factor split'] = all(sp.cancel(
            l-(fixed-y)/(2*(a-rho))-B) == 0
            for l, fixed, a, B in zip(L, data['b'], data['a'], data['B']))
        field = canonical_system(x, Q, beta)
        u, v, H, U, V = (field[k] for k in ('u', 'v', 'H', 'U', 'V'))
        bracket = lambda A, B: sum(sp.diff(A, a)*sp.diff(B, b)-sp.diff(A, b)*sp.diff(B, a) for a, b in zip(u, v))
        powers = [a*b for a, b in zip(u, v)]
        C = [beta[-1]*powers[j]-beta[j]*powers[-1] for j in range(len(beta)-1)]
        checks[tag+' weighted power rates'] = all(sp.expand(du*b+a*dv-be*(V-U)) == 0
            for a, b, be, du, dv in zip(u, v, beta, field['du'], field['dv']))
        checks[tag+' conserved powers commute with H'] = all(sp.expand(bracket(c, H)) == 0 for c in C)
        checks[tag+' mutual involution'] = all(sp.expand(bracket(a, b)) == 0 for a in C for b in C)
        witness = dict(zip(u+v, range(1, 2*len(beta)+1)))
        checks[tag+' independence witness'] = sp.Matrix([[sp.diff(c, a) for a in u+v] for c in C+[H]]).subs(witness).rank() == len(beta)
    data = example()
    checks['example smooth genus two'] = data['squarefree'] and data['genus'] == 2
    checks['example nonbranch fixed points'] = all(b != 0 for b in data['b'])
    checks['example affine gauges'] = data['holomorphic_gauges']
    checks['example leading coefficient four'] = sp.LC(sp.Poly(data['F'], data['rho'])) == 4
    checks['example weighted centre'] = sum(data['gamma']) == 0
    checks['example repeated intersection at double mode'] = sp.diff(data['F']-data['q']**2, data['rho']).subs(data['rho'], data['a'][0]) == 0
    for b, a, fixed in zip(data['beta'], data['a'], data['b']):
        checks[f'example contact order {b} at {a}'] = all(sp.diff(data['F']-data['q']**2, data['rho'], k).subs(data['rho'], a) == 0 for k in range(b))
    # A low scalar degree alone does not bound individual phases.
    x1, x2, x3 = sp.symbols('x1 x2 x3')
    high = scalar_data((x1, x2, x3), 1-(x1+x2+x3)/3+(x1-x2+1)*x3**4,
                       (1, 1, 1), (-1, 0, 1), rho)
    checks['reject high-degree off-line phase despite low q'] = not high['holomorphic_gauges']
    g3 = genus_three_example()
    checks['genus three smooth degree-seven curve'] = g3['squarefree'] and sp.degree(g3['F'], g3['rho']) == 7
    checks['genus three leading coefficient four'] = sp.LC(sp.Poly(g3['F'], g3['rho'])) == 4
    checks['genus three nonbranch fixed points'] = all(b != 0 for b in g3['b'])
    checks['genus three quadratic gauges'] = g3['holomorphic_gauges'] and any(sp.degree(B, g3['rho']) == 2 for B in g3['B'])
    # Regression for the actual form-feed corruption in an example's LaTeX.
    import ast
    source = Path(__file__).with_name('notebook.py').read_text()
    markdown = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == 'md' and node.args:
            argument = node.args[0]
            if isinstance(argument, ast.Constant) and isinstance(argument.value, str):
                markdown.append(argument.value)
            elif isinstance(argument, ast.JoinedStr):
                markdown.append(''.join(value.value for value in argument.values if isinstance(value, ast.Constant) and isinstance(value.value, str)))
    checks['notebook markdown has no escaped control characters'] = all(
        ord(character) >= 32 or character in '\n\t' for block in markdown for character in block)
    checks['notebook display math delimiters balanced'] = all(block.count(r'\[') == block.count(r'\]') for block in markdown)
    checks['notebook example fraction uses braced LaTeX'] = any(
        r'\frac{1}{3}' in block for block in markdown)
    checks = {k: bool(v) for k, v in checks.items()}
    return dict(passed=all(checks.values()), checks=checks)


class SigmaExample:
    """Explicit genus-two/three sigma point factors on normalized curves."""
    def __init__(self, dps=30, genus=2):
        from mpmath import mp
        from generapy import Curve, kleinian_sigma_jet
        self.ctx = c = mp.clone(); c.dps = dps
        if genus not in (2, 3):
            raise ValueError('Explicit sigma implementation supports genus 2 or 3')
        self.data = d = example() if genus == 2 else genus_three_example()
        self.genus = genus; self.N = len(d['beta'])
        self.beta = d['beta']; self.gamma = [c.mpf(str(a)) for a in d['gamma']]
        self.a = [c.mpf(str(a)) for a in d['a']]
        assert sp.LC(sp.Poly(d['F'], d['rho'])) == 4
        self.sharp = tuple(int(k == 1) for k in range(genus))
        self.flat = (0, 0) if genus == 2 else (0, 0, 1)
        self.rho0 = c.mpf('4.5') if genus == 2 else c.mpf('5.5')
        self.sample_rho = [self.rho0+c.mpf(value) for value in ('.05', '.1', '.15')]
        # Polynomial arithmetic uses this instance's precision context.
        self.coefficients = [c.mpf(str(d['F'].coeff(d['rho'], k))) for k in range(2*genus+2)]
        self.F = lambda rho: sum(a*rho**k for k, a in enumerate(self.coefficients))
        qcoeff = [c.mpf(str(d['q'].coeff(d['rho'], k))) for k in range(sp.degree(d['q'], d['rho'])+1)]
        self.q = lambda rho: sum(a*rho**k for k, a in enumerate(qcoeff))
        self.qprime = lambda rho: sum(k*a*rho**(k-1) for k, a in enumerate(qcoeff) if k)
        self.b = [self.q(a) for a in self.a]
        self.B = [[c.mpf(str(sp.expand(B).coeff(d['rho'], k))) for k in range(genus)] for B in d['B']]
        self.curve = Curve({(0, 2): 1}|{(k, 0): -a for k, a in enumerate(self.coefficients) if a}, ctx=c)
        first = self.curve.periods_kind_1(); second = self.curve.periods_kind_2()
        options = dict(omega=first.omega, tau=first.tau, kappa=second.kappa,
                       characteristic=self.curve.riemann_constant().characteristic,
                       normalization='hyperelliptic', ctx=c)
        self.jet = lambda w, order=2: kleinian_sigma_jet(w, order=order, **options)
        self.eta = [self.curve.abel_map_kind_1((a, b), reduce=False).value for a, b in zip(self.a, self.b)]
        self.normals = [self.jet(eta)[self.sharp] for eta in self.eta]
        self.h = fixed_point_h(self)
        self.nodes, self.weights = c.gauss_quadrature(32, 'legendre')

    def evaluate(self, rho, sheet=1):
        c = self.ctx; y = sheet*c.sqrt(self.F(rho))
        w = self.curve.abel_map_kind_1((rho, y), reduce=False).value
        jet = self.jet(w)
        def rate(value, index):
            derivatives = []
            for k in range(self.genus):
                extended = tuple(a+int(i == k) for i, a in enumerate(index))
                derivatives.append(rho**k*value[extended])
            return sum(derivatives)/value[index]
        normalrate = rate(jet, self.sharp)
        raw_u = []; raw_v = []; rates = []; raw_v_rates = []; raw_power = []
        for be, eta, norm in zip(self.beta, self.eta, self.normals):
            plus = self.jet(w+eta); minus = self.jet(w-eta)
            denominator = jet[self.sharp]*norm
            raw_u.append(c.sqrt(be)*plus[self.flat]/denominator)
            raw_v.append((-1)**(self.genus+1)*c.sqrt(be)*minus[self.flat]/denominator)
            raw_power.append(raw_u[-1]*raw_v[-1])
            rates.append(rate(plus, self.flat)-normalrate)
            raw_v_rates.append(rate(minus, self.flat)-normalrate)
        inversion = -jet[(1, 0)]/jet[(0, 1)] if self.genus == 2 else -jet[(1, 0, 1)]/jet[(0, 1, 1)]
        return dict(rho=rho, y=y, w=w, raw_u=raw_u, raw_v=raw_v, rates=rates,
                    raw_v_rates=raw_v_rates, raw_power=raw_power,
                    inversion=inversion, sigma=jet[(0,)*self.genus],
                    sigma3=jet[(0, 0, 1)] if self.genus == 3 else c.zero)

    def bare(self, item, j):
        return (self.b[j]-item['y'])/(2*(self.a[j]-item['rho']))

    def integral(self, function, start, end):
        half = (end-start)/2; mid = (start+end)/2
        return half*sum(weight*function(mid+half*node) for node, weight in zip(self.nodes, self.weights))

    def point_at_z(self, z_minus_z0, rho0, sheet=1):
        """Invert the scalar quadrature on this short, nonbranch local chart.

        The predictor fixes the local branch. Paths must remain in the same
        regular chart; this is not a turning-point continuation algorithm.
        """
        c = self.ctx
        y = lambda rho: sheet*c.sqrt(self.F(rho))
        if z_minus_z0 == 0:
            return self.evaluate(rho0, sheet)
        predictor = rho0+y(rho0)*z_minus_z0
        rho = c.findroot(lambda coordinate: self.integral(lambda value: 1/y(value), rho0, coordinate)-z_minus_z0,
                         predictor, solver='newton', df=lambda coordinate: 1/y(coordinate))
        return self.evaluate(rho, sheet)

    def partials(self, powers):
        rho = -sum(powers)/sum(self.beta)
        common = -self.qprime(rho)/sum(self.beta)
        constraint = powers[0]/2-powers[1]+self.a[1]-self.a[0]
        if self.genus == 2:
            return [common+powers[2]/6, common-powers[2]/3, common+constraint/3]+[common]*(self.N-3)
        return [common+powers[2]**2/14, common-powers[2]**2/7,
                common+2*constraint*powers[2]/7]+[common]*(self.N-3)

    def Q(self, powers):
        rho = -sum(powers)/sum(self.beta)
        correction = powers[2]/3 if self.genus == 2 else powers[2]**2/7
        return self.q(rho)+(powers[0]/2-powers[1]+self.a[1]-self.a[0])*correction

    def field(self, state):
        c = self.ctx; u, v = state[:self.N], state[self.N:]
        partials = self.partials([a*b for a, b in zip(u, v)])
        # Polynomial field: never divide by a mode, including zero modes.
        du = [be*c.fprod(v[k]**(self.beta[k]-(k == j)) for k in range(self.N))-partials[j]*u[j]
              for j, be in enumerate(self.beta)]
        dv = [-be*c.fprod(u[k]**(self.beta[k]-(k == j)) for k in range(self.N))+partials[j]*v[j]
              for j, be in enumerate(self.beta)]
        return du+dv

    def initial(self, start):
        c = self.ctx
        u = [c.mpc('1.1', '.2'), c.mpc('.8', '-.1'), c.mpc('1.2', '.1')]+[c.mpc(1+c.mpf(j)/10, '.1') for j in range(self.N-3)]
        u[-1] = (self.q(start['rho'])+start['y'])/(2*c.fprod(u[j]**self.beta[j] for j in range(self.N-1)))
        v = [(gamma-be*start['rho'])/value for gamma, be, value in zip(self.gamma, self.beta, u)]
        return u+v

    def solution(self, start, end, initial):
        c = self.ctx
        dw = (end['w']-start['w'])
        point_factors = [last/first*c.exp(-sum(hk*increment for hk, increment in zip(h, dw)))
                         for first, last, h in zip(start['raw_u'], end['raw_u'], self.h)]
        polynomial_integrals = [sum(ck*increment for ck, increment in zip(B, dw)) for B in self.B]
        u = [value*factor*c.exp(integral) for value, factor, integral in
             zip(initial[:self.N], point_factors, polynomial_integrals)]
        v = [(gamma-be*end['rho'])/value for gamma, be, value in zip(self.gamma, self.beta, u)]
        return u+v


def rk4(field, initial, length, steps):
    state = list(initial); h = length/steps
    for _ in range(steps):
        k1 = field(state)
        k2 = field([a+h*b/2 for a, b in zip(state, k1)])
        k3 = field([a+h*b/2 for a, b in zip(state, k2)])
        k4 = field([a+h*b for a, b in zip(state, k3)])
        state = [a+h*(b+2*c+2*d+e)/6 for a, b, c, d, e in zip(state, k1, k2, k3, k4)]
    return state


def numerical_checks(dps=30, genus=2):
    model = SigmaExample(dps, genus); c = model.ctx; N = model.N
    rows = []; tol = c.mpf(10)**(-(dps-8))
    def record(name, actual, expected, tolerance=tol):
        error = abs(actual-expected)/max(c.one, abs(actual), abs(expected))
        rows.append(dict(check=name, scaled_error=c.nstr(error, 8), passed=bool(error < tolerance)))
    trajectories = []
    # Three independent realizations of the same fixed-point coefficients:
    # short sigma primitives, the old moving-reference formula, and direct
    # regularized second-kind Abel integrals with the same unreduced lifts.
    reference_h = explicit_h_at_zero(model)
    for j, (actual, expected) in enumerate(zip(model.h, reference_h)):
        second_kind = model.curve.abel_map_kind_2((model.a[j], model.b[j]), reduce=False).value
        for k, (a, b) in enumerate(zip(actual, expected)):
            record(f'h_{j+1},{k+1}: fixed-point sigma versus reference formula', a, b)
            record(f'h_{j+1},{k+1}: fixed-point sigma versus second-kind Abel integral', a, -second_kind[k])
    # Independent finite Taylor composition checks the expanded zeta formulas.
    for reference in (c.zero, model.rho0):
        derived = extract_h(model, reference=reference)
        for j, (actual, expected) in enumerate(zip(derived, model.h)):
            for k, (a, b) in enumerate(zip(actual, expected)):
                record(f'h_{j+1},{k+1}: finite sigma jets at rho={reference}', a, b)
    for sheet in (1, -1):
        start = model.evaluate(model.rho0, sheet); initial = model.initial(start)
        for coordinate in model.sample_rho:
            end = model.evaluate(coordinate, sheet)
            tag = f'sheet {sheet}, rho={coordinate}'
            record(tag+' inversion', end['inversion'], end['rho'])
            record(tag+' sigma divisor', end['sigma'], 0)
            if genus == 3:
                record(tag+' sigma_3 one-point constraint', end['sigma3'], 0)
            dw = (end['w']-start['w'])
            for k in range(genus):
                direct = model.integral(lambda rho: rho**k/(sheet*c.sqrt(model.F(rho))), start['rho'], end['rho'])
                record(tag+f' Abel increment {k+1}', dw[k], direct)
            state = model.solution(start, end, initial)
            field = model.field(state); u, v = state[:N], state[N:]
            powers = [a*b for a, b in zip(u, v)]; partials = model.partials(powers)
            U = c.fprod(a**b for a, b in zip(u, model.beta)); V = c.fprod(a**b for a, b in zip(v, model.beta))
            record(tag+' weighted U', U, (model.q(end['rho'])+end['y'])/2)
            record(tag+' weighted V', V, (model.q(end['rho'])-end['y'])/2)
            record(tag+' Hamiltonian', U+V-model.Q(powers), 0)
            for j in range(N):
                record(tag+f' raw sigma mode power {j+1}', end['raw_power'][j], model.gamma[j]-model.beta[j]*end['rho'])
                record(tag+f' paired raw sigma mode rates {j+1}', end['rates'][j]+end['raw_v_rates'][j], -model.beta[j]*end['y']/powers[j])
                correction = sum(h*end['rho']**k for k, h in enumerate(model.h[j]))
                record(tag+f' holomorphic discrepancy {j+1}', end['rates'][j]-model.bare(end, j), correction)
                B = sum(a*end['rho']**k for k, a in enumerate(model.B[j]))
                rate = end['rates'][j]+B-correction
                target = model.beta[j]*(model.q(end['rho'])-end['y'])/(2*powers[j])-partials[j]
                record(tag+f' reconstructed log rate {j+1}', rate, target)
                record(tag+f' canonical u {j+1}', u[j]*rate, field[j])
                record(tag+f' canonical v {j+1}', v[j]*(-model.beta[j]*end['y']/powers[j]-rate), field[j+N])
                direct = model.integral(lambda rho: sum(a*rho**k for k, a in enumerate(model.B[j]))/(sheet*c.sqrt(model.F(rho))), start['rho'], end['rho'])
                record(tag+f' polynomial gauge integral {j+1}', direct, sum(a*increment for a, increment in zip(model.B[j], dw)))
                bare_integral = model.integral(lambda rho: (model.b[j]-sheet*c.sqrt(model.F(rho)))/(2*(model.a[j]-rho)*sheet*c.sqrt(model.F(rho))), start['rho'], end['rho'])
                h_integral = sum(h*increment for h, increment in zip(model.h[j], dw))
                primitive = c.log(end['raw_u'][j]/start['raw_u'][j])-h_integral
                record(tag+f' third-kind primitive {j+1}', primitive, bare_integral)
                third_integral = model.integral(lambda rho: model.b[j]/((model.a[j]-rho)*sheet*c.sqrt(model.F(rho))), start['rho'], end['rho'])
                sigma_ratio = (end['raw_u'][j]/end['raw_v'][j])/(start['raw_u'][j]/start['raw_v'][j])
                record(tag+f' paired third-kind sigma identity {j+1}', c.log(sigma_ratio)/2-h_integral, third_integral/2)
                record(tag+f' split integral equals point primitive {j+1}', c.log((model.a[j]-end['rho'])/(model.a[j]-start['rho']))/2+third_integral/2, primitive)
        end = model.evaluate(model.sample_rho[-1], sheet)
        exact = model.solution(start, end, initial)
        length = model.integral(lambda rho: 1/(sheet*c.sqrt(model.F(rho))), start['rho'], end['rho'])
        inverted = model.point_at_z(length, start['rho'], sheet)
        record(f'sheet {sheet} scalar quadrature inversion', inverted['rho'], end['rho'])
        inverted_state = model.solution(start, inverted, initial)
        for j, (actual, expected) in enumerate(zip(inverted_state, exact)):
            record(f'sheet {sheet} mode {j+1} from prescribed z', actual, expected)
        errors = []
        for steps in (16, 32, 64):
            result = rk4(model.field, initial, length, steps)
            error = max(abs(a-b)/max(c.one, abs(b)) for a, b in zip(result, exact))
            errors.append(error)
            powers = [a*b for a, b in zip(result[:N], result[N:])]
            U = c.fprod(a**b for a, b in zip(result[:N], model.beta)); V = c.fprod(a**b for a, b in zip(result[N:], model.beta))
            record(f'sheet {sheet} RK4 {steps} Hamiltonian', U+V-model.Q(powers), 0, c.mpf('1e-7'))
            for j in range(N-1):
                record(f'sheet {sheet} RK4 {steps} invariant {j+1}', powers[j]-model.beta[j]*powers[-1], model.gamma[j]-model.beta[j]*model.gamma[-1], c.mpf('1e-7'))
        ratios = [errors[k]/errors[k+1] for k in range(2)]
        record(f'sheet {sheet} RK4 finest agreement', errors[-1], 0, c.mpf('1e-8'))
        for k, ratio in enumerate(ratios):
            rows.append(dict(check=f'sheet {sheet} RK4 fourth-order ratio {k+1}',
                             scaled_error=c.nstr(abs(ratio-16)/16, 8), passed=bool(12 < ratio < 20)))
        trajectories.append(dict(sheet=sheet, z_length=c.nstr(length, 12), steps=[16, 32, 64],
                                 errors=[c.nstr(e, 8) for e in errors], ratios=[c.nstr(e, 8) for e in ratios]))
    # Full polynomial field remains finite at zero modes, outside our log chart.
    zero_state = list(initial); zero_state[0] = c.zero
    rows.append(dict(check='polynomial field finite at a zero mode', scaled_error='0', passed=all(c.isfinite(a) for a in model.field(zero_state))))
    return dict(passed=all(row['passed'] for row in rows), dps=dps, genus=genus, rows=rows, trajectories=trajectories,
                scope=f'Fresh genus-{genus} degenerate {N}-pair example, two sheets, regular short paths; no turning-point or global continuation claim.')


def run_all(dps=30):
    report = dict(symbolic=symbolic_checks(), numerical=numerical_checks(dps), genus_three=numerical_checks(dps, genus=3))
    report['passed'] = report['symbolic']['passed'] and report['numerical']['passed'] and report['genus_three']['passed']
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--dps', type=int, default=30)
    args = parser.parse_args(); report = run_all(args.dps)
    target = Path(__file__).with_name('verification_results.json')
    target.write_text(json.dumps(report, indent=2)+'\n')
    print('Exact checks:', len(report['symbolic']['checks']), report['symbolic']['passed'])
    for key in ('numerical', 'genus_three'):
        print(f"Genus {report[key]['genus']} numerical checks:", len(report[key]['rows']), report[key]['passed'])
        for row in report[key]['rows']:
            if not row['passed']: print('FAILED:', row)
        print('Trajectories:', report[key]['trajectories'])
    print('Report:', target)
    raise SystemExit(0 if report['passed'] else 1)
