"""Leading-four curve normalization and corrected sigma point factors.

The canonical fields remain in the original Hamiltonian scaling. Only the
curve ordinate and its Abel differential are normalized for Generapy.
"""

import sympy as sp
from sympy.printing.pycode import MpmathPrinter
from mpmath import mp
from generapy import Curve, kleinian_sigma_jet


def normalized_curve(data):
    rho = data["rho"]
    leading = sp.Poly(data["F"], rho).LC()
    scale = sp.sqrt(leading/4)
    normalized = sp.expand(data["F"]/(scale*scale))
    q_normalized = data["q"]/scale
    checks = {
        "positive leading coefficient": bool(leading > 0),
        "leading four": sp.Poly(normalized, rho).LC() == 4,
        "ordinate scaling": sp.expand(scale**2*normalized-data["F"]) == 0,
        "normalized curve norm": sp.expand(q_normalized**2-normalized-4*data["r"]/(scale*scale)) == 0,
        "special corrections": all(value == 0 for value in data["B"]),
        "smooth odd degree": data["squarefree"] and data["degree"] % 2 == 1,
    }
    for j, (point, ordinate) in enumerate(zip(data["a"], data["ordinate"])):
        quotient = sp.cancel((data["F"]-ordinate**2)/(rho-point))
        normalized_quotient = sp.cancel((normalized-(ordinate/scale)**2)/(rho-point))
        checks[f"mode {j+1} physical potential scaling"] = sp.cancel(
            sp.diff(quotient, rho)-scale**2*sp.diff(normalized_quotient, rho)) == 0
    if not all(checks.values()):
        raise ValueError(f"Unsupported normalization: {checks}")
    return dict(scale=scale, F=normalized, checks=checks)


class SigmaModes:
    def __init__(self, data, dps=35, progress=None):
        self.data = data
        self.normalization = normalized_curve(data)
        self.ctx = c = mp.clone()
        c.dps = dps
        self.genus = g = data["genus"]
        if g not in (2, 3):
            raise ValueError("The explicit fixed-point sigma primitives cover genera two and three")
        self.n = len(data["beta"])
        self.sharp = tuple(int(k == 1) for k in range(g))
        self.flat = (0, 0) if g == 2 else (0, 0, 1)
        self.scale = self.number(self.normalization["scale"])
        self.a = [self.number(x) for x in data["a"]]
        self.b_original = [self.number(x) for x in data["ordinate"]]
        self.b = [x/self.scale for x in self.b_original]
        self.rho_start = self.number(data["rho_start"])
        self.rho_end = self.number(data["rho_end"])
        rho = data["rho"]
        self.F_original = self.compile(rho, data["F"])
        self.F = self.compile(rho, self.normalization["F"])
        self.q = self.compile(rho, data["q"])
        coefficients = [self.number(self.normalization["F"].coeff(rho, k))
                        for k in range(data["degree"]+1)]
        self.curve = Curve({(0, 2): 1} | {(k, 0): -a for k, a in enumerate(coefficients) if a}, ctx=c)
        if progress:
            progress("computing periods")
        first = self.curve.periods_kind_1()
        second = self.curve.periods_kind_2()
        self.periods = first
        self.sigma_options = dict(omega=first.omega, tau=first.tau, kappa=second.kappa,
                                  characteristic=self.curve.riemann_constant().characteristic,
                                  normalization="hyperelliptic", ctx=c)
        if progress:
            progress("computing fixed-point Abel lifts and second-kind coefficients")
        # Both kinds of differential follow ONE prescribed path. Independent
        # automatic endpoint calls can choose different homotopy classes at
        # a target equidistant from conjugate branch points.
        self.anchor_Y = c.sqrt(self.F(self.rho_start))
        anchor = (self.rho_start, self.anchor_Y)
        self.anchor_w = self.curve.abel_map_kind_1(anchor, reduce=False).value
        self.anchor_second = self.curve.abel_map_kind_2(anchor, reduce=False).value
        roots = sp.Poly(self.normalization["F"], rho).nroots(n=dps+5, maxsteps=1000)
        self.roots = [c.mpc(str(sp.re(root)), str(sp.im(root))) for root in roots]
        self.height = 1+max(abs(root.imag) for root in self.roots)
        self.second_forms = [self.compile(rho, sum(
            (m+1-k)*self.normalization["F"].coeff(rho, m+1+k)*rho**m/4
            for m in range(k, 2*g+2-k))) for k in range(1, g+1)]
        self.eta, self.second_at_fixed, self.fixed_paths = [], [], []
        for a, b in zip(self.a, self.b):
            first_increment, second_increment, endpoint_y = self.fixed_path(a)
            sign = 1 if abs(endpoint_y-b) < abs(endpoint_y+b) else -1
            if abs(sign*endpoint_y-b) > c.mpf(10)**(-(dps-8))*max(c.one, abs(b)):
                raise ValueError("Prescribed fixed-point path failed to reach its sheet")
            self.eta.append(sign*(self.anchor_w+first_increment))
            self.second_at_fixed.append(sign*(self.anchor_second+second_increment))
            self.fixed_paths.append(dict(sheet_sign=sign, height=c.nstr(self.height, 18),
                                         endpoint_ordinate_error=c.nstr(abs(sign*endpoint_y-b), 12)))
        self.normals = []
        self.h = []
        for eta in self.eta:
            jet = self.jet(eta, order=2)
            normal = jet[self.sharp]
            self.normals.append(normal)
            zeta = [jet[self.extend(self.sharp, k)]/normal for k in range(g)]
            h2 = zeta[1]/2
            sigma1 = jet[(1,)+(0,)*(g-1)]
            h = [zeta[0]-sigma1*h2/normal, h2]
            if g == 3:
                h.append(zeta[2])
            self.h.append(h)

    def fixed_path(self, target):
        """Transport both differential kinds along the same three segments.

        On each segment each root factor starts at sqrt(1)=1. A line from
        1 cannot meet the square-root cut unless the segment crosses a root;
        those crossings are explicitly rejected. This fixes the sheet without
        stateful square roots inside adaptive quadrature.
        """
        c, g = self.ctx, self.genus
        vertices = [c.mpc(self.rho_start), self.rho_start+c.j*self.height,
                    target+c.j*self.height, c.mpc(target)]
        first = c.matrix(g, 1)
        second = c.matrix(g, 1)
        y_start = self.anchor_Y
        for start, end in zip(vertices, vertices[1:]):
            delta = end-start
            if not delta:
                continue
            breaks = [c.zero, c.one]
            for root in self.roots:
                parameter = (root-start)/delta
                if abs(parameter.imag) < c.mpf(10)**(-(c.dps//2)) and 0 <= parameter.real <= 1:
                    raise ValueError("Chosen fixed-point contour crosses a branch point")
                if 0 < parameter.real < 1:
                    # Resolve narrow peaks near almost-collinear branch
                    # points, rather than trusting an unsplit global rule.
                    for multiplier in (-16, -4, -1, 0, 1, 4, 16):
                        node = parameter.real+multiplier*abs(parameter.imag)
                        if 0 < node < 1:
                            breaks.append(node)
            breaks = sorted(set(breaks))
            ratios = [delta/(start-root) for root in self.roots]
            ordinate = lambda parameter: y_start*c.fprod(c.sqrt(1+ratio*parameter) for ratio in ratios)
            for k in range(g):
                first[k] += c.quad(lambda parameter: delta*(start+delta*parameter)**k/ordinate(parameter), breaks)
                second[k] += c.quad(lambda parameter: delta*self.second_forms[k](start+delta*parameter)/ordinate(parameter), breaks)
            y_start = ordinate(c.one)
        return first, second, y_start

    def number(self, expression):
        return self.ctx.mpf(str(sp.N(expression, self.ctx.dps+5)))

    def compile(self, arguments, expression, **options):
        return sp.lambdify(arguments, expression, modules=self.ctx,
                           printer=MpmathPrinter({"fully_qualified_modules": False}), **options)

    @staticmethod
    def extend(index, coordinate):
        return tuple(value+int(k == coordinate) for k, value in enumerate(index))

    def jet(self, w, order=3):
        return kleinian_sigma_jet(w, order=order, **self.sigma_options)

    def log_jet_along_curve(self, jet, index, rho, Y):
        """First and second derivatives w.r.t. t=s*z, with dw_k/dt=rho^k."""
        f = jet[index]
        first = [jet[self.extend(index, k)]/f for k in range(self.genus)]
        rate = sum(rho**k*value for k, value in enumerate(first))
        acceleration = Y*sum(k*rho**(k-1)*first[k] for k in range(1, self.genus))
        for i in range(self.genus):
            for j in range(self.genus):
                hessian = jet[self.extend(self.extend(index, i), j)]/f-first[i]*first[j]
                acceleration += rho**(i+j)*hessian
        return rate, acceleration

    def evaluate(self, rho, sheet=1, w=None):
        c = self.ctx
        Y = sheet*c.sqrt(self.F(rho))
        if w is None:
            w = self.curve.abel_map_kind_1((rho, Y), reduce=False).value
        denominator = self.jet(w)
        base_rate, base_second = self.log_jet_along_curve(denominator, self.sharp, rho, Y)
        raw_u, raw_v, rates_u, rates_v, potentials_u, potentials_v = [], [], [], [], [], []
        for weight, eta, normal, h in zip(self.data["beta"], self.eta, self.normals, self.h):
            plus, minus = self.jet(w+eta), self.jet(w-eta)
            factor = c.sqrt(weight)/(denominator[self.sharp]*normal)
            raw_u.append(factor*plus[self.flat])
            raw_v.append((-1)**(self.genus+1)*factor*minus[self.flat])
            correction = sum(value*rho**k for k, value in enumerate(h))
            correction_prime = sum(k*h[k]*rho**(k-1) for k in range(1, self.genus))
            for numerator, sign, rates, potentials in (
                (plus, -1, rates_u, potentials_u), (minus, 1, rates_v, potentials_v)
            ):
                raw_rate, raw_second = self.log_jet_along_curve(numerator, self.flat, rho, Y)
                rate = self.scale*(raw_rate-base_rate+sign*correction)
                log_second = self.scale**2*(raw_second-base_second+sign*Y*correction_prime)
                rates.append(rate)
                potentials.append(log_second+rate*rate)
        inversion = (-denominator[(1, 0)]/denominator[(0, 1)] if self.genus == 2
                     else -denominator[(1, 0, 1)]/denominator[(0, 1, 1)])
        return dict(rho=rho, Y=Y, y=self.scale*Y, w=w, raw_u=raw_u, raw_v=raw_v,
                    rates_u=rates_u, rates_v=rates_v, potentials_u=potentials_u,
                    potentials_v=potentials_v, inversion=inversion,
                    sigma=denominator[(0,)*self.genus],
                    sigma3=denominator[(0, 0, 1)] if self.genus == 3 else c.zero)

    def initial(self, point):
        c = self.ctx
        U = (self.q(point["rho"])+point["y"])/2
        u = [c.one]*self.n
        u[0] = c.root(U, self.data["beta"][0])
        v = [weight*(a-point["rho"])/value
             for weight, a, value in zip(self.data["beta"], self.a, u)]
        return u+v

    def solution(self, start, end, initial):
        """Evaluate u and v from their OWN sigma ratios and opposite corrections."""
        dw = end["w"]-start["w"]
        u, v = [], []
        for j, h in enumerate(self.h):
            integral = sum(value*delta for value, delta in zip(h, dw))
            u.append(initial[j]*end["raw_u"][j]/start["raw_u"][j]*self.ctx.exp(-integral))
            v.append(initial[self.n+j]*end["raw_v"][j]/start["raw_v"][j]*self.ctx.exp(integral))
        return u+v

    def abel_increment(self, start, end, sheet):
        c = self.ctx
        return c.matrix([c.quad(lambda x: x**k/(sheet*c.sqrt(self.F(x))), [start, end])
                         for k in range(self.genus)])

    def quadrature_solution(self, start, end, initial, sheet):
        c = self.ctx
        amplitude = lambda a: c.sqrt((a-end)/(a-start))
        u, v = [], []
        for j, (a, b) in enumerate(zip(self.a, self.b_original)):
            exponent = c.quad(lambda x: b/(2*(a-x)*sheet*c.sqrt(self.F_original(x))), [start, end])
            u.append(initial[j]*amplitude(a)*c.exp(exponent))
            v.append(initial[self.n+j]*amplitude(a)*c.exp(-exponent))
        return u+v
