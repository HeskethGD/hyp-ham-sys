"""Branch-point endpoint transport and sigma-mode jets."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "code"))

import sympy as sp
from hyperelliptic_scalar_first.special_q.model import scalar_data
from hyperelliptic_scalar_first.sigma_modes.model import SigmaModes


class BranchSigmaModes(SigmaModes):
    def fixed_path(self, target):
        """Transport both kinds together, removing the endpoint square root.

        The last segment ends at a simple root a. Set x=a-delta*t**2;
        cancel t analytically in dx/y before integrating. Other segments
        use continuously transported root factors, as in the main code.
        """
        c, g = self.ctx, self.genus
        vertices = [c.mpc(self.rho_start), self.rho_start+c.j*self.height,
                    target+c.j*self.height, c.mpc(target)]
        first, second = c.matrix(g, 1), c.matrix(g, 1)
        y_start = self.anchor_Y
        for segment, (start, end) in enumerate(zip(vertices, vertices[1:])):
            delta = end-start
            if not delta:
                continue
            if segment == 2:
                branch = min(range(len(self.roots)), key=lambda j: abs(self.roots[j]-target))
                assert abs(self.roots[branch]-target) < c.mpf(10)**(-(c.dps-8))
                ratios = [delta/(start-root) for j, root in enumerate(self.roots) if j != branch]
                def regular(t):
                    return c.fprod(c.sqrt(1+ratio*(1-t*t)) for ratio in ratios)
                def x(t):
                    return end-delta*t*t
                for k in range(g):
                    first[k] += c.quad(lambda t: 2*delta*x(t)**k/(y_start*regular(t)), [0, 1])
                    second[k] += c.quad(lambda t: 2*delta*self.second_forms[k](x(t))/(y_start*regular(t)), [0, 1])
                y_start = c.zero
            else:
                ratios = [delta/(start-root) for root in self.roots]
                def ordinate(t):
                    return y_start*c.fprod(c.sqrt(1+ratio*t) for ratio in ratios)
                for k in range(g):
                    first[k] += c.quad(lambda t: delta*(start+delta*t)**k/ordinate(t), [0, 1])
                    second[k] += c.quad(lambda t: delta*self.second_forms[k](start+delta*t)/ordinate(t), [0, 1])
                y_start = ordinate(c.one)
        return first, second, y_start


def example(genus):
    rho = sp.Symbol("rho")
    n = 2*genus+1
    variables = sp.symbols(f"I1:{n+1}")
    data = scalar_data(0, variables, (1,)*n, tuple(range(-genus, genus+1)), rho)
    data.update(name=f"branch_g{genus}", rho_start=sp.Integer(genus+1),
                rho_end=sp.Integer(genus+1)+sp.Rational(1, 50))
    return data


def quotient_jet(model, w, rho, Y, eta, normal, h, sign):
    """Value and two propagation derivatives without division by a mode.

    The moving denominator sigma_sharp is nonzero here. The shifted
    numerator may vanish, notably at the fixed branch point itself.
    """
    c, g = model.ctx, model.genus
    denominator = model.jet(w)
    numerator = model.jet(w+sign*eta)
    f, d = numerator[model.flat], denominator[model.sharp]
    def derivatives(jet, index):
        first = sum(rho**k*jet[model.extend(index, k)] for k in range(g))
        second = Y*sum(k*rho**(k-1)*jet[model.extend(index, k)] for k in range(1, g))
        second += sum(rho**(i+j)*jet[model.extend(model.extend(index, i), j)]
                      for i in range(g) for j in range(g))
        return first, second
    f1, f2 = derivatives(numerator, model.flat)
    d1, d2 = derivatives(denominator, model.sharp)
    raw = f/(d*normal)
    raw1 = (f1-f*d1/d)/(d*normal)
    raw2 = (f2-2*f1*d1/d+f*(2*(d1/d)**2-d2/d))/(d*normal)
    correction = -sign*sum(h[k]*rho**k for k in range(g))
    correction1 = -sign*Y*sum(k*h[k]*rho**(k-1) for k in range(1, g))
    # Initial-value normalization cancels the omitted fixed exponential.
    factor = c.exp(-sign*sum(h[k]*w[k] for k in range(g)))
    parity = 1 if sign == 1 else (-1)**(g+1)
    return tuple(parity*factor*v for v in (
        raw, model.scale*(raw1+correction*raw),
        model.scale**2*(raw2+2*correction*raw1+(correction1+correction**2)*raw)))


