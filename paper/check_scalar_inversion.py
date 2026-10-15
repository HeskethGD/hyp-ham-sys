"""Targeted checks of the new scalar coordinate ratios, separate from the notebook.

Run from the project root: .venv/bin/python paper/check_scalar_inversion.py
The genus-four curve is a coordinate-recovery example, not a canonical mode example.
"""
import json
from pathlib import Path
import sys

import sympy as sp
from mpmath import mp
from generapy import Curve, kleinian_sigma_jet

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'code'))
from hyperelliptic_scalar_first.model import genus_three_example


def check():
    reports = []
    for genus in (3, 4):
        ctx = mp.clone()
        ctx.dps = 25
        rho = sp.Symbol('rho')
        F = (genus_three_example()['F'] if genus == 3
             else 4 * sp.prod(rho - k for k in range(1, 10)))
        coefficients = [ctx.mpf(str(F.expand().coeff(rho, k)))
                        for k in range(2 * genus + 2)]
        curve = Curve({(0, 2): 1} | {(k, 0): -v for k, v in enumerate(coefficients) if v}, ctx=ctx)
        first = curve.periods_kind_1()
        second = curve.periods_kind_2()
        options = dict(omega=first.omega, tau=first.tau, kappa=second.kappa,
                       characteristic=curve.riemann_constant().characteristic,
                       normalization='hyperelliptic', ctx=ctx)
        for sheet in (1, -1):
            for coordinate in (('5.5', '5.6') if genus == 3 else ('10', '10.2')):
                value = ctx.mpf(coordinate)
                ordinate = sheet * ctx.sqrt(sum(v * value**k for k, v in enumerate(coefficients)))
                w = curve.abel_map_kind_1((value, ordinate), reduce=False).value
                jet = kleinian_sigma_jet(w, order=3, **options)
                # D_g takes derivatives in the even coordinates 4,6,... .
                base = [int(k >= 3 and (k + 1) % 2 == 0) for k in range(genus)]
                numerator, denominator = base.copy(), base.copy()
                numerator[0] += 1
                denominator[1] += 1
                general = -jet[tuple(numerator)] / jet[tuple(denominator)]
                numerator, denominator = [0] * genus, [0] * genus
                numerator[0], denominator[1] = 1, 1
                numerator[-1] = denominator[-1] = genus - 2
                published = -jet[tuple(numerator)] / jet[tuple(denominator)]
                errors = [abs(answer - value) / max(ctx.one, abs(value))
                          for answer in (general, published)]
                assert max(errors) < ctx.mpf('1e-18'), (genus, sheet, coordinate, errors)
                reports.append(dict(genus=genus, sheet=sheet, rho=coordinate,
                                    general_scaled_error=ctx.nstr(errors[0], 10),
                                    published_scaled_error=ctx.nstr(errors[1], 10)))
        print(f'Genus {genus}: targeted coordinate checks passed', flush=True)
    return dict(passed=True, dps=25, tolerance='1e-18', rows=reports,
                scope='Coordinate recovery at two regular points on each sheet in genera three and four; no genus-four canonical trajectory check.')


if __name__ == '__main__':
    report = check()
    target = Path(__file__).with_name('scalar_inversion_checks.json')
    target.write_text(json.dumps(report, indent=2) + '\n')
    print('Report:', target)
