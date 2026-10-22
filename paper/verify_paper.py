"""Check printed example algebra and numerical tables against the notebook.

Run from the project root: .venv/bin/python paper/verify_paper.py
The small LaTeX reader accepts only the polynomial notation used in this paper;
it does not interpret arbitrary LaTeX or evaluate Python code.
"""
from pathlib import Path
import json
import re
import sys

import sympy as sp
from sympy.parsing.sympy_parser import (
    parse_expr, standard_transformations, implicit_multiplication_application,
)

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'code'))
from hyperelliptic_scalar_first.model import example, genus_three_example
from hyperelliptic_scalar_first.verification import symbolic_checks

SOURCE = (ROOT / 'paper/main.tex').read_text()


def group(text, pos):
    if text[pos] != '{':
        return text[pos], pos + 1
    depth = 1
    end = pos + 1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}')
        end += 1
    return text[pos + 1:end - 1], end


def fractions(text):
    result = ''
    while text:
        if text.startswith(r'\frac'):
            numerator, pos = group(text, 5)
            denominator, end = group(text, pos)
            result += f'(({fractions(numerator)})/({fractions(denominator)}))'
            text = text[end:]
        else:
            result += text[0]
            text = text[1:]
    return result


def polynomial(text, symbols):
    text = re.sub(r'\\(?:begin|end)\{aligned\}', '', text)
    for command in (r'\left', r'\right', r'\quad', r'\\', '&', '{}', '$'):
        text = text.replace(command, '')
    text = text.replace(r'\times', '*').replace(r'\rho', 'rho')
    text = text.replace('S(x)', 'S').replace('A_2(x)', 'A2').replace('A_3(x)', 'A3')
    text = re.sub(r'x_(\d)', r'x\1', text)
    text = fractions(text).replace('^', '**').replace('{', '(').replace('}', ')')
    text = text.rstrip('., \n').replace('\n', ' ')
    assert re.fullmatch(r'[A-Za-z0-9+*/().\s-]+', text), text
    return parse_expr(text, local_dict=symbols,
                      transformations=standard_transformations +
                      (implicit_multiplication_application,))


def equation(label):
    start = r'\begin{equation}\label{' + label + '}'
    return SOURCE.split(start, 1)[1].split(r'\end{equation}', 1)[0]


def equal(actual, expected):
    assert sp.cancel(actual - expected) == 0, (actual, expected)


def check_special_examples():
    from hyperelliptic_scalar_first.special_q.model import example_models

    models = example_models()
    for label, data, expected in (
        ('eq:g2-F', models[0], models[0]['F']),
        ('eq:g3-F', models[1], models[1]['F']),
        ('eq:nonconstant-q', models[2], models[2]['q']),
        ('eq:nonconstant-F', models[2], models[2]['F']),
    ):
        actual = polynomial(equation(label).split('=', 1)[1], {'rho': data['rho']})
        equal(actual, expected)
        if label != 'eq:nonconstant-q':
            assert sp.degree(sp.gcd(expected, sp.diff(expected, data['rho'])), data['rho']) == 0
    data = models[0]
    rho = data['rho']
    quotient, remainder = sp.div(data['F'] - data['ordinate'][0]**2,
                                 rho - data['a'][0], rho)
    assert remainder == 0
    equal(polynomial(equation('eq:g2-potential').split('=', 1)[1], {'rho': rho}),
          sp.diff(quotient, rho)/4)
    power = sp.Symbol('I1')
    self_coefficient = equation('eq:special-g2-self').split('=', 1)[1].split(r'\psi_1', 1)[0]
    equal(polynomial(self_coefficient.replace('I_1', 'I1'), {'I1': power}),
          (sp.diff(quotient, rho)/4).subs(rho, data['a'][0] - power/data['beta'][0]))
    data = models[1]
    rho = data['rho']
    quotient, remainder = sp.div(data['F'] - data['ordinate'][0]**2,
                                 rho - data['a'][0], rho)
    assert remainder == 0
    equal(polynomial(equation('eq:special-g3-potential').split('=', 1)[1], {'rho': rho}),
          sp.diff(quotient, rho)/4)
    data = models[2]
    rho = data['rho']
    body = re.sub(r'\\(?:begin|end)\{aligned\}', '', equation('eq:special-nonconstant-potentials'))
    rows = body.split(r'\\')
    assert len(rows) == 2
    for j, row in enumerate(rows):
        quotient, remainder = sp.div(data['F'] - data['ordinate'][j]**2,
                                     rho - data['a'][j], rho)
        assert remainder == 0
        equal(polynomial(row.split('=', 1)[1], {'rho': rho}), sp.diff(quotient, rho)/4)
    algebra = json.loads((ROOT / 'code/hyperelliptic_scalar_first/special_q/results.json').read_text())
    sigma = json.loads((ROOT / 'code/hyperelliptic_scalar_first/sigma_modes/results.json').read_text())
    assert algebra['exact']['all_passed'] and algebra['numerical']['all_passed']
    assert algebra['exact']['check_count'] == 154
    assert sigma['all_passed'] and sigma['check_count'] == 680
    assert sum(len(row['exact_normalization_checks']) for row in sigma['examples']) == 30
    bounds = (sp.Rational('1.25e-33'), sp.Rational('1.75e-31'), sp.Rational('9.27e-33'))
    for row, bound in zip(sigma['examples'], bounds):
        assert sp.Rational(row['maximum_scaled_error']) <= bound


def check_examples():
    for tag, data in (('g2', example()), ('g3', genus_three_example())):
        rho = data['rho']
        symbols = {'rho': rho, 'S': sum(data['x'])}
        symbols.update({str(x): x for x in data['x']})
        Q = equation(f'ex:{tag}-Q').split('&=', 1)[1].split(r'\qquad', 1)[0]
        equal(polynomial(Q, symbols), data['Q'])
        genus = data['genus']
        common = equation(f'ex:{tag}-common-phase').split('=', 1)[1].split(r'\qquad', 1)[0]
        common = polynomial(common, symbols)
        equal(common, sp.diff(data['q'], rho).subs(rho, -symbols['S'] / sum(data['beta'])) / sum(data['beta']))
        symbols[f'A{genus}'] = common
        phases = equation(f'ex:{tag}-phases')
        for j, variable in enumerate(data['x'], 1):
            if j <= 3:
                phase = phases.split(rf'\Phi_{j}^{{({genus})}}(x)&=', 1)[1].split(r'\\', 1)[0]
                actual = polynomial(phase, symbols)
            else:
                actual = common
            equal(actual, -sp.diff(data['Q'], variable))
        system = equation(f'ex:{tag}-coupled')
        assert rf'u_1\prod_{{k=2}}^{len(data["x"])}u_k' in system
        assert rf'+2v_1\prod_{{k=2}}^{len(data["x"])}v_k' in system
        assert r'+v_1^2\prod_' in system and r'-u_1^2\prod_' in system
        for key in ('q', 'r', 'F'):
            text = equation(f'ex:{tag}-{key}').split('=', 1)[1]
            equal(polynomial(text, symbols), data[key])
        if tag == 'g2':
            for j, expected in enumerate(data['B'], 1):
                text = SOURCE.split(rf'B_{j}(\rho)&=', 1)[1].split(r'\label', 1)[0]
                equal(polynomial(text, symbols), expected)
        else:
            table = SOURCE.split('$j$ & $c_{j1}$ & $c_{j2}$ & $c_{j3}$', 1)[1]
            table = table.split(r'\bottomrule', 1)[0]
            for j, expected in enumerate(data['B'], 1):
                row = re.search(rf'^\s*{j} & (.*?)\\\\', table, re.M).group(1)
                coefficients = [polynomial(cell, symbols) for cell in row.split('&')]
                equal(sum(c * rho**k for k, c in enumerate(coefficients)), expected)


def check_results():
    report = json.loads((ROOT / 'code/hyperelliptic_scalar_first/verification_results.json').read_text())
    assert report['passed']
    assert len(report['symbolic']['checks']) == 46
    for key, count in (('numerical', 441), ('genus_three', 637)):
        assert len(report[key]['rows']) == count
    extended = json.loads((ROOT / 'code/hyperelliptic_scalar_first/genus_two_example/results.json').read_text())
    assert extended['passed'] and extended['dps'] == 30
    assert extended['analytic_samples'] == 33 and extended['plot_samples'] == 513
    assert extended['z_duration'] in SOURCE
    for row in extended['convergence']:
        assert f"{row['steps']:,}" in SOURCE
        mantissa, exponent = f"{float(row['max_mode_error']):.7e}".split('e')
        assert mantissa + r'\times10^{' + str(int(exponent)) + '}' in SOURCE
    for artifact in ('modes_combined.png',):
        assert 'genus_two_example/' + artifact in SOURCE
        assert (ROOT / 'code/hyperelliptic_scalar_first/genus_two_example' / artifact).is_file()
    oscillatory = json.loads((ROOT / 'code/hyperelliptic_scalar_first/genus_three_example/results.json').read_text())
    assert oscillatory['passed'] and oscillatory['dps'] == 30 and oscillatory['cycles'] == 3
    assert oscillatory['analytic_samples'] == 97 and oscillatory['plot_samples'] == 1537
    assert oscillatory['turning_and_zero_observations'] == 18
    assert oscillatory['z_duration'] in SOURCE and oscillatory['scalar_period'] in SOURCE
    for example in (extended, oscillatory):
        for key in ('sigma_inversion_error', 'paired_power_error', 'U_product_error', 'V_product_error'):
            mantissa, exponent = f"{float(example['analytic_maxima'][key]):.6e}".split('e')
            assert mantissa + r'\times10^{' + str(int(exponent)) + '}' in SOURCE
    for row in oscillatory['convergence']:
        assert f"{row['steps']:,}" in SOURCE
        for key in ('step', 'max_mode_error'):
            mantissa, exponent = f"{float(row[key]):.7e}".split('e')
            assert mantissa + r'\times10^{' + str(int(exponent)) + '}' in SOURCE
    for artifact in ('scalar_powers_combined.png', 'modes_combined.png'):
        assert 'genus_three_example/' + artifact in SOURCE
        assert (ROOT / 'code/hyperelliptic_scalar_first/genus_three_example' / artifact).is_file()
    determinants = json.loads((ROOT / 'code/hyperelliptic_scalar_first/onishi_products/results.json').read_text())
    assert determinants['passed'] and determinants['dps'] == 30
    assert len(determinants['rows']) == 180 and all(row['passed'] for row in determinants['rows'])
    assert sum(len(row['checks']) for row in determinants['exact']) == 40
    assert len(determinants['lower_family_interpolation']) == 4
    for row in determinants['exact']:
        for key in ('Vandermonde_inverse_recovers_q', 'Schur_block_determinant', 'Onishi_Schur_column_sign'):
            assert row['checks'][key]
    assert r'\mathsf W_A^{-1}' in equation('eq:onishi-vandermonde-q')
    assert r'\varepsilon D_n(z)' in equation('eq:onishi-schur-block')
    assert all(all(row['checks'].values()) for row in determinants['exact'])
    assert r"\rho'/2" in SOURCE
    assert r'\prod_{1\le i<j\le4}' in equation('eq:onishi-g2-full-sigma')
    assert r'\sigma_2(w(z))^5' in equation('eq:onishi-g2-full-sigma')
    assert r'\sum_{j=1}^{N}\beta_j(c_{jk}-h_{jk})=0' in equation('eq:onishi-weighted-gauge')
    assert r'U(z)=\varepsilon\frac{D_n^*(z)}{\Delta_A}' in equation('eq:onishi-mode-products')
    assert '320 & ' not in SOURCE and '1280 & ' not in SOURCE
    assert 'numerical checks do not test passage through them' not in SOURCE
    coordinate_report = json.loads((ROOT / 'paper/scalar_inversion_checks.json').read_text())
    assert coordinate_report['passed'] and coordinate_report['dps'] == 25
    assert len(coordinate_report['rows']) == 8
    largest = max(float(row[field]) for row in coordinate_report['rows'] if row['genus'] == 4
                  for field in ('general_scaled_error', 'published_scaled_error'))
    assert abs(largest - 3.304667643e-22) < 1e-30


if __name__ == '__main__':
    check_examples()
    check_special_examples()
    check_results()
    assert symbolic_checks()['passed']
    assert r"\frac{\beta_jV}{I_j}-Q_j(I)" in equation('eq:mode-log')
    assert r"\frac{q(\rho)-\rho'}{2[\gamma_j/\beta_j-\rho]}" in equation('eq:mode-log')
    assert r"-\frac{\beta_jU}{I_j}+Q_j(I)" in equation('eq:mode-log')
    assert r'\beta_j^2' in SOURCE.split(r'\begin{proof}', 1)[1]
    assert r"\beta_j\rho'(\mu_j)" not in SOURCE
    assert r'\boxed' not in SOURCE and 'R_j^' not in SOURCE
    assert r'\bar u_j(z)' in SOURCE and r'\bar v_j(z)' in SOURCE
    assert r'\sigma_{144}(w(z))' in SOURCE and r'\sigma_{244}(w(z))' in SOURCE
    for body in re.split(r'\\section\{[^}]+\}(?:\\label\{[^}]+\})?', SOURCE)[1:]:
        intro = body.split(r'\subsection', 1)[0].strip()
        assert intro and not intro.startswith(r'\subsection')
    assert 'enolski-g3' not in SOURCE
    assert r'\rho(\mu_j)' in SOURCE and r"\rho'(\mu_j)" in SOURCE
    labels = re.findall(r'\\label\{([^}]+)\}', SOURCE)
    assert len(labels) == len(set(labels))
    bibkeys = set(re.findall(r'\\bibitem\{([^}]+)\}', SOURCE))
    cited = set()
    for group_text in re.findall(r'\\cite(?:\[[^]]*\])?\{([^}]+)\}', SOURCE):
        cited.update(group_text.split(','))
    assert cited == bibkeys, (cited - bibkeys, bibkeys - cited)
    print('PASS: printed Hamiltonians, coupled-system polynomials, curve and gauge data, numerical tables, special-polynomial examples and saved sigma checks, inversion report and references')
