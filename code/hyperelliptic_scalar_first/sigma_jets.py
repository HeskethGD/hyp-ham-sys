"""Finite Taylor algebra for exact coefficient prescriptions from sigma jets.

Lists hold ordinary Taylor coefficients, not derivatives. No sampling,
finite differences, interpolation or fitted parameters are used here.
"""
from math import factorial, comb
from itertools import product


def multiply(left, right, degree):
    return [sum(left[k]*right[n-k] for k in range(n+1)
                if k < len(left) and n-k < len(right)) for n in range(degree+1)]


def divide(numerator, denominator, degree):
    if not denominator[0]:
        raise ValueError('The reference point must have nonzero sigma normalizers')
    result = []
    for n in range(degree+1):
        value = numerator[n] if n < len(numerator) else 0
        value -= sum(denominator[k]*result[n-k] for k in range(1, n+1) if k < len(denominator))
        result.append(value/denominator[0])
    return result


def compose_sigma(jet, derivative, increments, degree):
    """Taylor expansion of sigma_derivative(w_* + increments(s))."""
    g = len(increments)
    result = [0]*(degree+1)
    powers = []
    for increment in increments:
        row = [[1]+[0]*degree]
        for _ in range(degree):
            row.append(multiply(row[-1], increment, degree))
        powers.append(row)
    for index in product(range(degree+1), repeat=g):
        if sum(index) > degree:
            continue
        term = [1]+[0]*degree
        for k, order in enumerate(index):
            term = multiply(term, powers[k][order], degree)
        coefficient = jet[tuple(a+b for a, b in zip(derivative, index))]
        coefficient /= product_factorials(index)
        result = [a+coefficient*b for a, b in zip(result, term)]
    return result


def product_factorials(index):
    value = 1
    for order in index:
        value *= factorial(order)
    return value


def logarithmic_derivative(series, degree):
    derivative = [(k+1)*series[k+1] for k in range(degree+1)]
    return divide(derivative, series, degree)


def zeta_derivative(jet, base, indices):
    """Derivatives of log(sigma_base), with one to three coordinate indices."""
    def value(coordinates):
        index = list(base)
        for coordinate in coordinates:
            index[coordinate] += 1
        return jet[tuple(index)]
    f = value(())
    if len(indices) == 1:
        return value(indices)/f
    i, j = indices[:2]
    if len(indices) == 2:
        return value((i, j))/f-value((i,))*value((j,))/f**2
    if len(indices) == 3:
        k = indices[2]
        return (value((i, j, k))/f
                -(value((i, j))*value((k,))+value((i, k))*value((j,))
                  +value((j, k))*value((i,)))/f**2
                +2*value((i,))*value((j,))*value((k,))/f**3)
    raise ValueError('Only first, second and third log derivatives are needed')


def fixed_point_h(model):
    """Infinity-normalized second-kind primitives at the fixed mode points.

    For g=2,3 the primitives are short sigma expressions on the one-point
    image. This is the genus-three formula in Enolski et al., (III.10),
    with regularized infinity base point; the genus-two restriction has
    the same first two expressions. No moving reference point is used.
    """
    result = []
    for eta in model.eta:
        jet = model.jet(eta, order=2)
        normalizer = jet[model.sharp]
        zeta = [zeta_derivative(jet, model.sharp, (k,)) for k in range(model.genus)]
        sigma1 = jet[(1,)+(0,)*(model.genus-1)]
        h2 = zeta[1]/2
        h = [zeta[0]-sigma1*h2/normalizer, h2]
        if model.genus == 3:
            h.append(zeta[2])
        result.append(h)
    return result


def explicit_h_at_zero(model):
    """The expanded genus-two/three formulas displayed in the notebook."""
    c = model.ctx
    y = c.sqrt(model.F(c.zero))
    d1 = model.coefficients[1]/(2*y)
    d2 = model.coefficients[2]/y-model.coefficients[1]**2/(4*y**3)
    w = model.curve.abel_map_kind_1((c.zero, y), reduce=False).value
    order = model.genus+max(sum(model.sharp), sum(model.flat))
    denominator = model.jet(w, order)
    result = []
    for a, b, eta in zip(model.a, model.b, model.eta):
        numerator = model.jet(w+eta, order)
        Z = lambda *index: zeta_derivative(numerator, model.flat, index)-zeta_derivative(denominator, model.sharp, index)
        h = [Z(0)-(b-y)/(2*a),
             Z(1)+Z(0, 0)/y+d1/(2*a)-(b-y)/(2*a**2)]
        if model.genus == 3:
            h.append(Z(2)+3*Z(0, 1)/(2*y)
                     +Z(0, 0, 0)/(2*y**2)
                     -d1*Z(0, 0)/(2*y**2)
                     +d2/(4*a)+d1/(2*a**2)-(b-y)/(2*a**3))
        result.append(h)
    return result


def extract_h(model, reference=0, sheet=1):
    """Evaluate h from finite sigma jets at one regular reference point.

    Expand G_j in s=rho-reference to order g-1 and translate its polynomial
    coefficients back to powers of rho. The sigma derivatives are evaluated
    to order g+max(|sharp|, |flat|), and the Abel increments follow dw/drho.
    """
    c = model.ctx; reference = c.mpf(reference); g = model.genus
    shifted_F = [sum(a*comb(k, n)*reference**(k-n)
                     for k, a in enumerate(model.coefficients) if k >= n) for n in range(g)]
    y = [sheet*c.sqrt(shifted_F[0])]
    if not y[0]:
        raise ValueError('The rho reference chart must be nonbranch')
    for n in range(1, g):
        y.append((shifted_F[n]-sum(y[k]*y[n-k] for k in range(1, n)))/(2*y[0]))
    inverse_y = divide([c.one], y, g-1)
    increments = []
    for k in range(g):
        power = [c.mpf(comb(k, n))*reference**(k-n) if n <= k else c.zero for n in range(g)]
        tangent = multiply(power, inverse_y, g-1)
        increments.append([c.zero]+[tangent[n]/(n+1) for n in range(g)])
    w = model.curve.abel_map_kind_1((reference, y[0]), reduce=False).value
    order = g+max(sum(model.sharp), sum(model.flat))
    denominator = compose_sigma(model.jet(w, order), model.sharp, increments, g)
    denominator_log = logarithmic_derivative(denominator, g-1)
    coefficients = []
    for a, b, eta in zip(model.a, model.b, model.eta):
        if a == reference:
            raise ValueError('The reference point must avoid fixed mode labels')
        numerator = compose_sigma(model.jet(w+eta, order), model.flat, increments, g)
        numerator_log = logarithmic_derivative(numerator, g-1)
        raw_z_rate = multiply(y, [u-v for u, v in zip(numerator_log, denominator_log)], g-1)
        bare = divide([b-y[0]]+[-value for value in y[1:]], [2*(a-reference), -2], g-1)
        G = [u-v for u, v in zip(raw_z_rate, bare)]
        coefficients.append([sum(G[m]*comb(m, k)*(-reference)**(m-k) for m in range(k, g)) for k in range(g)])
    return coefficients
