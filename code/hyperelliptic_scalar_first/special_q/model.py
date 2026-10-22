"""Exact algebra for special polynomial Hamiltonians and mode potentials."""

from itertools import product

import sympy as sp


def validate_weights(beta):
    beta = tuple(beta)
    if not beta or any(not isinstance(b, int) or isinstance(b, bool) or b <= 0 for b in beta):
        raise ValueError("Weights must be a nonempty tuple of positive integers")
    return beta


def special_polynomial(variables, beta):
    """H_beta in the original, unnormalised power variables.

    Return None for odd total weight: there is no nonconstant solution.
    Truncated convolution avoids expanding an infinite formal power series.
    """
    beta = validate_weights(beta)
    if len(variables) != len(beta):
        raise ValueError("One independent power variable is required per weight")
    if sum(beta) % 2:
        return None
    m = sum(beta) // 2
    coefficients = [sp.Integer(1)] + [sp.Integer(0)] * m
    for variable, weight in zip(variables, beta):
        factor = [sp.binomial(sp.Rational(weight, 2), k)
                  * (-variable / sp.Integer(weight))**k for k in range(m + 1)]
        coefficients = [sp.expand(sum(coefficients[k] * factor[n-k]
                                     for k in range(n + 1))) for n in range(m + 1)]
    return coefficients[m]


def cancellation_residuals(Q, variables, beta):
    """Denominator-cleared universal identities, in original power variables."""
    beta = validate_weights(beta)
    Q = sp.sympify(Q)
    return tuple(sp.expand(Q - Q.xreplace(
        {x: x - sp.Rational(b, beta[j]) * variables[j]
         for x, b in zip(variables, beta)})
        - 2 * variables[j] * sp.diff(Q, variables[j]) / beta[j])
        for j in range(len(beta)))


def homogeneous_kernel(beta, degree):
    """Solve coefficient constraints without using the proposed H_beta formula.

    This finite calculation checks completeness only at the specified degree;
    the arbitrary-degree argument is in proof_scope.md and the notebook.
    """
    beta = validate_weights(beta)
    variables = sp.symbols(f"I1:{len(beta)+1}")
    powers = [p for p in product(range(degree + 1), repeat=len(beta)) if sum(p) == degree]
    monomials = [sp.prod(x**k for x, k in zip(variables, p)) for p in powers]
    columns = []
    for monomial in monomials:
        column = []
        for residual in cancellation_residuals(monomial, variables, beta):
            polynomial = sp.Poly(residual, *variables)
            column.extend(polynomial.coeff_monomial(p) for p in powers)
        columns.append(column)
    matrix = sp.Matrix.hstack(*(sp.Matrix(c) for c in columns))
    return variables, monomials, matrix.nullspace()


def scalar_data(Q, variables, beta, gamma, rho):
    """Differentiate Q BEFORE restricting to its invariant power line."""
    beta = validate_weights(beta)
    if not len(variables) == len(beta) == len(gamma):
        raise ValueError("Variables, weights and invariants must have equal lengths")
    if sp.simplify(sum(gamma)) != 0:
        raise ValueError("The paper's scalar convention requires sum(gamma)=0")
    Q = sp.sympify(Q)
    powers = tuple(c - b*rho for c, b in zip(gamma, beta))
    line = dict(zip(variables, powers))
    q = sp.expand(Q.subs(line, simultaneous=True))
    r = sp.prod(power**b for power, b in zip(powers, beta))
    F = sp.Poly(sp.expand(q*q - 4*r), rho)
    a = tuple(sp.sympify(c)/b for c, b in zip(gamma, beta))
    ordinate = tuple(q.subs(rho, point) for point in a)
    partials = tuple(sp.diff(Q, x).subs(line, simultaneous=True) for x in variables)
    corrections = tuple(sp.cancel((q-value)/(2*(point-rho)) - partial)
                        for point, value, partial in zip(a, ordinate, partials))
    smooth = not F.is_zero and sp.degree(sp.gcd(F, F.diff())) == 0
    degree = None if F.is_zero else int(F.degree())
    return dict(Q=Q, variables=variables, beta=beta, gamma=gamma, rho=rho,
                powers=powers, q=q, r=r, F=F.as_expr(), a=a, ordinate=ordinate,
                partials=partials, B=corrections, degree=degree,
                squarefree=bool(smooth),
                genus=(degree-1)//2 if smooth and degree >= 3 else None)


def canonical_field(data):
    """Differentiate the full Hamiltonian in independent canonical variables."""
    beta = data["beta"]
    u = sp.symbols(f"u1:{len(beta)+1}")
    v = sp.symbols(f"v1:{len(beta)+1}")
    U = sp.prod(x**b for x, b in zip(u, beta))
    V = sp.prod(x**b for x, b in zip(v, beta))
    Q = data["Q"].subs(dict(zip(data["variables"], [x*y for x, y in zip(u, v)])),
                        simultaneous=True)
    H = U + V - Q
    return u + v, tuple(sp.diff(H, x) for x in v) + tuple(-sp.diff(H, x) for x in u), H


def example_models():
    """Smooth odd-degree examples; leading coefficient need not yet equal four."""
    rho = sp.Symbol("rho")
    specifications = [
        ("constant_weighted_g2", (2, 1, 1, 1), (-3, -1, 1, 3), 1, 4),
        ("constant_weighted_g3", (2, 1, 1, 1, 1, 1), (-5, -3, -1, 1, 3, 5), 1, 6),
        ("nonconstant_weighted_g2", (5, 5), (-1, 1), None, 1),
    ]
    models = []
    for name, beta, gamma, Q, start in specifications:
        variables = sp.symbols(f"I1:{len(beta)+1}")
        if Q is None:
            # C=prod beta_i**beta_i=5**10; kappa=2 sqrt(C).
            Q = 6250 * special_polynomial(variables, beta) + 1
        data = scalar_data(Q, variables, beta, gamma, rho)
        data.update(name=name, rho_start=sp.Integer(start), rho_end=sp.Integer(start)+sp.Rational(1, 50))
        models.append(data)
    return models
