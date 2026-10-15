"""Scalar-first polynomial Hamiltonians, in the notation of gen-int-sys.tex."""
import sympy as sp


def scalar_data(x, Q, beta, gamma, rho):
    """Restrict the actual polynomial Q, including its individual partials."""
    if not (len(x) == len(beta) == len(gamma)):
        raise ValueError("x, beta and gamma must have the same length")
    if any(not isinstance(b, int) or b <= 0 for b in beta):
        raise ValueError("This chart assumes positive integer beta_j")
    if sp.simplify(sum(gamma)) != 0:
        raise ValueError("The paper's mean-power convention requires sum(gamma)=0")
    powers = tuple(c-b*rho for b, c in zip(beta, gamma))
    line = dict(zip(x, powers))
    q = sp.expand(Q.subs(line, simultaneous=True))
    r = sp.prod(power**b for power, b in zip(powers, beta))
    F = sp.expand(q**2-4*r)
    a = tuple(sp.sympify(c)/b for c, b in zip(gamma, beta))
    ordinate = tuple(q.subs(rho, point) for point in a)
    partials = tuple(sp.expand(sp.diff(Q, variable).subs(line, simultaneous=True)) for variable in x)
    B = tuple(sp.cancel((q-value)/(2*(point-rho))-partial)
              for point, value, partial in zip(a, ordinate, partials))
    degree = sp.degree(F, rho)
    genus = (degree-1)//2
    return dict(x=x, Q=Q, beta=beta, gamma=gamma, rho=rho, powers=powers,
                q=q, r=r, F=F, a=a, b=ordinate, partials=partials, B=B,
                genus=genus, squarefree=sp.degree(sp.gcd(F, sp.diff(F, rho)), rho) == 0,
                holomorphic_gauges=all(sp.degree(value, rho) < genus for value in B))


def canonical_system(x, Q, beta):
    N = len(x)
    u = sp.symbols(f'u1:{N+1}')
    v = sp.symbols(f'v1:{N+1}')
    line = dict(zip(x, [a*b for a, b in zip(u, v)]))
    U = sp.prod(a**b for a, b in zip(u, beta))
    V = sp.prod(a**b for a, b in zip(v, beta))
    H = U+V-Q.subs(line, simultaneous=True)
    return dict(u=u, v=v, U=U, V=V, H=H,
                du=tuple(sp.diff(H, a) for a in v),
                dv=tuple(-sp.diff(H, a) for a in u))


def example():
    """Five pairs with a double mode, F leading with 4, and affine gauges."""
    rho = sp.Symbol('rho')
    x = sp.symbols('x1:6')
    beta = (2, 1, 1, 1, 1)
    gamma = tuple(map(sp.Rational, ('-3/5', '-1/5', '-1/10', '2/5', '1/2')))
    mean = -sum(x)/sum(beta)
    Q = 4*(mean+sp.Rational(27,40))*(mean+sp.Rational(1,20))*(mean-sp.Rational(3,5))
    Q += (x[0]/2-x[1]+sp.Rational(1,10))*x[2]/3
    return scalar_data(x, Q, beta, gamma, rho)


def genus_three_example():
    """Seven pairs with a double mode, F leading with 4, and quadratic gauges."""
    rho = sp.Symbol('rho')
    x = sp.symbols('x1:8')
    beta = (2, 1, 1, 1, 1, 1, 1)
    gamma = tuple(map(sp.Rational, ('-3/5', '-1/4', '-1/5', '-1/10', '1/5', '2/5', '11/20')))
    mean = -sum(x)/sum(beta)
    Q = 4*(mean+sp.Rational(37,40))*(mean+sp.Rational(3,20))*(mean-sp.Rational(3,10))*(mean-sp.Rational(13,20))
    Q += (x[0]/2-x[1]+sp.Rational(1,20))*x[2]**2/7
    return scalar_data(x, Q, beta, gamma, rho)
