import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import sys as _sys
    from pathlib import Path

    example_dir = Path(__file__).resolve().parent
    _root = example_dir.parent.parent
    if str(_root) not in _sys.path:
        _sys.path.insert(0, str(_root))

    import json
    import importlib.metadata as metadata
    from time import perf_counter
    import marimo as mo
    import sympy as sp
    from mpmath import mp
    import matplotlib.pyplot as plt
    from generapy import Curve, kleinian_sigma_jet

    return Curve, example_dir, json, kleinian_sigma_jet, metadata, mo, mp, perf_counter, plt, sp


@app.cell
def _(mo):
    mo.md(r"""
    # Genus two: from the polynomial Hamiltonian to five mode trajectories

    This is the **first worked example**, isolated for review. Every formula,
    generapy call, scalar inversion step and RK4 step used here appears in this
    notebook. We retain the five pairs, the double weight, and the polynomial
    Hamiltonian of the paper. We choose a symmetric set of initial mode
    amplitudes and propagate on the positive sheet from $\rho=4.5$ to $\rho=10$.

    The calculation follows three levels:

    | Level | What happens here | Is anything fitted? |
    |---|---|---|
    | Exact algebra | Define $Q$ with rational coefficients; derive $q,r,F$, $\partial Q/\partial x_j$, and $B_j$ with SymPy. | No. |
    | Analytic formulas, evaluated numerically | Use the published sigma ratios and fixed-point expressions for $h_{jk}$. Generapy evaluates periods, Abel integrals, sigma and sigma derivatives. | No coefficients are inferred from a trajectory. |
    | Numerical propagation and comparison | Invert a known scalar quadrature; independently integrate the full ten-component polynomial field with RK4. | Root finding determines the curve coordinate for a specified $z$; it does not fit a model. |

    **30 dps is the arithmetic precision, not a promise that finite-step RK4
    has 30 correct digits.** We measure step-halving convergence rather than
    assume the two accuracies are the same. Sigma evaluations are sparse;
    the RK4 integration grid is much finer. All displayed modes are actual
    amplitudes, without rescaling each curve to make agreement look better.

    The original combined notebook, the second example, and the paper are
    independent of this review notebook. This path is regular: it does not
    test turning-point passage or global continuation around periods.
    """)
    return


@app.cell
def _(mo, mp):
    # These are numerical experiment settings, not fitted Hamiltonian parameters.
    dps = 30
    ctx = mp.clone()
    ctx.dps = dps
    rho_start = ctx.mpf("4.5")
    rho_end = ctx.mpf("10")
    analytic_samples = 33
    plot_samples = 513
    rk4_steps = (32768, 65536, 131072)
    quadrature_orders = (64, 96)
    mo.md(f"""
    ## 1. Choose the experiment settings

    We evaluate the analytic modes at **{analytic_samples}** propagation values.
    The finest RK4 run takes **{rk4_steps[-1]:,}** steps, storing only the
    analytic comparison nodes and **{plot_samples}** points for the dashed lines.
    The two other grids check convergence. Edit the constants in this cell to
    change the experiment; no analytic evaluation is performed at every RK4 step.
    """)
    return analytic_samples, ctx, dps, plot_samples, quadrature_orders, rho_end, rho_start, rk4_steps


@app.cell
def _(mo):
    mo.md(r"""
    ## 2. Specify the exact Hamiltonian

    This is input to the construction, not output from generapy. The variables
    $x_j$ are independent powers; only after differentiating $Q$ do we restrict
    them to the scalar line. Put $S(x)=\sum_{j=1}^5x_j$ and
    \[
    Q=4\left(-\frac{S(x)}6+\frac{27}{40}\right)
    \left(-\frac{S(x)}6+\frac1{20}\right)
    \left(-\frac{S(x)}6-\frac35\right)
    +\frac13\left(\frac{x_1}2-x_2+\frac1{10}\right)x_3.
    \]
    The weights are $(2,1,1,1,1)$ and the offsets are
    $(-3/5,-1/5,-1/10,2/5,1/2)$. The canonical Hamiltonian is
    $H=u_1^2u_2u_3u_4u_5+v_1^2v_2v_3v_4v_5-Q(u_1v_1,\ldots,u_5v_5)$.
    The term vanishing on the power line still changes the individual modes.
    """)
    return


@app.cell
def _(mo, sp):
    rho = sp.Symbol("rho")
    x = sp.symbols("x1:6")
    beta = (2, 1, 1, 1, 1)
    gamma = tuple(map(sp.Rational, ("-3/5", "-1/5", "-1/10", "2/5", "1/2")))
    S = sum(x)
    Q = 4*(-S/6+sp.Rational(27, 40))*(-S/6+sp.Rational(1, 20))*(-S/6-sp.Rational(3, 5))
    Q += (x[0]/2-x[1]+sp.Rational(1, 10))*x[2]/3
    powers = tuple(offset-weight*rho for weight, offset in zip(beta, gamma))
    power_line = dict(zip(x, powers))
    q = sp.expand(Q.subs(power_line, simultaneous=True))
    r = sp.expand(sp.prod(power**weight for power, weight in zip(powers, beta)))
    F = sp.expand(q*q-4*r)
    Q_partials = tuple(sp.diff(Q, variable) for variable in x)
    restricted_partials = tuple(sp.expand(value.subs(power_line, simultaneous=True)) for value in Q_partials)
    mode_coordinates = tuple(offset/weight for offset, weight in zip(gamma, beta))
    mode_ordinates = tuple(q.subs(rho, coordinate) for coordinate in mode_coordinates)
    B = tuple(sp.cancel((q-value)/(2*(coordinate-rho))-partial)
              for coordinate, value, partial in zip(mode_coordinates, mode_ordinates, restricted_partials))

    exact_checks = {
        "centred offsets": sum(gamma) == 0,
        "degree five": sp.degree(F, rho) == 5,
        "leading coefficient four": sp.LC(sp.Poly(F, rho)) == 4,
        "squarefree": sp.degree(sp.gcd(F, sp.diff(F, rho)), rho) == 0,
        "holomorphic gauges": all(sp.degree(value, rho) <= 1 for value in B),
        "fixed points nonbranch": all(value != 0 for value in mode_ordinates),
    }
    assert all(exact_checks.values())
    mo.md(r"The exact restriction gives \[F(\rho)="+sp.latex(F)+r".\]")
    return B, F, Q, Q_partials, beta, exact_checks, gamma, mode_coordinates, mode_ordinates, q, r, rho, x


@app.cell
def _(B, Q_partials, mo, rho, sp):
    mo.ui.table([
        {"mode": j+1, "full partial of Q": str(partial), "B_j(rho)": str(value),
         "c_j1 (exact)": str(value.coeff(rho, 0)), "c_j2 (exact)": str(value.coeff(rho, 1))}
        for j, (partial, value) in enumerate(zip(Q_partials, B))
    ])
    return


@app.cell
def _(F, Q, Q_partials, beta, mo, rho, sp, x):
    # Differentiate the full polynomial Hamiltonian to check the field below.
    exact_u = sp.symbols("u1:6")
    exact_v = sp.symbols("v1:6")
    exact_I = tuple(a*b for a, b in zip(exact_u, exact_v))
    exact_H = sp.prod(a**b for a, b in zip(exact_u, beta)) + sp.prod(a**b for a, b in zip(exact_v, beta))
    exact_H -= Q.subs(dict(zip(x, exact_I)), simultaneous=True)
    exact_field = tuple(sp.diff(exact_H, value) for value in exact_v) + tuple(-sp.diff(exact_H, value) for value in exact_u)
    phi_common = sum(x)**2/18-sum(x)/36-sp.Rational(107, 400)
    phi = (phi_common-x[2]/6, phi_common+x[2]/3,
           phi_common-(x[0]/2-x[1]+sp.Rational(1, 10))/3, phi_common, phi_common)
    assert all(sp.expand(value+partial) == 0 for value, partial in zip(phi, Q_partials))
    largest_real_root_interval = sp.Poly(F, rho).intervals(eps=sp.Rational(1, 10**12))[-1][0]
    mo.md(r"""
    ## 3. Read off the complete polynomial field

    The common multiplicative term is $S(x)^2/18-S(x)/36-107/400$.
    The first three modes add $-x_3/6$, $+x_3/3$, and
    $-(x_1/2-x_2+1/10)/3$, respectively. These equal $-\partial_{x_j}Q$
    exactly, even off the power line. The double mode contributes
    $2v_1v_2v_3v_4v_5$ to $u_1'$ and the analogous negative product to $v_1'$;
    unit-weight modes contain $v_1^2$ or $u_1^2$ instead.

    The symbolic check differentiates $H$ directly. This is useful because the
    independent RK4 integration must use the full field, including deviations
    from the power line, rather than substitute the analytic scalar trajectory.
    """)
    return exact_field, exact_u, exact_v, largest_real_root_interval, phi


@app.cell
def _(B, F, ctx, gamma, mode_coordinates, mode_ordinates, q, rho, sp):
    # Convert exact coefficients to arbitrary-precision numbers, never via float.
    def _number(value):
        _numerator, _denominator = sp.fraction(value)
        return ctx.mpf(str(_numerator))/ctx.mpf(str(_denominator))

    F_coefficients = tuple(_number(F.coeff(rho, k)) for k in range(6))
    q_coefficients = tuple(_number(q.coeff(rho, k)) for k in range(4))
    gauge_coefficients = tuple(tuple(_number(value.coeff(rho, k)) for k in range(2)) for value in B)
    numeric_gamma = tuple(_number(value) for value in gamma)
    fixed_rho = tuple(_number(value) for value in mode_coordinates)
    fixed_ordinate = tuple(_number(value) for value in mode_ordinates)

    def F_numeric(coordinate):
        return sum(value*coordinate**k for k, value in enumerate(F_coefficients))

    def q_numeric(coordinate):
        return sum(value*coordinate**k for k, value in enumerate(q_coefficients))

    return F_coefficients, F_numeric, fixed_ordinate, fixed_rho, gauge_coefficients, numeric_gamma, q_numeric


@app.cell
def _(Curve, F_coefficients, ctx, kleinian_sigma_jet, mo, perf_counter):
    mo.md(r"""
    ## 4. What generapy supplies

    We now give the **already derived curve** to generapy. Its calls supply
    first- and second-kind period data, the characteristic of the Riemann
    constant, forward Abel integrals from infinity, and numerical sigma jets
    (sigma and partial derivatives). The prescribed normalization has leading
    polynomial $w_1-w_2^3/3$.

    The curve, Hamiltonian, mode weights, correction polynomials and equations
    for the amplitudes were supplied or derived above. They are not inferred by
    a generapy call. In this workflow we also implement scalar **inverse**
    quadrature, the exponential coefficient formula and RK4 ourselves. A
    forward Abel-map call takes a known curve point; it does not itself produce
    $\rho(z)$ from a prescribed propagation value.

    Analytic sigma identities remain analytic identities even though their
    periods, integrals and function values are computed numerically. Neither
    theta-series evaluation nor numerical quadrature is coefficient fitting.
    """)
    _started = perf_counter()
    curve = Curve({(0, 2): 1} | {(k, 0): -value for k, value in enumerate(F_coefficients) if value}, ctx=ctx)
    first_periods = curve.periods_kind_1()
    second_periods = curve.periods_kind_2()
    characteristic = curve.riemann_constant().characteristic
    sigma_options = dict(omega=first_periods.omega, tau=first_periods.tau,
                         kappa=second_periods.kappa, characteristic=characteristic,
                         normalization="hyperelliptic", ctx=ctx)

    def sigma_jet(abel_vector, order=2):
        return kleinian_sigma_jet(abel_vector, order=order, **sigma_options)

    period_seconds = perf_counter()-_started
    return curve, first_periods, period_seconds, second_periods, sigma_jet


@app.cell
def _(ctx, curve, fixed_ordinate, fixed_rho, mo, sigma_jet):
    mo.md(r"""
    ## 5. Evaluate the fixed-point coefficients without fitting

    The fixed point for mode $j$ is
    $(\rho(\mu_j),\rho'(\mu_j))=(\gamma_j/\beta_j,q(\gamma_j/\beta_j))$.
    We calculate its unreduced Abel lift with generapy and insert sigma
    derivatives into the **specified formula**
    \[
    h_{j2}=\frac{\sigma_{22}(w(\mu_j))}{2\sigma_2(w(\mu_j))},\qquad
    h_{j1}=\frac{\sigma_{12}(w(\mu_j))}{\sigma_2(w(\mu_j))}
    -\frac{\sigma_1(w(\mu_j))}{\sigma_2(w(\mu_j))}h_{j2}.
    \]
    These are constants, not parameters estimated from sampled mode data.
    We compare them with the negatives of regularized second-kind Abel
    primitives, evaluated by a different generapy call. That checks the sigma
    expression, sign and lift convention; both evaluations still share period
    infrastructure, so they are not wholly independent software implementations.
    """)
    fixed_w = []
    fixed_sigma2 = []
    h_coefficients = []
    second_kind_errors = []
    for _coordinate, _ordinate in zip(fixed_rho, fixed_ordinate):
        _w = curve.abel_map_kind_1((_coordinate, _ordinate), reduce=False).value
        _jet = sigma_jet(_w, order=2)
        _normalizer = _jet[(0, 1)]
        _h2 = _jet[(0, 2)]/(2*_normalizer)
        _h1 = _jet[(1, 1)]/_normalizer-_jet[(1, 0)]*_h2/_normalizer
        _second = curve.abel_map_kind_2((_coordinate, _ordinate), reduce=False).value
        fixed_w.append(_w)
        fixed_sigma2.append(_normalizer)
        h_coefficients.append((_h1, _h2))
        second_kind_errors.append(max(abs(a+b)/max(1, abs(a), abs(b)) for a, b in zip((_h1, _h2), _second)))
    assert max(second_kind_errors) < ctx.mpf("1e-23")
    mo.ui.table([{"mode": j+1, "h_j1 (numerical evaluation)": str(value[0]),
                  "h_j2 (numerical evaluation)": str(value[1]),
                  "second-kind comparison": str(error), "fitted coefficients": 0}
                 for j, (value, error) in enumerate(zip(h_coefficients, second_kind_errors))])
    return fixed_sigma2, fixed_w, h_coefficients, second_kind_errors


@app.cell
def _(F_numeric, ctx, curve, largest_real_root_interval, mo, quadrature_orders, rho_end, rho_start):
    mo.md(r"""
    ## 6. Choose a longer regular path and invert the scalar quadrature

    The largest real branch point is about $1.9255$. Our path stays at
    $\rho\ge4.5$, on the positive sheet, and never approaches the fixed mode
    coordinates (all at most $0.5$). Exact real-root isolation checks the branch
    point bound. These choices avoid a period or sheet transition.

    On this branch $\rho'(z)=+\sqrt{F(\rho(z))}$. We obtain the endpoint
    propagation length from the first Abel increment and verify it by scalar
    quadrature. For each prescribed analytic sample value of $z$, a safeguarded
    Newton iteration solves the known equation
    \[
    \int_{4.5}^{\rho(z)}\frac{d\rho}{\sqrt{F(\rho)}}=z-z_0.
    \]
    This is numerical inversion of a specified integral, with derivative
    $1/\sqrt{F}$; it does not fit $F$, $B_j$, $h_{jk}$ or mode parameters.
    Two Gauss--Legendre orders and the forward Abel map check the inversion.
    """)
    assert ctx.mpf(str(largest_real_root_interval[1])) < rho_start < rho_end
    assert F_numeric(rho_start) > 0 and F_numeric(rho_end) > 0
    start_w = curve.abel_map_kind_1((rho_start, ctx.sqrt(F_numeric(rho_start))), reduce=False).value
    end_w = curve.abel_map_kind_1((rho_end, ctx.sqrt(F_numeric(rho_end))), reduce=False).value
    _length = end_w[0]-start_w[0]
    assert abs(_length.imag) < ctx.mpf("1e-25")
    duration = _length.real
    quadratures = {order: ctx.gauss_quadrature(order, "legendre") for order in quadrature_orders}

    def scalar_quadrature(coordinate, order=96):
        _nodes, _weights = quadratures[order]
        _half = (coordinate-rho_start)/2
        _middle = (coordinate+rho_start)/2
        return _half*sum(weight/ctx.sqrt(F_numeric(_middle+_half*node))
                         for node, weight in zip(_nodes, _weights))

    quadrature_endpoint_error = abs(scalar_quadrature(rho_end)-duration)
    quadrature_order_error = abs(scalar_quadrature(rho_end, quadrature_orders[0])-scalar_quadrature(rho_end, quadrature_orders[1]))
    assert max(quadrature_endpoint_error, quadrature_order_error) < ctx.mpf("1e-25")

    def coordinate_at_z(time):
        # Monotonicity supplies a bracket; Newton steps leaving it use bisection.
        if time == 0:
            return rho_start
        if time == duration:
            return rho_end
        _lower, _upper = rho_start, rho_end
        _coordinate = rho_start+(rho_end-rho_start)*time/duration
        for _iteration in range(100):
            _residual = scalar_quadrature(_coordinate)-time
            if abs(_residual) < ctx.mpf("1e-29"):
                return _coordinate
            if _residual > 0:
                _upper = _coordinate
            else:
                _lower = _coordinate
            _proposal = _coordinate-_residual*ctx.sqrt(F_numeric(_coordinate))
            _coordinate = _proposal if _lower < _proposal < _upper else (_lower+_upper)/2
        raise RuntimeError("Scalar inverse did not meet the requested residual")

    return coordinate_at_z, duration, quadrature_endpoint_error, quadrature_order_error, scalar_quadrature, start_w


@app.cell
def _(F_numeric, beta, ctx, fixed_sigma2, fixed_w, mo, numeric_gamma, q_numeric, rho_start, sigma_jet, start_w):
    mo.md(r"""
    ## 7. Specify initial modes on the selected invariant level

    The Hamiltonian only constrains paired powers and the complete weighted
    product. We start from five declared complex seed amplitudes and multiply
    **all five** by a common sixth root to enforce
    $\prod_j u_j(z_0)^{\beta_j}=(q(4.5)+\rho'(z_0))/2$.
    The partner modes then satisfy $v_j(z_0)=I_j(z_0)/u_j(z_0)$.

    This algebraic choice of initial data is not a trajectory fit. It keeps
    initial amplitudes more evenly sized than assigning the whole product
    normalization to a single mode. The seeds and root branch are visible in
    the cell below; they differ from the earlier plot's amplitude choices,
    while the Hamiltonian, scalar curve and invariant offsets are the same.
    """)
    seeds = [ctx.mpc("1.1", ".2"), ctx.mpc(".8", "-.1"), ctx.mpc("1.2", ".1"),
             ctx.mpc("1", ".1"), ctx.mpc("1.1", ".1")]
    initial_product_target = (q_numeric(rho_start)+ctx.sqrt(F_numeric(rho_start)))/2
    initial_scale = ctx.exp(ctx.log(initial_product_target/ctx.fprod(value**weight for value, weight in zip(seeds, beta)))/sum(beta))
    initial_u = [initial_scale*value for value in seeds]
    initial_v = [(offset-weight*rho_start)/value for offset, weight, value in zip(numeric_gamma, beta, initial_u)]
    initial_state = initial_u+initial_v
    _start_sigma2 = sigma_jet(start_w, order=1)[(0, 1)]
    initial_raw_u = []
    initial_raw_v = []
    for _weight, _fixed, _normalizer in zip(beta, fixed_w, fixed_sigma2):
        _denominator = _start_sigma2*_normalizer
        initial_raw_u.append(ctx.sqrt(_weight)*sigma_jet(start_w+_fixed, order=0)[(0, 0)]/_denominator)
        initial_raw_v.append(-ctx.sqrt(_weight)*sigma_jet(start_w-_fixed, order=0)[(0, 0)]/_denominator)
    mo.ui.table([{"mode": j+1, "u_j(z0)": str(a), "v_j(z0)": str(b)} for j, (a, b) in enumerate(zip(initial_u, initial_v))])
    return initial_raw_u, initial_raw_v, initial_state


@app.cell
def _(F_numeric, beta, ctx, curve, fixed_sigma2, fixed_w, gauge_coefficients, h_coefficients, initial_raw_u, initial_raw_v, initial_state, mo, sigma_jet, start_w):
    mo.md(r"""
    ## 8. Evaluate both complete sigma solutions

    The analytic formulas are specified before comparison with RK4. For genus
    two, define the uncorrected factors
    \[
    \bar u_j=\sqrt{\beta_j}\frac{\sigma(w+w(\mu_j))}{\sigma_2(w)\sigma_2(w(\mu_j))},\qquad
    \bar v_j=-\sqrt{\beta_j}\frac{\sigma(w-w(\mu_j))}{\sigma_2(w)\sigma_2(w(\mu_j))}.
    \]
    With $\Delta w=w(z)-w(z_0)$, the full $u_j$ is its initial value times the
    normalized $\bar u_j$ ratio and
    $\exp[(c_{j1}-h_{j1})\Delta w_1+(c_{j2}-h_{j2})\Delta w_2]$.
    The $v_j$ formula uses its own normalized sigma factor and the negative
    exponent. We evaluate **both** formulas: reconstructing $v_j=I_j/u_j$
    would make their paired-power check automatic.

    The exponent contains the exact polynomial correction and the analytic
    second-kind coefficient formula. No least squares, interpolation of mode
    data, numerical derivative fitting, or fitted gauge constants are used.
    """)
    def analytic_at_coordinate(coordinate):
        _ordinate = ctx.sqrt(F_numeric(coordinate))
        _w = curve.abel_map_kind_1((coordinate, _ordinate), reduce=False).value
        _jet = sigma_jet(_w, order=1)
        _dw = _w-start_w
        _u, _v = [], []
        for _j, (_weight, _fixed, _normalizer, _B, _h) in enumerate(zip(beta, fixed_w, fixed_sigma2, gauge_coefficients, h_coefficients)):
            _denominator = _jet[(0, 1)]*_normalizer
            _raw_u = ctx.sqrt(_weight)*sigma_jet(_w+_fixed, order=0)[(0, 0)]/_denominator
            _raw_v = -ctx.sqrt(_weight)*sigma_jet(_w-_fixed, order=0)[(0, 0)]/_denominator
            _exponent = sum((a-b)*increment for a, b, increment in zip(_B, _h, _dw))
            _u.append(initial_state[_j]*_raw_u/initial_raw_u[_j]*ctx.exp(_exponent))
            _v.append(initial_state[_j+5]*_raw_v/initial_raw_v[_j]*ctx.exp(-_exponent))
        return dict(rho=coordinate, ordinate=_ordinate, w=_w, state=_u+_v,
                    inversion=-_jet[(1, 0)]/_jet[(0, 1)], sigma_divisor=_jet[(0, 0)])

    return (analytic_at_coordinate,)


@app.cell
def _(ctx, exact_field, exact_u, exact_v, initial_state, mo, sp):
    # A specialized polynomial implementation makes the fine RK4 grids affordable.
    # No scalar analytic trajectory or division by a mode enters this function.
    _constant = ctx.mpf(107)/400
    _offset = ctx.mpf(1)/10

    def canonical_field(state, constant=_constant, offset=_offset):
        u1, u2, u3, u4, u5, v1, v2, v3, v4, v5 = state
        x1, x2, x3, x4, x5 = u1*v1, u2*v2, u3*v3, u4*v4, u5*v5
        S_value = x1+x2+x3+x4+x5
        common = S_value*S_value/18-S_value/36-constant
        p1, p2 = common-x3/6, common+x3/3
        p3 = common-(x1/2-x2+offset)/3
        vv, uu = v1*v1, u1*u1
        return [p1*u1+2*v1*v2*v3*v4*v5, p2*u2+vv*v3*v4*v5,
                p3*u3+vv*v2*v4*v5, common*u4+vv*v2*v3*v5, common*u5+vv*v2*v3*v4,
                -p1*v1-2*u1*u2*u3*u4*u5, -p2*v2-uu*u3*u4*u5,
                -p3*v3-uu*u2*u4*u5, -common*v4-uu*u2*u3*u5, -common*v5-uu*u2*u3*u4]

    # Check that the handwritten fast field really is dH/dv, -dH/du.
    # The same field function accepts exact constants for this symbolic audit.
    # No numerical coefficient recovery or rational approximation is involved.
    _symbolic_fast = canonical_field(list(exact_u)+list(exact_v),
                                     constant=sp.Rational(107, 400), offset=sp.Rational(1, 10))
    assert all(sp.expand(actual-expected) == 0 for actual, expected in zip(_symbolic_fast, exact_field))
    _zero_probe = list(initial_state)
    _zero_probe[0] = ctx.zero
    assert all(ctx.isfinite(value) for value in canonical_field(_zero_probe))
    mo.md("The fast RK4 field passes the exact Hamiltonian-derivative audit and remains finite at a zero mode.")
    return (canonical_field,)


@app.cell
def _(analytic_at_coordinate, analytic_samples, beta, coordinate_at_z, ctx, duration, numeric_gamma, perf_counter, q_numeric):
    _started = perf_counter()
    analytic_times = [duration*k/(analytic_samples-1) for k in range(analytic_samples)]
    analytic_values = [analytic_at_coordinate(coordinate_at_z(time)) for time in analytic_times]
    analytic_seconds = perf_counter()-_started
    analytic_checks = []
    for _time, _item in zip(analytic_times, analytic_values):
        _state = _item["state"]
        _paired = [_state[j]*_state[j+5] for j in range(5)]
        _U = ctx.fprod(value**weight for value, weight in zip(_state[:5], beta))
        _V = ctx.fprod(value**weight for value, weight in zip(_state[5:], beta))
        _q = q_numeric(_item["rho"])
        _expected = [(offset-weight*_item["rho"]) for offset, weight in zip(numeric_gamma, beta)]
        analytic_checks.append(dict(
            abel_time_error=abs((_item["w"][0]-analytic_values[0]["w"][0]).real-_time),
            sigma_inversion_error=abs(_item["inversion"]-_item["rho"])/max(1, abs(_item["rho"])),
            paired_power_error=max(abs(a-b)/max(1, abs(a), abs(b)) for a, b in zip(_paired, _expected)),
            U_product_error=abs(_U-(_q+_item["ordinate"])/2)/max(1, abs(_U)),
            V_product_error=abs(_V-(_q-_item["ordinate"])/2)/max(1, abs(_V)),
        ))
    assert max(value for row in analytic_checks for value in row.values()) < ctx.mpf("1e-23")
    return analytic_checks, analytic_seconds, analytic_times, analytic_values


@app.cell
def _(canonical_field, ctx):
    def integrate_rk4(initial, duration, steps, observations):
        """Integrate a uniform fine grid and save only requested integer indices."""
        state = list(initial)
        h = duration/steps
        saved = {0: list(state)}
        for index in range(1, steps+1):
            k1 = canonical_field(state)
            k2 = canonical_field([a+h*b/2 for a, b in zip(state, k1)])
            k3 = canonical_field([a+h*b/2 for a, b in zip(state, k2)])
            k4 = canonical_field([a+h*b for a, b in zip(state, k3)])
            state = [a+h*(b+2*c+d*2+e)/6 for a, b, c, d, e in zip(state, k1, k2, k3, k4)]
            if index in observations:
                saved[index] = list(state)
            if index % 32768 == 0:
                print(f"RK4 {steps:,} steps: {index:,} completed", flush=True)
        assert all(ctx.isfinite(value) for value in state)
        return saved

    return (integrate_rk4,)


@app.cell
def _(analytic_samples, analytic_values, beta, ctx, duration, initial_state, integrate_rk4, mo, numeric_gamma, perf_counter, plot_samples, q_numeric, rk4_steps):
    mo.md(r"""
    ## 9. Independently propagate the full canonical equations

    The RK4 routine above performs four polynomial-field evaluations per step.
    The first four weighted invariants and the Hamiltonian are checked on the
    saved nodes. No sigma, Abel, or scalar inverse calculation is used inside
    the integration loop.

    Observation indices are exact integer multiples of the RK4 step: there is
    **no interpolation** when comparing with the 33 analytic values. The finest
    grid also stores 513 states for the dashed curves. RK4 lines are thin and
    black; analytic sigma values are black markers.
    """)
    convergence_rows = []
    rk4_runs = {}
    for _steps in rk4_steps:
        assert _steps % (analytic_samples-1) == 0 and _steps % (plot_samples-1) == 0
        _comparison_indices = [_steps*k//(analytic_samples-1) for k in range(analytic_samples)]
        _plot_indices = [_steps*k//(plot_samples-1) for k in range(plot_samples)] if _steps == rk4_steps[-1] else []
        _started = perf_counter()
        _saved = integrate_rk4(initial_state, duration, _steps, set(_comparison_indices+_plot_indices))
        _seconds = perf_counter()-_started
        _errors, _H_errors, _C_errors = [], [], []
        for _index, _exact in zip(_comparison_indices, analytic_values):
            _state = _saved[_index]
            _errors.append(max(abs(a-b)/max(1, abs(b)) for a, b in zip(_state, _exact["state"])))
            _I = [_state[j]*_state[j+5] for j in range(5)]
            _rho = -sum(_I)/sum(beta)
            _Q = q_numeric(_rho)+(_I[0]/2-_I[1]+ctx.mpf(1)/10)*_I[2]/3
            _U = ctx.fprod(value**weight for value, weight in zip(_state[:5], beta))
            _V = ctx.fprod(value**weight for value, weight in zip(_state[5:], beta))
            _H_errors.append(abs(_U+_V-_Q)/max(1, abs(_U), abs(_V), abs(_Q)))
            _C_errors.append(max(abs((_I[j]-beta[j]*_I[4])-(numeric_gamma[j]-beta[j]*numeric_gamma[4]))
                                 /max(1, abs(numeric_gamma[j]-beta[j]*numeric_gamma[4])) for j in range(4)))
        convergence_rows.append(dict(steps=_steps, step=duration/_steps, seconds=_seconds,
                                     max_mode_error=max(_errors), max_scaled_H_drift=max(_H_errors),
                                     max_scaled_C_drift=max(_C_errors), node_errors=_errors))
        rk4_runs[_steps] = _saved
        print(f"RK4 {_steps:,}: max scaled mode error {ctx.nstr(max(_errors), 10)}", flush=True)
    convergence_ratios = [convergence_rows[k]["max_mode_error"]/convergence_rows[k+1]["max_mode_error"] for k in range(2)]
    assert all(12 < ratio < 20 for ratio in convergence_ratios)
    assert convergence_rows[-1]["max_mode_error"] < ctx.mpf("1e-20")
    finest_states = [rk4_runs[rk4_steps[-1]][rk4_steps[-1]*k//(plot_samples-1)] for k in range(plot_samples)]
    plot_times = [duration*k/(plot_samples-1) for k in range(plot_samples)]
    mo.ui.table([{"RK4 steps": row["steps"], "step size": ctx.nstr(row["step"], 10),
                  "max scaled mode error (all 33 nodes)": ctx.nstr(row["max_mode_error"], 12),
                  "scaled Hamiltonian drift": ctx.nstr(row["max_scaled_H_drift"], 12),
                  "scaled weighted-invariant drift": ctx.nstr(row["max_scaled_C_drift"], 12),
                  "seconds": round(row["seconds"], 2)} for row in convergence_rows])
    return convergence_ratios, convergence_rows, finest_states, plot_times


@app.cell
def _(analytic_times, analytic_values, convergence_rows, ctx, duration, example_dir, finest_states, mo, plot_times, plt, rk4_steps):
    _plot_z = [float(value) for value in plot_times]
    _marker_z = [float(value) for value in analytic_times]
    mode_figures = {}
    for _tag, _offset in (("u", 0), ("v", 5)):
        _figure, _axes = plt.subplots(2, 3, figsize=(11, 7.8), constrained_layout=True)
        for _j, _axis in enumerate(_axes.flat):
            if _j >= 5:
                _axis.set_visible(False)
                continue
            _axis.plot(_plot_z, [float(state[_j+_offset].real) for state in finest_states],
                       color="black", linestyle="--", linewidth=.8,
                       label=f"RK4 ({rk4_steps[-1]:,} steps)")
            _axis.plot(_marker_z, [float(item["state"][_j+_offset].real) for item in analytic_values],
                       color="black", marker="o", markersize=3.2, linestyle="None",
                       label=f"analytic sigma ({len(analytic_values)} samples)")
            _axis.set_xlabel(r"$z-z_0$", fontsize=13)
            _axis.tick_params(labelsize=12)
            _axis.legend([], [], title=rf"$\mathrm{{Re}}\,{_tag}_{{{_j+1}}}(z)$", title_fontsize=15, loc="best", framealpha=.9)
            _axis.grid(alpha=.15)
        _handles, _labels = _axes.flat[0].get_legend_handles_labels()
        _legend_axis = _axes.flat[-1]
        _legend_axis.set_visible(True)
        _legend_axis.set_axis_off()
        _legend_axis.legend(_handles, _labels, loc="center", fontsize=13, frameon=False)
        _figure.savefig(example_dir/f"{_tag}_trajectories.png", dpi=180)
        mode_figures[_tag] = _figure
    mo.vstack([mode_figures["u"], mode_figures["v"]])
    return (mode_figures,)


@app.cell
def _(analytic_samples, analytic_times, analytic_values, convergence_rows, ctx, example_dir, mo, plt):
    error_figure, _axis = plt.subplots(figsize=(9, 4), constrained_layout=True)
    for _row, _style in zip(convergence_rows, (":", "--", "-")):
        _axis.semilogy([float(time) for time in analytic_times[1:]],
                       [float(value) for value in _row["node_errors"][1:]],
                       color="black", linestyle=_style, linewidth=.9,
                       label=f"RK4 {_row['steps']:,} steps")
    _axis.set(xlabel=r"$z-z_0$", ylabel="max scaled mode error at analytic sample",
              title="RK4 refinement against the same sparse sigma solution")
    _axis.grid(alpha=.2)
    _axis.legend()
    error_figure.savefig(example_dir/"rk4_refinement.png", dpi=180)
    mo.vstack([error_figure, mo.md(r"""
    ## 10. Interpret the comparison

    The endpoint errors are only one part of the comparison: the table and error
    plot use the maximum over all ten $u,v$ components at all analytic sample
    values. Halving the step should reduce the error by approximately $16$ for
    fourth-order RK4 until evaluation or rounding errors dominate.

    The Hamiltonian error is scaled by the largest Hamiltonian term, since
    $H=0$ is a cancellation between much larger values. Weighted-invariant
    errors use their conserved offsets. Those scales are stated explicitly;
    they are not the same as an unscaled absolute Hamiltonian residual.

    More digits in the arithmetic alone cannot remove RK4 discretization error.
    The analytic samples are unchanged between the three RK4 grids. The plots
    show numerical RK4 states joined by line segments, with independently
    evaluated sigma markers. No smoothing or interpolation is applied to the
    analytic amplitudes to manufacture agreement.
    """)])
    return (error_figure,)


@app.cell
def _(analytic_checks, analytic_samples, analytic_seconds, analytic_times, analytic_values, convergence_ratios, convergence_rows, ctx, dps, duration, exact_checks, example_dir, json, metadata, mo, period_seconds, plot_samples, quadrature_endpoint_error, quadrature_order_error, rho_end, rho_start, second_kind_errors):
    # Measure visible nonlinearity against the chord through the two endpoints.
    # This is a geometry diagnostic, not a least-squares line or a fitted curve.
    curvature_rows = []
    for _j in range(5):
        _first = analytic_values[0]["state"][_j].real
        _last = analytic_values[-1]["state"][_j].real
        _span = abs(_last-_first)
        _deviation = max(abs(item["state"][_j].real-(_first+(_last-_first)*time/duration))
                         for time, item in zip(analytic_times, analytic_values))
        curvature_rows.append(dict(mode=_j+1, endpoint_relative_change=abs((_last-_first)/_first),
                                   max_chord_deviation_fraction=_deviation/_span))

    def _serial(value):
        if isinstance(value, dict):
            return {key: _serial(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [_serial(item) for item in value]
        if isinstance(value, (bool, int, float, str)) or value is None:
            return value
        return ctx.nstr(value, 16)

    report = dict(passed=True, dps=dps, rho_start=rho_start, rho_end=rho_end, sheet=1,
                  z_duration=duration, analytic_samples=analytic_samples, plot_samples=plot_samples,
                  exact_checks=exact_checks, fitting="None: model, gauges and second-kind formulas are prescribed.",
                  initial_data="All five seed amplitudes receive one common sixth-root scale.",
                  second_kind_max_error=max(second_kind_errors),
                  quadrature_endpoint_error=quadrature_endpoint_error,
                  quadrature_order_error=quadrature_order_error,
                  analytic_maxima={key: max(row[key] for row in analytic_checks) for key in analytic_checks[0]},
                  convergence=convergence_rows, step_halving_ratios=convergence_ratios,
                  curvature=curvature_rows, period_seconds=period_seconds, analytic_seconds=analytic_seconds,
                  versions={package: metadata.version(package) for package in ("generapy", "mpmath", "sympy", "marimo", "matplotlib")},
                  scope="Genus-two example only; positive-sheet regular path, no turning-point or global period-continuation test.")
    (example_dir/"results.json").write_text(json.dumps(_serial(report), indent=2)+"\n")
    mo.vstack([
        mo.md(f"""
        ## 11. Summary and saved artifacts

        Propagation length: **{ctx.nstr(duration, 12)}**; arithmetic: **{dps} dps**.
        Finest-grid maximum scaled mode error: **{ctx.nstr(convergence_rows[-1]['max_mode_error'], 10)}**.
        Step-halving reduction factors: **{', '.join(ctx.nstr(value, 8) for value in convergence_ratios)}**.

        Computing the {analytic_samples} analytic samples took **{analytic_seconds:.2f} s**;
        the individual RK4 grid timings appear above. No analytic evaluation runs
        on the {convergence_rows[-1]['steps']:,}-step integration grid.

        `results.json`, `u_trajectories.png`, `v_trajectories.png`, and
        `rk4_refinement.png` are saved beside this notebook. The table below
        quantifies departure from the straight endpoint chord, as a fraction
        of each real $u_j$'s total endpoint variation. It is a diagnostic;
        no chord is fitted to the data or used to calculate a mode.
        """),
        mo.ui.table([{"mode": row["mode"], "endpoint change (%)": ctx.nstr(100*row["endpoint_relative_change"], 8),
                      "max departure from endpoint chord (%)": ctx.nstr(100*row["max_chord_deviation_fraction"], 8)}
                     for row in curvature_rows])
    ])
    return (report,)


@app.cell
def _(analytic_times, analytic_values, beta, example_dir, finest_states, mo, plot_times, rk4_steps):
    from hyperelliptic_scalar_first.paper_plots import save_paper_plots
    paper_figures = save_paper_plots(example_dir, plot_times, finest_states,
                                    analytic_times, analytic_values, beta, rk4_steps)
    mo.vstack(paper_figures)
    return (paper_figures,)


if __name__ == "__main__":
    app.run()
