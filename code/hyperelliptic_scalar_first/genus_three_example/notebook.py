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
    # Genus three: an oscillatory branch of the second worked example

    We retain the paper's seven-pair Hamiltonian, double weight and invariant
    offsets. The change is **initial data and real branch**, not the model.
    We follow three oscillations of the scalar coordinate and paired powers.
    Individual modes can acquire a multiplier each cycle; periodic powers do
    not imply that every amplitude is periodic. This is a complex canonical
    system: we do not impose $v_j=\overline{u_j}$ or interpret signed powers as
    positive optical intensities.

    | Stage | What is analytic or exact? | What is computed numerically? |
    |---|---|---|
    | Hamiltonian | Rational coefficients, canonical field, curve, corrections. | Nothing is fitted. SymPy performs exact algebra. |
    | Geometry | Prescribed sigma identities and second-kind coefficient formulas. | Generapy evaluates periods, characteristics, forward Abel maps and sigma jets. |
    | Turning points | A smooth angular parametrization of the real oval. | Root finding and quadrature invert a specified scalar integral. |
    | Validation | Full polynomial equations specified independently. | Fine-grid RK4 compared with sparse sigma samples, without interpolation or fitting. |

    Unlike the exterior genus-two path, this calculation crosses turning
    points and modal zeros. We integrate the Abel lift along the actual path,
    rather than resetting it with a principal Abel-map call at every sample.
    The first example and the LaTeX paper are unchanged by this experiment.
    """)
    return


@app.cell
def _(mp, mo):
    dps = 30
    ctx = mp.clone()
    ctx.dps = dps
    cycles = 3
    analytic_samples = 97
    plot_samples = 1537
    rk4_steps = (24576, 49152, 98304)
    quadrature_orders = (64, 96)
    # Each grid is divisible by 96 observation intervals and 1536 plot intervals.
    mo.md(f"""
    ## 1. Precision and independent observation grids

    Arithmetic uses **{dps} decimal digits**. We sample sigma at
    **{analytic_samples}** uniform propagation values over **{cycles} cycles**.
    RK4 uses **{', '.join(f'{n:,}' for n in rk4_steps)} steps**;
    the finest run stores **{plot_samples}** plotting states. All observation
    times are exact RK4 grid nodes. Arithmetic precision does not determine
    finite-step truncation error: the measured refinement results do.
    """)
    return analytic_samples, ctx, cycles, dps, plot_samples, quadrature_orders, rk4_steps


@app.cell
def _(mo, sp):
    rho = sp.Symbol("rho")
    x = sp.symbols("x1:8")
    beta = (2, 1, 1, 1, 1, 1, 1)
    gamma = tuple(map(sp.Rational, ("-3/5", "-1/4", "-1/5", "-1/10", "1/5", "2/5", "11/20")))
    S = sum(x)
    Q = 4*(-S/8+sp.Rational(37, 40))*(-S/8+sp.Rational(3, 20))*(-S/8-sp.Rational(3, 10))*(-S/8-sp.Rational(13, 20))
    Q += (x[0]/2-x[1]+sp.Rational(1, 20))*x[2]**2/7
    powers = tuple(offset-weight*rho for weight, offset in zip(beta, gamma))
    line = dict(zip(x, powers))
    q = sp.expand(Q.subs(line, simultaneous=True))
    r = sp.expand(sp.prod(power**weight for power, weight in zip(powers, beta)))
    F = sp.expand(q*q-4*r)
    Q_partials = tuple(sp.diff(Q, variable) for variable in x)
    fixed_coordinates = tuple(offset/weight for offset, weight in zip(gamma, beta))
    fixed_ordinates = tuple(q.subs(rho, value) for value in fixed_coordinates)
    B = tuple(sp.cancel((q-ordinate)/(2*(coordinate-rho))-partial.subs(line, simultaneous=True))
              for coordinate, ordinate, partial in zip(fixed_coordinates, fixed_ordinates, Q_partials))
    exact_u = sp.symbols("u1:8")
    exact_v = sp.symbols("v1:8")
    exact_H = sp.prod(value**weight for value, weight in zip(exact_u, beta))+sp.prod(value**weight for value, weight in zip(exact_v, beta))
    exact_H -= Q.subs(dict(zip(x, [a*b for a, b in zip(exact_u, exact_v)])), simultaneous=True)
    exact_field = tuple(sp.diff(exact_H, value) for value in exact_v)+tuple(-sp.diff(exact_H, value) for value in exact_u)
    root_intervals = sp.Poly(F, rho).intervals(eps=sp.Rational(1, 10**35))
    exact_checks = dict(centred_offsets=sum(gamma) == 0,
                        degree_seven=sp.degree(F, rho) == 7,
                        leading_coefficient_four=sp.LC(sp.Poly(F, rho)) == 4,
                        squarefree=sp.degree(sp.gcd(F, sp.diff(F, rho)), rho) == 0,
                        quadratic_corrections=all(sp.degree(value, rho) <= 2 for value in B),
                        seven_simple_real_roots=len(root_intervals) == 7 and all(m == 1 for _, m in root_intervals),
                        fixed_points_nonbranch=all(value != 0 for value in fixed_ordinates))
    assert all(exact_checks.values())
    mo.vstack([mo.md(r"""
    ## 2. Derive the exact model before choosing a trajectory

    Set $S(x)=\sum x_j$. The canonical system is $u_j'=\partial H/\partial v_j$,
    $v_j'=-\partial H/\partial u_j$, with
    $H=u_1^2u_2\cdots u_7+v_1^2v_2\cdots v_7-Q$.
    The displayed $Q$ is the paper's polynomial in independent powers;
    we differentiate before substituting $x_j=\gamma_j-\beta_j\rho$.
    The off-line term vanishes on that line but changes individual phase terms.
    """), mo.md(r"\[Q="+sp.latex(Q)+r"\]\[F(\rho)=q(\rho)^2-4r(\rho)="+sp.latex(F)+r".\]"),
        mo.ui.table([{"mode": j+1, "phase = -Q_j": str(-partial), "B_j(rho)": str(value)}
                     for j, (partial, value) in enumerate(zip(Q_partials, B))])])
    return B, F, Q, beta, exact_checks, exact_field, exact_u, exact_v, fixed_coordinates, fixed_ordinates, gamma, q, rho, root_intervals


@app.cell
def _(B, F, ctx, fixed_coordinates, fixed_ordinates, gamma, q, rho, root_intervals, sp):
    def rational_number(value):
        numerator, denominator = sp.fraction(value)
        return ctx.mpf(str(numerator))/ctx.mpf(str(denominator))
    F_coefficients = tuple(rational_number(F.coeff(rho, k)) for k in range(8))
    q_coefficients = tuple(rational_number(q.coeff(rho, k)) for k in range(5))
    gauge_coefficients = tuple(tuple(rational_number(value.coeff(rho, k)) for k in range(3)) for value in B)
    numeric_gamma = tuple(rational_number(value) for value in gamma)
    fixed_rho = tuple(rational_number(value) for value in fixed_coordinates)
    fixed_ordinate = tuple(rational_number(value) for value in fixed_ordinates)
    def F_numeric(coordinate):
        return ctx.polyval(list(reversed(F_coefficients)), coordinate)
    def q_numeric(coordinate):
        return ctx.polyval(list(reversed(q_coefficients)), coordinate)
    roots = tuple(ctx.findroot(F_numeric, (rational_number(interval[0])+rational_number(interval[1]))/2)
                  for interval, _ in root_intervals)
    assert max(abs(F_numeric(value)) for value in roots) < ctx.mpf("1e-25")
    return F_coefficients, F_numeric, fixed_ordinate, fixed_rho, gauge_coefficients, numeric_gamma, q_numeric, rational_number, roots


@app.cell
def _(ctx, mo, roots):
    left, right = roots[4:6]
    centre = (left+right)/2
    radius = (right-left)/2
    theta_start = ctx.pi/2
    # Factor out the two turning roots to remove endpoint singularities.
    def G_numeric(coordinate):
        return -4*ctx.fprod(coordinate-value for j, value in enumerate(roots) if j not in (4, 5))
    def rho_at_theta(theta):
        return centre-radius*ctx.cos(theta)
    def ordinate_at_theta(theta):
        return radius*ctx.sin(theta)*ctx.sqrt(G_numeric(rho_at_theta(theta)))
    assert all(G_numeric(left+(right-left)*k/100) > 0 for k in range(101))
    mo.md(r"""
    ## 3. Select a bounded real oval and regularize turning points

    All seven branch points are real. We select the **fifth and sixth** roots,
    $L\approx0.3505953511$, $R\approx0.5739427200$, where $F\ge0$.
    Write $F(\rho)=(\rho-L)(R-\rho)G(\rho)$ with $G>0$ here, and put
    \[
    \rho=c-d\cos\theta,\quad c=(L+R)/2,\quad d=(R-L)/2,
    \quad \rho'=d\sin\theta\sqrt{G(\rho)},\quad
    \frac{dz}{d\theta}=\frac1{\sqrt{G(\rho)}}.
    \]
    Increasing $\theta$ takes us smoothly through both turning points; the
    ordinate changes sign automatically. The scalar coordinate repeats after
    $2\pi$ in $\theta$, while its speed is not constant. We start at
    $\theta_0=\pi/2$, at the midpoint on the positive sheet.

    This parametrization removes the square-root endpoint singularity from
    scalar and Abel quadrature. It does not add a fitted oscillation or replace
    the curve by a sinusoidal approximation. The powers of modes 6 and 7 cross
    zero at $\rho=0.4$ and $0.55$; full polynomial equations and sigma factors
    can pass through these points even though logarithmic mode rates cannot.
    """)
    return G_numeric, centre, left, ordinate_at_theta, radius, rho_at_theta, right, theta_start


@app.cell
def _(Curve, F_coefficients, ctx, kleinian_sigma_jet, mo, perf_counter):
    mo.md(r"""
    ## 4. Evaluate the curve geometry with generapy

    Generapy receives the curve already derived above. It supplies first- and
    second-kind periods, the Riemann characteristic, forward Abel maps and
    numerical sigma jets. We use hyperelliptic normalization, whose genus-three
    leading polynomial is $w_1w_3-w_2^2-w_2w_3^3/3+w_3^6/45$.
    We implement the real-oval parametrization, scalar inverse, path continuation
    and RK4 below. No generapy call estimates the Hamiltonian or fits mode data.
    """)
    _started = perf_counter()
    curve = Curve({(0, 2): 1} | {(k, 0): -value for k, value in enumerate(F_coefficients) if value}, ctx=ctx)
    first_periods = curve.periods_kind_1()
    second_periods = curve.periods_kind_2()
    sigma_options = dict(omega=first_periods.omega, tau=first_periods.tau,
                         kappa=second_periods.kappa, characteristic=curve.riemann_constant().characteristic,
                         normalization="hyperelliptic", ctx=ctx)
    def sigma_jet(vector, order=2):
        return kleinian_sigma_jet(vector, order=order, **sigma_options)
    period_seconds = perf_counter()-_started
    return curve, period_seconds, sigma_jet


@app.cell
def _(ctx, curve, fixed_ordinate, fixed_rho, mo, sigma_jet):
    mo.md(r"""
    ## 5. Compute prescribed exponential coefficients; fit nothing

    Each fixed point is $(\rho(\mu_j),\rho'(\mu_j))=(\gamma_j/\beta_j,q(\gamma_j/\beta_j))$.
    On its unreduced Abel lift we evaluate
    \[
    h_{j2}=\frac{\sigma_{22}}{2\sigma_2},\qquad
    h_{j1}=\frac{\sigma_{12}}{\sigma_2}-\frac{\sigma_1}{\sigma_2}h_{j2},\qquad
    h_{j3}=\frac{\sigma_{23}}{\sigma_2}.
    \]
    These constants are computed from fixed geometric data, not sampled
    trajectories. We compare them with the negatives of regularized
    second-kind primitives from generapy. That is a convention and identity
    check; both evaluations share the same period infrastructure.
    """)
    fixed_w, fixed_sigma2, h_coefficients, second_kind_errors = [], [], [], []
    for _coordinate, _ordinate in zip(fixed_rho, fixed_ordinate):
        _w = curve.abel_map_kind_1((_coordinate, _ordinate), reduce=False).value
        _jet = sigma_jet(_w, 2)
        _den = _jet[(0, 1, 0)]
        _h2 = _jet[(0, 2, 0)]/(2*_den)
        _h = (_jet[(1, 1, 0)]/_den-_jet[(1, 0, 0)]*_h2/_den, _h2, _jet[(0, 1, 1)]/_den)
        _second = curve.abel_map_kind_2((_coordinate, _ordinate), reduce=False).value
        fixed_w.append(_w)
        fixed_sigma2.append(_den)
        h_coefficients.append(_h)
        second_kind_errors.append(max(abs(a+b)/max(1, abs(a), abs(b)) for a, b in zip(_h, _second)))
    assert max(second_kind_errors) < ctx.mpf("1e-23")
    mo.ui.table([{"mode": j+1, "h_j1": str(value[0]), "h_j2": str(value[1]), "h_j3": str(value[2]),
                  "second-kind error": str(error)} for j, (value, error) in enumerate(zip(h_coefficients, second_kind_errors))])
    return fixed_sigma2, fixed_w, h_coefficients, second_kind_errors


@app.cell
def _(G_numeric, ctx, curve, cycles, mo, ordinate_at_theta, quadrature_orders, rho_at_theta, theta_start):
    mo.md(r"""
    ## 6. Continue the Abel lift around the oval

    Initialize $w_0$ with a single forward Abel-map call. Thereafter use
    \[
    w_k(\theta)-w_k(\theta_0)=
      \int_{\theta_0}^{\theta}\frac{\rho(t)^{k-1}}{\sqrt{G(\rho(t))}}\,dt.
    \]
    The first increment is $z-z_0$. We never discard accumulated periods or
    choose a fresh principal lift after a turning point. Quadrature is split
    into quarter-turn panels and checked at two orders. To sample uniform $z$,
    safeguarded Newton inversion finds $\theta$ from this prescribed integral.
    Root finding is coordinate inversion, not fitting.
    """)
    start_w = curve.abel_map_kind_1((rho_at_theta(theta_start), ordinate_at_theta(theta_start)), reduce=False).value
    quadratures = {order: ctx.gauss_quadrature(order, "legendre") for order in quadrature_orders}
    panel_width = ctx.pi/2
    panels = 4*cycles
    def panel_integral(a, b, order):
        if a == b:
            return [ctx.zero]*3
        nodes, weights = quadratures[order]
        half, middle = (b-a)/2, (b+a)/2
        total = [ctx.zero]*3
        for node, weight in zip(nodes, weights):
            coordinate = rho_at_theta(middle+half*node)
            factor = half*weight/ctx.sqrt(G_numeric(coordinate))
            for k in range(3):
                total[k] += factor*coordinate**k
        return total
    cumulative = {}
    for _order in quadrature_orders:
        _rows = [[ctx.zero]*3]
        for _index in range(panels):
            _increment = panel_integral(theta_start+_index*panel_width, theta_start+(_index+1)*panel_width, _order)
            _rows.append([a+b for a, b in zip(_rows[-1], _increment)])
        cumulative[_order] = _rows
    def abel_increment(theta, order=quadrature_orders[-1]):
        index = min(panels, max(0, int(ctx.floor((theta-theta_start)/panel_width))))
        increment = panel_integral(theta_start+index*panel_width, theta, order)
        return ctx.matrix([a+b for a, b in zip(cumulative[order][index], increment)])
    duration = cumulative[quadrature_orders[-1]][-1][0]
    scalar_period = cumulative[quadrature_orders[-1]][4][0]
    period_increment = ctx.matrix(cumulative[quadrature_orders[-1]][4])
    quadrature_order_error = max(abs(a-b) for row_a, row_b in zip(cumulative[quadrature_orders[0]], cumulative[quadrature_orders[1]]) for a, b in zip(row_a, row_b))
    assert quadrature_order_error < ctx.mpf("1e-25")
    def theta_at_z(time):
        if time == 0:
            return theta_start
        if time == duration:
            return theta_start+panels*panel_width
        times = [row[0] for row in cumulative[quadrature_orders[-1]]]
        index = next(j for j in range(panels) if time <= times[j+1])
        lower, upper = theta_start+index*panel_width, theta_start+(index+1)*panel_width
        theta = lower+(upper-lower)*(time-times[index])/(times[index+1]-times[index])
        for iteration in range(80):
            residual = abel_increment(theta)[0]-time
            if abs(residual) < ctx.mpf("1e-28"):
                return theta
            if residual > 0:
                upper = theta
            else:
                lower = theta
            proposal = theta-residual*ctx.sqrt(G_numeric(rho_at_theta(theta)))
            theta = proposal if lower < proposal < upper else (lower+upper)/2
        raise RuntimeError("Angular inverse failed its residual check")
    mo.md(f"Scalar period: **{ctx.nstr(scalar_period, 16)}**; total propagation: **{ctx.nstr(duration, 16)}**. Quadrature-order difference: **{ctx.nstr(quadrature_order_error, 5)}**.")
    return abel_increment, duration, panel_width, period_increment, quadrature_order_error, scalar_period, start_w, theta_at_z


@app.cell
def _(beta, ctx, mo, numeric_gamma, ordinate_at_theta, q_numeric, rho_at_theta, theta_start):
    mo.md(r"""
    ## 7. Declare initial amplitudes symmetrically

    Choose seven positive real seeds and multiply all by one common eighth
    root so that $U(z_0)=[q(\rho_0)+\rho'_0]/2$. Then set
    $v_j(z_0)=(\gamma_j-\beta_j\rho_0)/u_j(z_0)$.
    This specifies an initial state on $H=0$ with the same invariant offsets.
    It is algebraic initial-data selection, not a fit or an imposed optical
    reality condition. In particular some initial paired powers are negative.
    """)
    initial_seeds = tuple(ctx.mpf(value) for value in ("1.1", ".8", "1.2", "1", "1.1", ".9", "1.05"))
    rho_start = rho_at_theta(theta_start)
    initial_target = (q_numeric(rho_start)+ordinate_at_theta(theta_start))/2
    assert initial_target > 0
    initial_scale = ctx.root(initial_target/ctx.fprod(value**weight for value, weight in zip(initial_seeds, beta)), sum(beta))
    initial_u = [initial_scale*value for value in initial_seeds]
    initial_v = [(offset-weight*rho_start)/value for offset, weight, value in zip(numeric_gamma, beta, initial_u)]
    initial_state = initial_u+initial_v
    mo.ui.table([{"mode": j+1, "seed": str(seed), "u_j(z0)": str(a), "v_j(z0)": str(b)}
                 for j, (seed, a, b) in enumerate(zip(initial_seeds, initial_u, initial_v))])
    return initial_seeds, initial_state, rho_start


@app.cell
def _(abel_increment, beta, ctx, fixed_sigma2, fixed_w, gauge_coefficients, h_coefficients, initial_state, mo, ordinate_at_theta, rho_at_theta, sigma_jet, start_w, theta_start):
    mo.md(r"""
    ## 8. Evaluate both mode formulas, including zeros and sheet changes

    For genus three the point factors are
    \[
    \bar u_j=\sqrt{\beta_j}\frac{\sigma_3(w+w(\mu_j))}{\sigma_2(w)\sigma_2(w(\mu_j))},\qquad
    \bar v_j=\sqrt{\beta_j}\frac{\sigma_3(w-w(\mu_j))}{\sigma_2(w)\sigma_2(w(\mu_j))}.
    \]
    Normalize each by its own initial value and multiply by
    $\exp[\pm\sum_{k=1}^3(c_{jk}-h_{jk})\Delta w_k]$.
    We calculate **both** sigma formulas directly. We never obtain one mode
    by dividing its paired power by the other; this would fail at zeros and
    make the paired-power check automatic. The angular coordinate also avoids
    division by $\rho'$ at turning points. The scalar check uses
    $\rho=-\sigma_1/\sigma_2$ on the one-point image.
    """)
    def raw_factors(w):
        jet = sigma_jet(w, 1)
        raw_u, raw_v = [], []
        for weight, fixed, normalizer in zip(beta, fixed_w, fixed_sigma2):
            denominator = jet[(0, 1, 0)]*normalizer
            raw_u.append(ctx.sqrt(weight)*sigma_jet(w+fixed, 1)[(0, 0, 1)]/denominator)
            raw_v.append(ctx.sqrt(weight)*sigma_jet(w-fixed, 1)[(0, 0, 1)]/denominator)
        return raw_u, raw_v, jet
    initial_raw_u, initial_raw_v, _initial_jet = raw_factors(start_w)
    def analytic_at_theta(theta):
        dw = abel_increment(theta)
        w = start_w+dw
        raw_u, raw_v, jet = raw_factors(w)
        u_values, v_values = [], []
        for j, (B_value, h_value) in enumerate(zip(gauge_coefficients, h_coefficients)):
            exponent = sum((a-b)*increment for a, b, increment in zip(B_value, h_value, dw))
            u_values.append(initial_state[j]*raw_u[j]/initial_raw_u[j]*ctx.exp(exponent))
            v_values.append(initial_state[j+7]*raw_v[j]/initial_raw_v[j]*ctx.exp(-exponent))
        return dict(theta=theta, rho=rho_at_theta(theta), ordinate=ordinate_at_theta(theta), w=w,
                    state=u_values+v_values, inversion=-jet[(1, 0, 0)]/jet[(0, 1, 0)],
                    sigma_divisor=jet[(0, 0, 0)], sigma3_divisor=jet[(0, 0, 1)])
    return (analytic_at_theta,)


@app.cell
def _(ctx, exact_field, exact_u, exact_v, mo, sp):
    # Pure polynomial Hamiltonian field. Prefix/suffix products avoid mode division.
    def canonical_field(state, constants=(ctx.mpf(3)/1024, ctx.mpf(11)/128, ctx.mpf(249)/6400, ctx.mpf(1)/20)):
        u, v = state[:7], state[7:]
        powers = [a*b for a, b in zip(u, v)]
        S_value = sum(powers)
        common = -S_value**3/256+constants[0]*S_value**2+constants[1]*S_value+constants[2]
        phases = [common-powers[2]**2/14, common+powers[2]**2/7,
                  common-2*(powers[0]/2-powers[1]+constants[3])*powers[2]/7]+[common]*4
        def product_derivatives(values):
            factors = [values[0]**2]+values[1:]
            prefix, suffix = [1], [1]*8
            for value in factors:
                prefix.append(prefix[-1]*value)
            for j in range(6, -1, -1):
                suffix[j] = suffix[j+1]*factors[j]
            return [2*values[0]*suffix[1]]+[prefix[j]*suffix[j+1] for j in range(1, 7)]
        u_derivatives, v_derivatives = product_derivatives(u), product_derivatives(v)
        return [phase*value+product for phase, value, product in zip(phases, u, v_derivatives)]+[-phase*value-product for phase, value, product in zip(phases, v, u_derivatives)]
    _symbolic = canonical_field(list(exact_u)+list(exact_v), constants=(sp.Rational(3, 1024), sp.Rational(11, 128), sp.Rational(249, 6400), sp.Rational(1, 20)))
    assert all(sp.expand(a-b) == 0 for a, b in zip(_symbolic, exact_field))
    mo.md(r"""
    ## 9. Audit the fast full-system field

    Put $A=-S^3/256+3S^2/1024+11S/128+249/6400$.
    The phase coefficients are
    $\Phi_1=A-x_3^2/14$, $\Phi_2=A+x_3^2/7$,
    $\Phi_3=A-2(x_1/2-x_2+1/20)x_3/7$, and $\Phi_j=A$ for $j\ge4$.
    The coupled equations are
    \[
    u_1'=\Phi_1u_1+2v_1\prod_{k=2}^7v_k,\quad
    v_1'=-\Phi_1v_1-2u_1\prod_{k=2}^7u_k,
    \]
    \[
    u_j'=\Phi_ju_j+v_1^2\prod_{\substack{k=2\\k\ne j}}^7v_k,\quad
    v_j'=-\Phi_jv_j-u_1^2\prod_{\substack{k=2\\k\ne j}}^7u_k\quad(j\ge2).
    \]
    The fast implementation matches direct differentiation of the full
    Hamiltonian **exactly**. Prefix/suffix products evaluate it without dividing
    by a mode. RK4 does not call sigma, the scalar inverse or an Abel map.
    """)
    return (canonical_field,)


@app.cell
def _(abel_increment, analytic_at_theta, analytic_samples, beta, ctx, duration, mo, numeric_gamma, panel_width, perf_counter, q_numeric, rho_at_theta, theta_at_z, theta_start, cycles):
    _started = perf_counter()
    analytic_times = [duration*k/(analytic_samples-1) for k in range(analytic_samples)]
    analytic_values = [analytic_at_theta(theta_at_z(time)) for time in analytic_times]
    # Extra exact-geometry observations cover every turning point and modal zero.
    event_thetas = [theta_start+j*panel_width for j in range(1, 4*cycles+1) if j % 2 == 1]
    for _coordinate in (ctx.mpf(".4"), ctx.mpf(".55")):
        _theta = ctx.acos((rho_at_theta(ctx.pi)-2*_coordinate+rho_at_theta(ctx.zero))/(rho_at_theta(ctx.pi)-rho_at_theta(ctx.zero)))
        for _turn in range(cycles+1):
            for _value in (_theta+2*ctx.pi*_turn, 2*ctx.pi-_theta+2*ctx.pi*_turn):
                if theta_start <= _value <= theta_start+2*ctx.pi*cycles:
                    event_thetas.append(_value)
    event_values = [analytic_at_theta(theta) for theta in sorted(event_thetas)]
    analytic_seconds = perf_counter()-_started
    def check_item(item):
        state = item["state"]
        paired = [state[j]*state[j+7] for j in range(7)]
        expected = [offset-weight*item["rho"] for offset, weight in zip(numeric_gamma, beta)]
        U = ctx.fprod(value**weight for value, weight in zip(state[:7], beta))
        V = ctx.fprod(value**weight for value, weight in zip(state[7:], beta))
        return dict(sigma_inversion_error=abs(item["inversion"]-item["rho"]),
                    paired_power_error=max(abs(a-b)/max(1, abs(a), abs(b)) for a, b in zip(paired, expected)),
                    U_product_error=abs(U-(q_numeric(item["rho"])+item["ordinate"])/2)/max(1, abs(U)),
                    V_product_error=abs(V-(q_numeric(item["rho"])-item["ordinate"])/2)/max(1, abs(V)),
                    real_mode_error=max(abs(value.imag)/max(1, abs(value)) for value in state))
    analytic_checks = [check_item(item) for item in analytic_values+event_values]
    abel_time_error = max(abs(abel_increment(item["theta"])[0]-time) for time, item in zip(analytic_times, analytic_values))
    assert max(value for row in analytic_checks for value in row.values()) < ctx.mpf("1e-22")
    assert abel_time_error < ctx.mpf("1e-26")
    print(f"Sigma samples and {len(event_values)} turning/zero events checked in {analytic_seconds:.2f}s", flush=True)
    mo.md(f"**{len(analytic_values)}** uniform sigma samples and **{len(event_values)}** additional turning/zero observations pass the analytic checks.")
    return abel_time_error, analytic_checks, analytic_seconds, analytic_times, analytic_values, event_values


@app.cell
def _(canonical_field, ctx):
    def integrate_rk4(initial, duration, steps, observations):
        state = list(initial)
        h = duration/steps
        saved = {0: list(state)}
        for index in range(1, steps+1):
            k1 = canonical_field(state)
            k2 = canonical_field([a+h*b/2 for a, b in zip(state, k1)])
            k3 = canonical_field([a+h*b/2 for a, b in zip(state, k2)])
            k4 = canonical_field([a+h*b for a, b in zip(state, k3)])
            state = [a+h*(b+2*c+2*d+e)/6 for a, b, c, d, e in zip(state, k1, k2, k3, k4)]
            if index in observations:
                saved[index] = list(state)
            if index % 65536 == 0:
                print(f"RK4 {steps:,}: {index:,} completed", flush=True)
        assert all(ctx.isfinite(value) for value in state)
        return saved
    return (integrate_rk4,)


@app.cell
def _(analytic_samples, analytic_values, beta, ctx, duration, initial_state, integrate_rk4, mo, numeric_gamma, perf_counter, plot_samples, q_numeric, rk4_steps):
    convergence_rows, runs = [], {}
    for _steps in rk4_steps:
        assert _steps % (analytic_samples-1) == 0 and _steps % (plot_samples-1) == 0
        _indices = [_steps*k//(analytic_samples-1) for k in range(analytic_samples)]
        _plot_indices = [_steps*k//(plot_samples-1) for k in range(plot_samples)] if _steps == rk4_steps[-1] else []
        _started = perf_counter()
        _saved = integrate_rk4(initial_state, duration, _steps, set(_indices+_plot_indices))
        _errors, _H_errors, _C_errors = [], [], []
        for _index, _exact in zip(_indices, analytic_values):
            _state = _saved[_index]
            _errors.append(max(abs(a-b)/max(1, abs(b)) for a, b in zip(_state, _exact["state"])))
            _I = [_state[j]*_state[j+7] for j in range(7)]
            _rho = -sum(_I)/sum(beta)
            _Q = q_numeric(_rho)+(_I[0]/2-_I[1]+ctx.mpf(1)/20)*_I[2]**2/7
            _U = ctx.fprod(value**weight for value, weight in zip(_state[:7], beta))
            _V = ctx.fprod(value**weight for value, weight in zip(_state[7:], beta))
            _H_errors.append(abs(_U+_V-_Q)/max(1, abs(_U), abs(_V), abs(_Q)))
            _C_errors.append(max(abs((_I[j]-beta[j]*_I[6])-(numeric_gamma[j]-beta[j]*numeric_gamma[6]))
                                 /max(1, abs(numeric_gamma[j]-beta[j]*numeric_gamma[6])) for j in range(6)))
        convergence_rows.append(dict(steps=_steps, step=duration/_steps, seconds=perf_counter()-_started,
                                     max_mode_error=max(_errors), max_scaled_H_drift=max(_H_errors),
                                     max_scaled_C_drift=max(_C_errors), node_errors=_errors))
        runs[_steps] = _saved
        print(f"RK4 {_steps:,}: mode error {ctx.nstr(max(_errors), 10)}", flush=True)
    convergence_ratios = [convergence_rows[j]["max_mode_error"]/convergence_rows[j+1]["max_mode_error"] for j in range(2)]
    assert all(12 < value < 20 for value in convergence_ratios)
    assert convergence_rows[-1]["max_mode_error"] < ctx.mpf("1e-12")
    finest_states = [runs[rk4_steps[-1]][rk4_steps[-1]*k//(plot_samples-1)] for k in range(plot_samples)]
    plot_times = [duration*k/(plot_samples-1) for k in range(plot_samples)]
    mo.vstack([mo.md(r"""
    ## 10. Measure RK4 accuracy over all three cycles

    The scaled mode error is $\max|s^{RK4}-s^\sigma|/\max(1,|s^\sigma|)$
    over all 14 components and 97 common observation nodes. Hamiltonian drift
    uses $\max(1,|U|,|V|,|Q|)$; invariant drift for $C_j=I_j-\beta_jI_7$
    uses $\max(1,|\gamma_j-\beta_j\gamma_7|)$.
    Halving the step should reduce fourth-order truncation error by about 16.
    We require the finest mode discrepancy below $10^{-12}$, not 30-digit RK4
    accuracy over this much longer interval. The achieved errors appear below.
    """), mo.ui.table([{key: ctx.nstr(value, 12) if key not in ("steps", "seconds") else value
                            for key, value in row.items() if key != "node_errors"} for row in convergence_rows])])
    return convergence_ratios, convergence_rows, finest_states, plot_times


@app.cell
def _(analytic_times, analytic_values, beta, ctx, example_dir, finest_states, mo, numeric_gamma, plot_times, plt, rk4_steps):
    _plot_z = [float(value) for value in plot_times]
    _marker_z = [float(value) for value in analytic_times]
    figures = {}
    for _tag, _offset in (("u", 0), ("v", 7), ("powers", None)):
        _figure, _axes = plt.subplots(3, 3, figsize=(11, 10.2), constrained_layout=True)
        for _j, _axis in enumerate(_axes.flat):
            if _j >= 7:
                _axis.set_visible(False)
                continue
            if _offset is None:
                _line = [float((state[_j]*state[_j+7]).real) for state in finest_states]
                _markers = [float((item["state"][_j]*item["state"][_j+7]).real) for item in analytic_values]
                _ylabel = rf"$I_{{{_j+1}}}=u_{{{_j+1}}}v_{{{_j+1}}}$"
            else:
                _line = [float(state[_j+_offset].real) for state in finest_states]
                _markers = [float(item["state"][_j+_offset].real) for item in analytic_values]
                _ylabel = rf"$\mathrm{{Re}}\,{_tag}_{{{_j+1}}}$"
            _axis.plot(_plot_z, _line, color="black", linestyle="--", linewidth=.8, label=f"RK4 ({rk4_steps[-1]:,} steps)")
            _axis.plot(_marker_z, _markers, color="black", marker="o", markersize=2.5, linestyle="None", label="analytic sigma (97 samples)")
            _axis.set_xlabel(r"$z-z_0$", fontsize=13)
            _axis.tick_params(labelsize=12)
            _axis.legend([], [], title=_ylabel, title_fontsize=15, loc="best", framealpha=.9)
            _axis.grid(alpha=.15)
        _handles, _labels = _axes.flat[0].get_legend_handles_labels()
        _legend_axis = _axes.flat[-1]
        _legend_axis.set_visible(True)
        _legend_axis.set_axis_off()
        _legend_axis.legend(_handles, _labels, loc="center", fontsize=13, frameon=False)
        _figure.savefig(example_dir/f"{_tag}_trajectories.png", dpi=180)
        figures[_tag] = _figure
    scalar_figure, _axis = plt.subplots(figsize=(10, 4), constrained_layout=True)
    _axis.plot(_plot_z, [float((-sum(state[j]*state[j+7] for j in range(7))/sum(beta)).real) for state in finest_states],
               color="black", linestyle="--", linewidth=.8, label="rho from full-system RK4")
    _axis.plot(_marker_z, [float(item["inversion"].real) for item in analytic_values],
               color="black", marker="o", markersize=2.5, linestyle="None", label="rho from sigma ratio")
    _axis.set(xlabel=r"$z-z_0$", ylabel=r"$\rho(z)$", title="Three scalar oscillations, with both sheet changes each cycle")
    _axis.grid(alpha=.15)
    _axis.legend()
    scalar_figure.savefig(example_dir/"scalar_oscillations.png", dpi=180)
    mo.vstack([scalar_figure, figures["powers"], figures["u"], figures["v"]])
    return figures, scalar_figure


@app.cell
def _(analytic_times, convergence_rows, example_dir, mo, plt):
    error_figure, _axis = plt.subplots(figsize=(10, 4), constrained_layout=True)
    for _row, _style in zip(convergence_rows, (":", "--", "-")):
        _axis.semilogy([float(value) for value in analytic_times[1:]], [float(value) for value in _row["node_errors"][1:]],
                       color="black", linestyle=_style, linewidth=.9, label=f"RK4 {_row['steps']:,}")
    _axis.set(xlabel=r"$z-z_0$", ylabel="max scaled mode discrepancy", title="Refinement against the same sparse sigma samples")
    _axis.grid(alpha=.2)
    _axis.legend()
    error_figure.savefig(example_dir/"rk4_refinement.png", dpi=180)
    mo.vstack([error_figure, mo.md(r"""
    ## 11. What oscillates, and what does not have to repeat?

    The scalar and signed paired powers are periodic on this oval. Mode 6 and
    mode 7 have power zeros, and individual modes can pass smoothly through
    zero. The $u_j,v_j$ formulas acquire reciprocal cycle multipliers: repeated
    powers do not force repeated amplitudes. We display actual real amplitudes,
    without removing an exponential envelope or normalizing each plotted curve.
    This three-cycle test demonstrates continuation on this particular real
    oval, not arbitrary global continuation on every branch of every curve.
    """)])
    return (error_figure,)


@app.cell
def _(abel_time_error, analytic_checks, analytic_samples, analytic_seconds, analytic_values, convergence_ratios, convergence_rows, ctx, cycles, dps, duration, event_values, exact_checks, example_dir, initial_seeds, initial_state, json, left, metadata, mo, period_increment, period_seconds, plot_samples, quadrature_order_error, right, roots, scalar_period, second_kind_errors):
    cycle_multipliers = [analytic_values[-1]["state"][j]/initial_state[j] for j in range(7)]
    power_repeat_error = max(abs(analytic_values[-1]["state"][j]*analytic_values[-1]["state"][j+7]-initial_state[j]*initial_state[j+7]) for j in range(7))
    assert power_repeat_error < ctx.mpf("1e-22")
    def serial(value):
        if isinstance(value, dict):
            return {key: serial(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [serial(item) for item in value]
        if isinstance(value, (bool, int, float, str)) or value is None:
            return value
        return ctx.nstr(value, 16)
    report = dict(passed=True, dps=dps, cycles=cycles, scalar_period=scalar_period, z_duration=duration,
                  real_roots=roots, oscillation_interval=[left, right], analytic_samples=analytic_samples,
                  plot_samples=plot_samples, turning_and_zero_observations=len(event_values), exact_checks=exact_checks,
                  fitting="None: model, corrections, h coefficients and curve are prescribed.",
                  initial_seeds=initial_seeds, initial_state=initial_state,
                  abel_period_increment=list(period_increment), quadrature_order_error=quadrature_order_error,
                  abel_time_error=abel_time_error, second_kind_max_error=max(second_kind_errors),
                  analytic_maxima={key: max(row[key] for row in analytic_checks) for key in analytic_checks[0]},
                  three_cycle_u_multipliers=cycle_multipliers, power_repeat_error=power_repeat_error,
                  convergence=convergence_rows, step_halving_ratios=convergence_ratios,
                  period_seconds=period_seconds, analytic_seconds=analytic_seconds,
                  versions={package: metadata.version(package) for package in ("generapy", "mpmath", "sympy", "marimo", "matplotlib")},
                  scope="Three cycles on one bounded real oval, including turning points and modal zeros; no optical conjugacy imposed.")
    (example_dir/"results.json").write_text(json.dumps(serial(report), indent=2)+"\n")
    mo.vstack([mo.md(f"""
    ## 12. Saved artifacts and interpretation

    Scalar period: **{ctx.nstr(scalar_period, 12)}**; propagation length:
    **{ctx.nstr(duration, 12)}**. Finest mode discrepancy:
    **{ctx.nstr(convergence_rows[-1]['max_mode_error'], 10)}**.
    Refinement ratios: **{', '.join(ctx.nstr(value, 8) for value in convergence_ratios)}**.
    Sigma samples plus turning/zero checks took **{analytic_seconds:.2f}s**.

    `results.json` records the settings, versions, errors, initial data and
    accumulated Abel increment. The saved figures show scalar oscillations,
    signed powers, actual mode amplitudes and RK4 refinement. The table gives
    total $u_j$ multipliers over three cycles; the partner multipliers are
    reciprocal. No multiplier has been fitted or removed from the plots.
    """), mo.ui.table([{"mode": j+1, "three-cycle u multiplier": str(value)} for j, value in enumerate(cycle_multipliers)])])
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
