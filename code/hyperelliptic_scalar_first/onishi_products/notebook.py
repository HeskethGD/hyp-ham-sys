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
    import marimo as mo
    import sympy as sp
    import importlib.metadata as metadata
    from mpmath import mp
    from generapy import Curve, kleinian_sigma_jet
    from hyperelliptic_scalar_first.model import example, genus_three_example
    return Curve, example, example_dir, genus_three_example, json, kleinian_sigma_jet, metadata, mo, mp, sp


@app.cell
def _(mo):
    mo.md(r"""
    # The complete weighted mode products are Onishi determinants

    The paper has already obtained $U=(q+\rho')/2$ and $V=(q-\rho')/2$.
    This notebook asks how to represent these known curve functions by Onishi
    determinants and connect their sigma expressions with our mode solutions.
    We know the genus $g$, the mode count $N$ and the scalar polynomial $q$;
    the matrix size is still to be chosen.

    We first define the allowed columns and determine the required size. We
    then fit the determinant to the product, examine all point factors in its
    sigma expression, and compare the complete mode formulas. The exact same
    two Hamiltonians are imported from `hyperelliptic_scalar_first/model.py`;
    all determinant and sigma checks below are implemented here. Nothing is
    fitted to trajectory data. The final section states the scope and remaining
    extensions of the argument.
    """)
    return


@app.cell
def _(example, genus_three_example, mp):
    dps = 30
    ctx = mp.clone()
    ctx.dps = dps
    models = (example(), genus_three_example())
    return ctx, dps, models


@app.cell
def _(mo):
    mo.md(r"""
    ## 1. Determine the allowed columns before choosing the matrix size

    An Onishi determinant evaluates a prescribed basis of polynomial functions
    of both curve coordinates in increasing pole order at infinity, one curve
    point per row. The standard basis includes pure coordinate powers and mixed
    coordinate--ordinate monomials. The number of rows and
    columns is $n$; it is not an input mode count. With our monic ordinate
    $\rho'/2$, the column sequence begins
    \[
    1,\rho,\ldots,\rho^g,\rho'/2,\rho^{g+1},\rho\rho'/2,\rho^{g+2},\ldots.
    \]
    This order follows from pole orders $2$ for $\rho$ and $2g+1$ for the
    ordinate. These are **initial segments**, not freely selected monomials.

    To obtain $q\pm\rho'$, we need the ordinate and all polynomial powers
    present in $q$. If $\deg q\le g$, the first suitable list has $g+2$
    columns. If $\deg q=g+1$, it has $g+3$. In our examples this gives size
    five for the cubic genus-two restriction and size six for the quartic
    genus-three restriction, although their Hamiltonians have five and seven
    mode pairs respectively.

    Both lists contain **exactly one ordinate column**. Its moving-row
    coefficient is a cofactor involving only fixed rows, so it is constant
    in the moving coordinate. We check below that it is a nonzero Vandermonde
    minor. This explains why we obtain an ordinate alone, not an unknown
    polynomial multiplying it. The next column, $\rho\rho'/2$, would generally
    introduce a linear ordinate multiplier. A higher-degree $q$ needs a longer
    prefix and cannot be handled by silently omitting that column.

    Having chosen $n$, we can choose $n-1$ fixed rows and append the moving
    row. We next count the zeros of the resulting product function; only then
    will we introduce the part of Onishi's sigma expression that represents
    points beyond those fixed rows.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 2. Count the complete intersection divisor, including weights

    Write $P_j=(\rho(\mu_j),\rho'(\mu_j))$ with
    $\rho(\mu_j)=\gamma_j/\beta_j$ and
    $\rho'(\mu_j)=q(\rho(\mu_j))$. From
    \[
    (q-\rho')(q+\rho')=q^2-F=4\prod_j
          [\gamma_j-\beta_j\rho]^{\beta_j}
    \]
    and $q(\rho(\mu_j))\ne0$, the zero of $q-\rho'$ at $P_j$
    has order **$\beta_j$**. This is a contact multiplicity, not a duplicate
    coordinate pair in the Hamiltonian. Its only pole is at infinity.

    Let $M=\sum_j\beta_j$. The two interpolation families are:

    | Restriction | Pole order / total multiplicity | Main determinant size | Fixed rows | Residual points |
    |---|---:|---:|---:|---:|
    | $\deg q\le g$ | $M=2g+1$ | $n=g+2$ | $g+1$ | $g$ |
    | $\deg q=g+1$ | $M=2g+2$ | $n=g+3$ | $g+2$ | $g$ |

    Both present examples lie in the second family. We choose the main list
    $A=(P_1,\ldots,P_{g+2})$ and the residual list
    $B=(P_1,P_{g+3},\ldots,P_N)$. Each list is internally distinct; their
    shared point accounts for the weight two of mode 1.

    Abel's theorem applied to the meromorphic function $q-\rho'$ gives
    $\sum_j\beta_j w(\mu_j)\equiv0$ modulo the period lattice. For the
    composite split below we choose compatible representatives with exact
    weighted sum zero. Both examples have unit-weight points, so any resulting
    period vector can be absorbed into one such point's lift. This choice is
    distinct from the scalar centring condition $\sum_j\gamma_j=0$.
    """)
    return


@app.cell
def _(models, mo, sp):
    def pole_basis(genus, size, coordinate, ordinate):
        # The monic ordinate is ordinate/2: its square is F/4.
        terms = [(2*k, coordinate**k) for k in range(size)]
        terms += [(2*genus+1+2*k, coordinate**k*ordinate/2) for k in range(size)]
        return tuple(term for _, term in sorted(terms, key=lambda pair: pair[0])[:size])
    exact_reports = []
    exact_matrices = []
    for _data in models:
        _g = int(_data['genus'])
        _rho = _data['rho']
        _ordinate = sp.Symbol('rho_prime')
        _A = list(range(_g+2))
        _B = [0]+list(range(_g+2, len(_data['beta'])))
        _coordinates = _data['a']
        _ordinates = _data['b']
        _basis = pole_basis(_g, _g+3, _rho, _ordinate)
        _Delta = sp.prod(_coordinates[j]-_coordinates[i] for k, i in enumerate(_A) for j in _A[k+1:])
        _rows = [[term.subs({_rho: _coordinates[j], _ordinate: _ordinates[j]}, simultaneous=True) for term in _basis] for j in _A]
        _matrix = sp.Matrix(_rows+[list(_basis)])
        _D = sp.expand(_matrix.det())
        _sheet_rows = [[term.subs({_rho: _coordinates[j], _ordinate: -_ordinates[j]}, simultaneous=True) for term in _basis] for j in _A]
        _sheet_D = sp.expand(sp.Matrix(_sheet_rows+[list(_basis)]).det())
        _interpolation = sp.interpolate([(_coordinates[j], _ordinates[j]) for j in _A], _rho)
        # Reorder the sole ordinate column to the end before taking a Schur complement.
        _W = sp.Matrix([[_coordinates[j]**k for k in range(_g+2)] for j in _A])
        _fixed_vector = sp.Matrix([_ordinates[j] for j in _A])
        _power_vector = sp.Matrix([_rho**k for k in range(_g+2)])
        _schur_coefficients = _W.inv()*_fixed_vector
        _q_from_inverse = sp.expand((_power_vector.T*_schur_coefficients)[0])
        _schur = (_ordinate-_q_from_inverse)/2
        _block = _W.row_join(_fixed_vector/2).col_join(_power_vector.T.row_join(sp.Matrix([[_ordinate/2]])))
        _auxiliary = sp.Matrix([[_coordinates[j]**k for k in range(_g+1)] for j in _B]+[[_rho**k for k in range(_g+1)]])
        _Delta_B = sp.prod(_coordinates[j]-_coordinates[i] for k, i in enumerate(_B) for j in _B[k+1:])
        _checks = dict(
            actual_q_recovered=sp.expand(_interpolation-_data['q']) == 0,
            divisor_norm=sp.expand(_data['q']**2-_data['F']-4*_data['r']) == 0,
            total_weight=sum(_data['beta']) == 2*_g+2,
            main_lists_distinct=len(set(_A)) == len(_A),
            residual_lists_distinct=len(set(_B)) == len(_B),
            multiplicities=[(_A+_B).count(j) for j in range(len(_data['beta']))] == list(_data['beta']),
            main_V_determinant=sp.expand(_D-_Delta*(_data['q']-_ordinate)/2) == 0,
            sheet_involution_U_determinant=sp.expand(-_sheet_D-_Delta*(_data['q']+_ordinate)/2) == 0,
            auxiliary_Vandermonde=sp.expand(_auxiliary.det()-_Delta_B*sp.prod(_rho-_coordinates[j] for j in _B)) == 0,
            coefficient_minor_nonzero=_Delta != 0,
            Vandermonde_inverse_recovers_q=sp.expand(_q_from_inverse-_data['q']) == 0,
            Schur_complement_is_half_difference=sp.expand(_schur-(_ordinate-_data['q'])/2) == 0,
            Schur_block_determinant=sp.expand(_block.det()-_W.det()*_schur) == 0,
            Onishi_Schur_column_sign=sp.expand(_D+_Delta*_schur) == 0,
        )
        for _j, _weight in enumerate(_data['beta']):
            _norm = _data['q']**2-_data['F']
            _checks[f'contact_order_mode_{_j+1}'] = all(sp.diff(_norm, _rho, k).subs(_rho, _coordinates[_j]) == 0 for k in range(_weight)) and sp.diff(_norm, _rho, _weight).subs(_rho, _coordinates[_j]) != 0
        assert all(_checks.values())
        exact_reports.append(dict(genus=_g, checks=_checks, Delta=str(_Delta), Delta_residual=str(_Delta_B), main_indices=[j+1 for j in _A], residual_indices=[j+1 for j in _B], schur_q_latex=sp.latex(_q_from_inverse), schur_coefficients=[str(value) for value in _schur_coefficients]))
        exact_matrices.append((_matrix, _auxiliary))
    # Check the lower-degree family's column orientation in arbitrary low genera.
    lower_family_checks = []
    for _g in range(2, 6):
        _rho, _ordinate = sp.symbols('rho rho_prime')
        _q = 2*_rho**_g+_rho+1
        _nodes = list(range(1, _g+2))
        _basis = pole_basis(_g, _g+2, _rho, _ordinate)
        _rows = [[term.subs({_rho: node, _ordinate: _q.subs(_rho, node)}, simultaneous=True) for term in _basis] for node in _nodes]
        _Delta = sp.prod(_nodes[j]-_nodes[i] for i in range(len(_nodes)) for j in range(i+1, len(_nodes)))
        _ok = sp.expand(sp.Matrix(_rows+[list(_basis)]).det()-_Delta*(_ordinate-_q)/2) == 0
        assert _ok
        lower_family_checks.append(dict(genus=_g, interpolation_columns_passed=_ok))
    mo.vstack([mo.md(r"""
    ## 3. Interpolation makes the products determinants

    Use Onishi's columns in pole order, with the **monic ordinate** $\rho'/2$.
    For our examples they are $1,\rho,\ldots,\rho^g,\rho'/2,\rho^{g+1}$.
    The fixed rows give interpolation conditions for the already prescribed $q$.
    Subtract $q/2$ from the ordinate column. All fixed entries vanish; expanding
    the bottom entry gives
    \[
    D_n(P_1,\ldots,P_{g+2},P(z))
      =\Delta_A\frac{q(\rho(z))-\rho'(z)}2.
    \]
    The minus sign relative to $\rho'-q$ comes from moving the ordinate column
    past $\rho^{g+1}$. Reversing **only the fixed** ordinates gives
    \[
    D_n(P_1^*,\ldots,P_{g+2}^*,P(z))
      =-\Delta_A\frac{q(\rho(z))+\rho'(z)}2.
    \]
    Hence $V=D_n/\Delta_A$, $U=-D_n^*/\Delta_A$ once the product normalization
    is established. These are literal Onishi alternants; the $1/2$ in the
    ordinate column reconciles our leading-four curve with his monic model.
    """)]+[mo.md(r"Genus "+str(row['genus'])+r": \[D=\det"+sp.latex(matrix, symbol_names={sp.Symbol("rho_prime"): r"\rho'(z)", sp.Symbol("rho"): r"\rho(z)"})+r"\]\[\Delta_A="+row['Delta']+r".\]") for row, (matrix, _) in zip(exact_reports, exact_matrices)])
    return exact_reports, lower_family_checks, pole_basis


@app.cell
def _(exact_reports, mo):
    mo.vstack([mo.md(r"""
    ## 4. Recover q by the Vandermonde inverse and identify the Schur complement

    The original columns for our examples are
    $(1,\rho,\ldots,\rho^g,\rho'/2,\rho^{g+1})$. Their upper-left block is
    **not yet Vandermonde**. Swap the final polynomial column with the ordinate
    column. There is one sign change; the reordered determinant is $-D_n$.
    More generally write this sign as $\varepsilon$, equal to $+1$ for the
    lower-degree family and $-1$ for our examples.

    Put $\mathbf p(\rho)=(1,\rho,\ldots,\rho^{n-2})^T$ and let
    $\mathsf W_A=(\rho(\mu_j)^k)_{j\in A,0\le k\le n-2}$.
    The vector $\boldsymbol{\rho}'_A=(\rho'(\mu_j))_{j\in A}$ contains the
    prescribed **unhalved** fixed ordinates. The fixed values of q give
    \[
    q(\rho)=\mathbf p(\rho)^T\mathsf W_A^{-1}\boldsymbol{\rho}'_A.
    \]
    Thus the inverse computes the known polynomial from its fixed-point values;
    it does not estimate Hamiltonian coefficients from a numerical trajectory.
    The exact rational calculation below must reproduce the original q.

    The reordered matrix and its Schur complement are
    \[
    \mathcal M_n=\begin{pmatrix}
    \mathsf W_A&\boldsymbol{\rho}'_A/2\\
    \mathbf p(\rho)^T&\rho'/2
    \end{pmatrix},\qquad
    \mathcal S=\frac{\rho'}2-
      \frac12\mathbf p(\rho)^T\mathsf W_A^{-1}\boldsymbol{\rho}'_A
      =\frac{\rho'-q(\rho)}2.
    \]
    Block elimination gives $\det\mathcal M_n=\det\mathsf W_A\,\mathcal S$.
    Consequently $D_n/\Delta_A=\varepsilon\mathcal S$. For our examples
    $V=D_n/\Delta_A=-\mathcal S$. Reversing the **fixed** ordinates gives
    the Schur complement $(\rho'+q)/2=U$, with the same column-order sign.
    The point $P_j^*=(\rho(\mu_j),-\rho'(\mu_j))$ is the image under the
    hyperelliptic sheet involution; the star does not denote complex conjugation.

    The exact checks compare the inverse expression with the original q,
    verify the block determinant formula and check the sign against the
    actual pole-ordered Onishi matrix. The numerical checks repeat these
    comparisons at every sampled moving point on both sheets.
    """)]+[mo.md(r"Genus "+str(row['genus'])+r": the inverse gives \[q(\rho)="+row['schur_q_latex']+r".\]") for row in exact_reports])
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 5. Resolve the composite sigma factor

    Put $E=\sum_{j\in A}w(\mu_j)=-\sum_{j\in B}w(\mu_j)$ and use exactly
    the paper's uncorrected factors $\bar u_j,\bar v_j$. Onishi's main identity
    can be grouped as
    \[
    D_n=K_A\frac{\sigma(w+E)}{\sigma_\sharp(w)}
                    \prod_{j\in A}\bar v_j,
    \]
    where $K_A$ depends only on fixed points. The composite factor must supply
    all residual points, rather than pretending to be one elementary mode.

    Apply the same identity to the $g$ **residual points on the opposite sheet** and the
    moving point. Its $g+1$ columns are $1,\rho,\ldots,\rho^g$, so its
    determinant is $\Delta_B\prod_{j\in B}(\rho-\rho(\mu_j))$.
    Using $\bar u_j\bar v_j=(-1)^g(\rho-\rho(\mu_j))$ and cancelling
    the common $\bar u_j$ factors proves
    \[
    \frac{\sigma(w+E)}{\sigma_\sharp(w)}
       =\Gamma_B\prod_{j\in B}\bar v_j,
    \qquad
    \Gamma_B=\frac{\Delta_B\prod_{j\in B}\sigma_\sharp(w(\mu_j))^g}
     {c_{g+1}\prod_{i<j\in B}\sigma_\flat(-w(\mu_i)+w(\mu_j))}.
    \]
    Cancellation is done on a regular open set; the resulting identity extends
    wherever the continued meromorphic expressions are defined. Consequently
    \[
    D_n=K\prod_j\bar v_j^{\beta_j},\qquad
    \prod_j\bar u_j^{\beta_j}=K_+(q+\rho'),
    \]
    with fixed nonzero constants. The second formula follows from multiplying
    the two-point identities and $q^2-F=4r$, or from the sheet-involution determinants.
    Shared membership of $P_1$ in $A$ and $B$ produces its squared point factor.

    The constants are prescribed, not calibrated. Here
    \[
    K_A=c_n(-1)^{n-1}
      \frac{\prod_{i<j\in A}\sigma_\flat(w(\mu_i)-w(\mu_j))}
           {\prod_{j\in A}\sigma_\sharp(w(\mu_j))^{n-1}},\quad K=K_A\Gamma_B.
    \]
    We use the sign corrected by Onishi in his published errata
    $c_n=(-1)^{(n-g-1)(n+g^2+2g)/2}$ for $n\ge g$.
    It applies without a moving normalization because the determinant ordinate
    is $\rho'/2$ and the Abel basis is $d\rho/\rho'=d\rho/(2(\rho'/2))$.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 6. Prove cancellation of the weighted exponential correction

    Directly from the scalar polynomial identity and the chain rule for the
    actual $Q$, we already have
    \[
    \frac{d}{dz}\log(q+\rho')=
    \sum_j\beta_j\left\{
       \frac{\rho'(\mu_j)-\rho'}{2[\rho(\mu_j)-\rho]}+B_j(\rho)\right\}.
    \]
    Differentiating the **Onishi-derived** complete $\bar u$ product instead
    gives the same left side with $B_j$ replaced by
    $\sum_k h_{jk}\rho^{k-1}$. Subtracting shows, on a nonconstant scalar branch,
    \[
    \sum_j\beta_j(c_{jk}-h_{jk})=0\quad(k=1,\ldots,g).
    \]
    This coefficient comparison is an analytic proof, not an inference from
    floating-point cancellation. It uses balanced fixed lifts; unbalanced
    lifts can introduce a holomorphic exponential multiplier in the raw product.

    Multiplying the paper's normalized mode formulas therefore gives
    \[
    U(z)=U(z_0)\frac{q(\rho(z))+\rho'(z)}{q(\rho(z_0))+\rho'(z_0)},\qquad
    V(z)=V(z_0)\frac{q(\rho(z))-\rho'(z)}{q(\rho(z_0))-\rho'(z_0)}.
    \]
    Their initial product conditions fix the constants. Combining with the interpolation in Section 3,
    \[
    U(z)=-\frac{D_n(P_j^*:j\in A;P(z))}{\Delta_A},\qquad
    V(z)=\frac{D_n(P_j:j\in A;P(z))}{\Delta_A}
    \]
    in the second interpolation family. In the lower-degree family the
    ordinate is the last column, so $U=+D_n^*/\Delta_A$ and $V=-D_n/\Delta_A$.
    This proves the determinant representation directly from the sigma modes,
    rather than simply restating the previously known $U,V$ scalar identities.
    """)
    return


@app.cell
def _(Curve, ctx, dps, kleinian_sigma_jet, mo, models, pole_basis):
    mo.md(r"""
    ## 7. Numerical checks of each link, with no fitting

    Fresh generapy curves are created for the two actual Hamiltonians. We check
    the weighted Abel sum against the period lattice **before** balancing one
    unit-weight lift. Solving for lattice coordinates identifies a discrete
    path class; its near-integer residual is checked independently. We do not
    alter a point's curve coordinates or fit a gauge. We then recompute its
    sigma normalizer and $h_{jk}$ on the chosen lift.

    Each test compares the matrix determinant with Onishi's corrected sigma
    expression, checks the auxiliary identity and composite factor, compares
    complete products with both $q\pm\rho'$ and the determinants, and checks
    the weighted exponential coefficients. Exterior points on both sheets
    are tested in each genus; genus three additionally uses the bounded branch.
    """)
    numerical_rows = []
    numerical_summaries = []
    tolerance = ctx.mpf(10)**(-(dps-8))
    def record(genus, check, actual, expected):
        error = abs(actual-expected)/max(1, abs(actual), abs(expected))
        numerical_rows.append(dict(genus=genus, check=check, scaled_error=ctx.nstr(error, 12), passed=bool(error < tolerance)))
    for _data in models:
        _g = int(_data['genus'])
        _N = len(_data['beta'])
        _M = sum(_data['beta'])
        _beta = _data['beta']
        def _number(value):
            _num, _den = value.as_numer_denom()
            return ctx.mpf(str(_num))/ctx.mpf(str(_den))
        _rho_symbol = _data['rho']
        _F_coefficients = [_number(_data['F'].coeff(_rho_symbol, k)) for k in range(2*_g+2)]
        _q_coefficients = [_number(_data['q'].coeff(_rho_symbol, k)) for k in range(_g+2)]
        def _F(coordinate):
            return ctx.polyval(list(reversed(_F_coefficients)), coordinate)
        def _q(coordinate):
            return ctx.polyval(list(reversed(_q_coefficients)), coordinate)
        _coordinates = [_number(value) for value in _data['a']]
        _ordinates = [_q(value) for value in _coordinates]
        _gamma = [_number(value) for value in _data['gamma']]
        _c = [[_number(value.coeff(_rho_symbol, k)) for k in range(_g)] for value in _data['B']]
        _curve = Curve({(0, 2): 1} | {(k, 0): -value for k, value in enumerate(_F_coefficients) if value}, ctx=ctx)
        _first = _curve.periods_kind_1()
        _second = _curve.periods_kind_2()
        _options = dict(omega=_first.omega, tau=_first.tau, kappa=_second.kappa,
                        characteristic=_curve.riemann_constant().characteristic,
                        normalization='hyperelliptic', ctx=ctx)
        def _jet(vector, order=2):
            return kleinian_sigma_jet(vector, order=order, **_options)
        _sharp = tuple(int(k == 1) for k in range(_g))
        _flat = (0, 0) if _g == 2 else (0, 0, 1)
        _zero = (0,)*_g
        _fixed_w = [_curve.abel_map_kind_1((a, b), reduce=False).value for a, b in zip(_coordinates, _ordinates)]
        _weighted_sum = sum((weight*w for weight, w in zip(_beta, _fixed_w)), ctx.zeros(_g, 1))
        _periods = ctx.matrix(_g, 2*_g)
        _periods[:, :_g], _periods[:, _g:] = 2*_first.omega, 2*_first.omega_prime
        _real_periods = ctx.matrix([[ctx.re(_periods[i, j]) for j in range(2*_g)] for i in range(_g)]+[[ctx.im(_periods[i, j]) for j in range(2*_g)] for i in range(_g)])
        _cycles = ctx.lu_solve(_real_periods, ctx.matrix([ctx.re(value) for value in _weighted_sum]+[ctx.im(value) for value in _weighted_sum]))
        _shift = ctx.matrix([ctx.nint(value) for value in _cycles])
        for _k in range(_g):
            record(_g, f'weighted divisor Abel lattice coordinate {_k+1}', (_weighted_sum-_periods*_shift)[_k], 0)
        assert _beta[-1] == 1
        _fixed_w[-1] -= _periods*_shift
        record(_g, 'balanced weighted Abel sum', ctx.norm(sum((weight*w for weight, w in zip(_beta, _fixed_w)), ctx.zeros(_g, 1))), 0)
        _normalizers = [_jet(w, 1)[_sharp] for w in _fixed_w]
        _h = []
        for _coordinate, _ordinate, _w in zip(_coordinates, _ordinates, _fixed_w):
            _j = _jet(_w, 2)
            _h2 = _j[(0, 2)+(0,)*(_g-2)]/(2*_j[_sharp])
            _row = [_j[(1, 1)+(0,)*(_g-2)]/_j[_sharp]-_j[(1,)+(0,)*(_g-1)]*_h2/_j[_sharp], _h2]
            if _g == 3:
                _row.append(_j[(0, 1, 1)]/_j[_sharp])
            _h.append(_row)
            # Second-kind primitives require the same changed path lift.
            # Their raw principal Abel-map values are therefore not substituted here.
        for _k in range(_g):
            record(_g, f'weighted exponential coefficient {_k+1}', sum(weight*(B[_k]-h[_k]) for weight, B, h in zip(_beta, _c, _h)), 0)
        _A = list(range(_g+2))
        _B = [0]+list(range(_g+2, _N))
        _n = _g+3
        _Delta_A = ctx.fprod(_coordinates[j]-_coordinates[i] for k, i in enumerate(_A) for j in _A[k+1:])
        _Delta_B = ctx.fprod(_coordinates[j]-_coordinates[i] for k, i in enumerate(_B) for j in _B[k+1:])
        _W_numeric = ctx.matrix([[_coordinates[j]**k for k in range(_n-1)] for j in _A])
        _fixed_vector_numeric = ctx.matrix([_ordinates[j] for j in _A])
        # Solving this fixed system is the numerical evaluation of W^{-1} times its data.
        _inverse_coefficients = ctx.lu_solve(_W_numeric, _fixed_vector_numeric)
        def _sign(size):
            return (-1)**((size-_g-1)*(size+_g*_g+2*_g)//2)
        _pair_A = ctx.fprod(_jet(_fixed_w[i]-_fixed_w[j], 1)[_flat] for k, i in enumerate(_A) for j in _A[k+1:])
        _pair_B = ctx.fprod(_jet(-_fixed_w[i]+_fixed_w[j], 1)[_flat] for k, i in enumerate(_B) for j in _B[k+1:])
        _Gamma = _Delta_B*ctx.fprod(_normalizers[j]**_g for j in _B)/(_sign(_g+1)*_pair_B)
        _K_A = _sign(_n)*(-1)**(_n-1)*_pair_A/ctx.fprod(_normalizers[j]**(_n-1) for j in _A)
        _K = _K_A*_Gamma
        _E = sum((_fixed_w[j] for j in _A), ctx.zeros(_g, 1))
        def _row(coordinate, ordinate, size):
            # Only polynomial columns are present in the auxiliary determinant.
            if size == _g+1:
                return [coordinate**k for k in range(size)]
            return [1]+[coordinate**k for k in range(1, _g+1)]+[ordinate/2, coordinate**(_g+1)]
        def _fs(points, abels):
            size = len(points)
            determinant = ctx.det(ctx.matrix([_row(a, b, size) for a, b in points]))
            numerator = _sign(size)*_jet(sum(abels, ctx.zeros(_g, 1)), 0)[_zero]
            numerator *= ctx.fprod(_jet(abels[i]-abels[j], 1)[_flat] for i in range(size) for j in range(i+1, size))
            sigma_expression = numerator/ctx.fprod(_jet(w, 1)[_sharp]**size for w in abels)
            return determinant, sigma_expression
        def _raw(w):
            normalizer = _jet(w, 1)[_sharp]
            return ([_jet(w+fixed, 1)[_flat]/(normalizer*fixed_normalizer) for fixed, fixed_normalizer in zip(_fixed_w, _normalizers)],
                    [_jet(w-fixed, 1)[_flat]/(normalizer*fixed_normalizer) for fixed, fixed_normalizer in zip(_fixed_w, _normalizers)])
        _positions = ('4.5', '4.6', '4.7') if _g == 2 else ('5.5', '5.6', '.462269035550823690056927809775')
        for _sheet in (1, -1):
            _start_rho = ctx.mpf(_positions[0])
            _start_y = _sheet*ctx.sqrt(_F(_start_rho))
            _start_w = _curve.abel_map_kind_1((_start_rho, _start_y), reduce=False).value
            _start_u, _start_v = _raw(_start_w)
            _target_U0 = (_q(_start_rho)+_start_y)/2
            _initial_u = [ctx.exp(ctx.log(_target_U0)/_M)]*_N
            _initial_v = [(offset-weight*_start_rho)/value for offset, weight, value in zip(_gamma, _beta, _initial_u)]
            for _text in _positions:
                _coordinate = ctx.mpf(_text)
                _ordinate = _sheet*ctx.sqrt(_F(_coordinate))
                _w = _curve.abel_map_kind_1((_coordinate, _ordinate), reduce=False).value
                _u_bar, _v_bar = _raw(_w)
                _tag = f'rho={_text}, sheet={_sheet:+d}'
                _D, _FS = _fs([(_coordinates[j], _ordinates[j]) for j in _A]+[(_coordinate, _ordinate)], [_fixed_w[j] for j in _A]+[_w])
                _D_star, _FS_star = _fs([(_coordinates[j], -_ordinates[j]) for j in _A]+[(_coordinate, _ordinate)], [-_fixed_w[j] for j in _A]+[_w])
                _aux, _aux_FS = _fs([(_coordinates[j], -_ordinates[j]) for j in _B]+[(_coordinate, _ordinate)], [-_fixed_w[j] for j in _B]+[_w])
                record(_g, f'main corrected Onishi identity: {_tag}', _D, _FS)
                record(_g, f'sheet-involution Onishi identity: {_tag}', _D_star, _FS_star)
                record(_g, f'auxiliary corrected Onishi identity: {_tag}', _aux, _aux_FS)
                _power_vector_numeric = ctx.matrix([_coordinate**k for k in range(_n-1)])
                _q_from_inverse_numeric = (_power_vector_numeric.T*_inverse_coefficients)[0]
                _schur_numeric = (_ordinate-_q_from_inverse_numeric)/2
                _block_numeric = ctx.matrix(
                    [[_coordinates[j]**k for k in range(_n-1)]+[_ordinates[j]/2] for j in _A]
                    +[[_coordinate**k for k in range(_n-1)]+[_ordinate/2]])
                record(_g, f'Vandermonde inverse gives prescribed q: {_tag}', _q_from_inverse_numeric, _q(_coordinate))
                record(_g, f'Schur block determinant: {_tag}', ctx.det(_block_numeric)/_Delta_A, _schur_numeric)
                record(_g, f'Onishi to Schur column sign: {_tag}', _D/_Delta_A, -_schur_numeric)
                record(_g, f'V interpolation: {_tag}', _D/_Delta_A, (_q(_coordinate)-_ordinate)/2)
                record(_g, f'U interpolation: {_tag}', -_D_star/_Delta_A, (_q(_coordinate)+_ordinate)/2)
                record(_g, f'composite residual split: {_tag}', _jet(_w+_E, 0)[_zero]/_jet(_w, 1)[_sharp], _Gamma*ctx.fprod(_v_bar[j] for j in _B))
                record(_g, f'raw weighted point product: {_tag}', _D, _K*ctx.fprod(value**weight for value, weight in zip(_v_bar, _beta)))
                _dw = _w-_start_w
                _u_values, _v_values = [], []
                for _j in range(_N):
                    _exponent = sum((a-b)*increment for a, b, increment in zip(_c[_j], _h[_j], _dw))
                    _u_values.append(_initial_u[_j]*_u_bar[_j]/_start_u[_j]*ctx.exp(_exponent))
                    _v_values.append(_initial_v[_j]*_v_bar[_j]/_start_v[_j]*ctx.exp(-_exponent))
                _U = ctx.fprod(value**weight for value, weight in zip(_u_values, _beta))
                _V = ctx.fprod(value**weight for value, weight in zip(_v_values, _beta))
                record(_g, f'sigma solution U = sheet-involution determinant: {_tag}', _U, -_D_star/_Delta_A)
                record(_g, f'sigma solution V = determinant: {_tag}', _V, _D/_Delta_A)
                record(_g, f'Hamiltonian product identity: {_tag}', _U+_V, _q(_coordinate))
                record(_g, f'paired mode powers: {_tag}', max(abs(a*b-(offset-weight*_coordinate)) for a, b, offset, weight in zip(_u_values, _v_values, _gamma, _beta)), 0)
        numerical_summaries.append(dict(genus=_g, main_size=_n, auxiliary_size=_g+1,
                                        lattice_shift=[int(value) for value in _shift],
                                        main_sign=_sign(_n), auxiliary_sign=_sign(_g+1),
                                        max_scaled_error=max(float(row['scaled_error']) for row in numerical_rows if row['genus'] == _g)))
    assert all(row['passed'] for row in numerical_rows)
    mo.ui.table(numerical_rows, page_size=15)
    return numerical_rows, numerical_summaries, tolerance


@app.cell
def _(mo):
    mo.md(r"""
    ## 8. What is complete, and what remains outside this proof?

    **Complete for both paper examples.** Interpolation, divisor multiplicity,
    main and auxiliary Onishi identities, the residual composite split, weighted
    gauge cancellation and exact initial normalization connect the sigma modes
    to the normalized $5\times5$ and $6\times6$ determinants. The double weight
    causes no missing residual mode and no repeated row in either determinant.

    **Scope conditions.** The two interpolation families require $\deg q\le g+1$.
    The fixed-point multiset must admit main and residual lists that are each
    internally distinct, with regular sigma normalizers and compatible balanced
    weighted Abel lifts. All these conditions hold in our examples. If the
    weight pattern forces repeated points inside either list, confluent Onishi
    limits and Hermite interpolation must be established separately. If the
    common divisor of all weights obstructs balancing identical mode lifts,
    quasi-period exponential multipliers must be retained instead.

    **Higher scalar degree.** The paper's mode-solution theorem is broader.
    For $\deg q>g+1$ the longer pole-order basis contains further $\rho^k y$
    columns. Divisor counting still suggests a generalized determinant, but
    nonzero interpolation rank and its sigma/composite normalization require
    another argument. This notebook does not quietly claim that extension.

    **What is still a modeling choice?** A determinant determines the scalar
    restriction $q$ on the invariant line. It does not determine the off-line
    polynomial $Q$ or its individual partial derivatives. Those phases remain
    supplied by the specified Hamiltonian, exactly as in the paper. This proof
    completes the product picture; it does not reconstruct arbitrary individual
    equations of motion from $q$ alone.

    Sources: [Onishi, Theorem 7.2 and Lemma 7.1](https://arxiv.org/pdf/math/0105189),
    using the author's [published errata](https://math-meijo-u.heteml.net/research/pub/correct_pem.pdf).
    This notebook derives and checks the weighted product result directly.
    """)
    return


@app.cell
def _(dps, exact_reports, example_dir, json, lower_family_checks, metadata, mo, numerical_rows, numerical_summaries, tolerance):
    report = dict(passed=all(row['passed'] for row in numerical_rows), dps=dps,
                  exact=exact_reports, lower_family_interpolation=lower_family_checks,
                  numerical=numerical_summaries, rows=numerical_rows, tolerance=str(tolerance),
                  scope='Complete weighted determinant proof for the two paper examples and stated interpolation families; higher degree and forced confluent lists excluded.',
                  fitting='None. FS signs, split constants and gauges are prescribed. Lattice rounding is verified discrete path identification.',
                  versions={package: metadata.version(package) for package in ('generapy', 'mpmath', 'sympy', 'marimo')})
    (example_dir/'results.json').write_text(json.dumps(report, indent=2)+'\n')
    mo.vstack([mo.md(f"**PASS:** {len(numerical_rows)} numerical checks at {dps} dps, plus exact algebra. Results are saved in `results.json`."), mo.ui.table(numerical_summaries)])
    return (report,)


if __name__ == '__main__':
    app.run()
