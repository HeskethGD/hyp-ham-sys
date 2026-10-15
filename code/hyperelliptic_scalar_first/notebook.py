import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import sys as _sys
    from pathlib import Path as _Path

    # Marimo can start with this folder, rather than the repository, on sys.path.
    # Resolve our sibling package from the notebook's location in either case.
    _project_root = str(_Path(__file__).resolve().parent.parent)
    if _project_root not in _sys.path:
        _sys.path.insert(0, _project_root)

    import marimo as mo
    import sympy as sp
    from hyperelliptic_scalar_first.model import example, genus_three_example
    from hyperelliptic_scalar_first.verification import run_all, SigmaExample, rk4

    return SigmaExample, example, genus_three_example, mo, rk4, run_all, sp


@app.cell
def _(mo):
    mo.md(r"""
    # Solve the scalar equation, then reconstruct the Hamiltonian modes

    This supporting notebook solves for $\rho(z)$, integrates the individual
    mode equations, and verifies the polynomial canonical system. The propagation
    variable is $z$. The worked-example notebooks provide the longer trajectories
    and figures used in the paper.

    We assume the length and mode-power normalization has already been made so
    that the scalar polynomial has leading coefficient $4$. We use $w$ throughout.
    A prime on $\rho$ or a mode means differentiation with respect to $z$.
    In fixed curve and polynomial formulas we may write $\rho$ without its argument.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 1. Set up the scalar equation

    **Aim:** specify the data we need before introducing sigma functions.
    Choose the actual polynomial $Q(x_1,\ldots,x_N)$, positive integer degeneracies
    $\beta_j$, and constants $\gamma_j$ with $\sum_j\gamma_j=0$. The powers are
    \[
    x_j(z)=u_j(z)v_j(z)=\gamma_j-\beta_j\rho(z).
    \]
    Restrict $Q$ and the weighted product to this line:
    \[
    q(\rho)=Q(\gamma_1-\beta_1\rho,\ldots,\gamma_N-\beta_N\rho),\qquad
    r(\rho)=\prod_j(\gamma_j-\beta_j\rho)^{\beta_j}.
    \]
    Our scalar equation is
    \[
    \boxed{[\rho'(z)]^2=F(\rho(z)),\qquad F(\rho)=q(\rho)^2-4r(\rho).}
    \]
    We assume $F$ is squarefree of degree $2g+1$ with leading coefficient $4$.
    The fixed values $\gamma_j/\beta_j$ are distinct and
    $q(\gamma_j/\beta_j)\ne0$. The polynomial corrections $B_j$ defined in
    Section 3 must have degree at most $g-1$. We check that condition directly;
    it also follows from $\deg Q\le g$, but higher-degree $Q$ can satisfy it by cancellation.

    The paper's centred power convention gives
    \[
    \rho(z)=-\frac{\sum_jx_j(z)}{\sum_j\beta_j}.
    \]
    When defining the examples' polynomials in independent powers, we write this
    linear function as $\rho(x)=-\sum_jx_j/\sum_j\beta_j$. Its value on the
    power line is $\rho(z)$; it is not another evolving unknown.

    Choose an initial point $(\rho(z_0),\rho'(z_0))$ on the curve and consistent
    integration paths. The formulas below apply locally where their denominators
    and modes are nonzero. Singular curves, branch fixed points and even-degree
    curves need additional charts. Positive integer $\beta_j$, including repeated
    modes, require no change to this construction under these assumptions.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 2. Solve for $\rho(z)$ using Abel integrals and sigma

    ### 2.1. Put the scalar quadrature into $w$

    **Aim:** choose coordinates whose first component measures propagation and
    whose other components integrate powers of $\rho(z)$.
    On the normalized curve
    \[
    y^2=F(\rho)=4\rho^{2g+1}+\cdots,
    \]
    the moving point is $(\rho(z),y(z))$, with $y(z)=\rho'(z)$.
    The standard holomorphic differentials and their integrals are
    \[
    \omega_k=\frac{\rho^{k-1}d\rho}{y},\qquad
    w_k(z)=\int_\infty^{(\rho(z),\rho'(z))}\omega_k,\qquad k=1,\ldots,g.
    \]
    “Holomorphic” means these differentials have no poles on the compact curve.
    This is the standard basis for a smooth odd-degree hyperelliptic curve
    ([Enolski et al., Section 3.1](https://arxiv.org/pdf/1011.6459#page=6)).
    Choosing paths gives actual values of $w(z)$ rather than values modulo periods;
    this is all that “lifted Abel integrals” means here.

    Substitution of $y(z)=\rho'(z)$ gives, directly,
    \[
    \boxed{w_1(z)-w_1(z_0)=z-z_0,\qquad w_k'(z)=\rho(z)^{k-1}.}
    \]
    Thus $w(z)$ is the Abel image of **one curve point**, not an arbitrary vector
    in $\mathbb C^g$. We only need this one-point image; no hierarchy of subsets
    is needed for the construction.

    ### 2.2. Recover $\rho(z)$ from sigma

    **Aim:** invert the quadrature while keeping $w$ on that one-point image.
    The explicit inversions used in our tests are
    \[
    \begin{array}{c|c|c}
    g & \rho(z) & \text{conditions on }w(z)\\ \hline
    2 & -\sigma_1(w(z))/\sigma_2(w(z)) & \sigma(w(z))=0\\
    3 & -\sigma_{13}(w(z))/\sigma_{23}(w(z)) & \sigma(w(z))=\sigma_3(w(z))=0.
    \end{array}
    \]
    Here sigma subscripts mean partial derivatives with respect to the indicated
    components of $w$. The first coordinate is fixed by $z-z_0$; the other
    coordinates follow the chosen curve branch. They do not all evolve linearly.
    The general scalar solution is the inverse of the displayed Abel quadrature;
    these sigma formulas make that inverse explicit in genera two and three.
    The arbitrary-genus derivative formula in the cited preprint is conjectural,
    so we do not use it as a theorem.

    For the mode ratios we use Ônishi's two sigma derivatives:
    \[
    \sigma_\sharp(w)=\left(\prod_{\substack{2\le k\le g\\k\ \mathrm{even}}}
    \partial_{w_k}\right)\sigma(w),\qquad
    \sigma_\flat(w)=\left(\prod_{\substack{3\le k\le g\\k\ \mathrm{odd}}}
    \partial_{w_k}\right)\sigma(w).
    \]
    An empty product means no derivative. In genus two these are
    $(\sigma_\sharp,\sigma_\flat)=(\sigma_2,\sigma)$; in genus three they are
    $(\sigma_2,\sigma_3)$. These are the functions from
    [Ônishi, Section 6](https://arxiv.org/pdf/math/0105189#page=20).
    [Enolski et al., Section 4.2](https://arxiv.org/pdf/1011.6459#page=12)
    explicitly credit Ônishi for the genus-three inversion. Generapy uses the
    corresponding hyperelliptic sigma normalization. We work where the displayed
    sigma denominators are nonzero.

    ### 2.3. Fix the mode points and their differentials

    **Aim:** identify now the fixed points and differentials needed to integrate
    the modes later. Choose continued propagation labels $\mu_j$ so that
    \[
    \boxed{\rho(\mu_j)=\gamma_j/\beta_j,\qquad
    \rho'(\mu_j)=q(\rho(\mu_j)).}
    \]
    These points lie on the curve because $r(\rho(\mu_j))=0$.
    Their Abel vectors $w(\mu_j)$ are constant in $z$. The labels refer to chosen
    analytic continuations; they can be complex. For calculation we can obtain
    $w(\mu_j)$ directly from the fixed curve points.

    In addition to the holomorphic differentials above, we will use
    \[
    \boxed{\omega_j^{\mathrm{mode}}=
    \frac{\rho'(\mu_j)-y}{2[\rho(\mu_j)-\rho]}\frac{d\rho}{y}.}
    \]
    It has a simple pole of residue $+1$ at
    $(\rho(\mu_j),-\rho'(\mu_j))$ and residue $-1$ at infinity.
    Such a differential is called a differential of the third kind.
    Along the moving point it is exactly
    \[
    \omega_j^{\mathrm{mode}}=
    \frac{\rho'(\mu_j)-\rho'(z)}{2[\rho(\mu_j)-\rho(z)]}\,dz.
    \]
    To continue the scalar solution through a turning point, use the regular system
    $\rho'(z)=y(z)$, $y'(z)=F_\rho(\rho(z))/2$. The numerical tests below use
    short paths away from turning points.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 3. Split the mode equation

    **Aim:** isolate the fixed-point ratio we can integrate by sigma, leaving a
    polynomial we can integrate using $w_k'(z)=\rho(z)^{k-1}$.
    The paper's equation (2.22), with $P=1$, gives
    \[
    \frac{d}{dz}\log u_j(z)=
    \frac{\beta_jq(\rho(z))}{2x_j(z)}
    -\left.\partial_{x_j}Q\right|_{x=x(z)}
    +\frac12\frac{d}{dz}\log x_j(z).
    \]
    Use $x_j(z)=\beta_j[\rho(\mu_j)-\rho(z)]$ and
    $\rho'(\mu_j)=q(\rho(\mu_j))$ to write this as
    \[
    \boxed{\frac{d}{dz}\log u_j(z)=
    \frac{\rho'(\mu_j)-\rho'(z)}{2[\rho(\mu_j)-\rho(z)]}
    +B_j(\rho(z)),}
    \]
    \[
    \boxed{B_j(\rho(z))=
    \frac{q(\rho(z))-q(\rho(\mu_j))}{2[\rho(\mu_j)-\rho(z)]}
    -\left.\partial_{x_j}Q\right|_{x=x(z)}
    =\sum_{k=1}^g c_{jk}\rho(z)^{k-1}.}
    \]
    The divided difference is a polynomial: its apparent singularity is removable.
    The $c_{jk}$ are its constant coefficients after including the actual partial
    of $Q$. Differentiate $Q$ before restricting to the power line; $q$ alone
    does not determine the individual phases.

    Since $v_j(z)=x_j(z)/u_j(z)$, its equation is
    \[
    \boxed{\frac{d}{dz}\log v_j(z)=
    \frac{-\rho'(\mu_j)-\rho'(z)}{2[\rho(\mu_j)-\rho(z)]}
    -B_j(\rho(z)).}
    \]
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 4. Integrate the two parts and obtain the modes

    ### 4.1. Integrate the difference ratio

    **Aim:** find the point factors defined by
    \[
    \boxed{\frac{d}{dz}\log\widehat u_j(z)=
    \frac{\rho'(\mu_j)-\rho'(z)}{2[\rho(\mu_j)-\rho(z)]},\qquad
    \frac{d}{dz}\log\widehat v_j(z)=
    \frac{-\rho'(\mu_j)-\rho'(z)}{2[\rho(\mu_j)-\rho(z)]}.}
    \]
    The standard third-kind integration formula supplies the exponential
    coefficients from the **fixed mode point**. For genera two and three,
    $\sigma_\sharp=\sigma_2$. Define
    \[
    \zeta_{\sharp,k}(w)=\partial_{w_k}\log\sigma_\sharp(w).
    \]
    Then the required constants have the short expressions
    \[
    \boxed{h_{j2}=\tfrac12\zeta_{\sharp,2}(w(\mu_j)),\qquad
    h_{j1}=\zeta_{\sharp,1}(w(\mu_j))
    -\frac{\sigma_1(w(\mu_j))}{\sigma_2(w(\mu_j))}h_{j2},}
    \]
    \[
    \boxed{h_{j3}=\zeta_{\sharp,3}(w(\mu_j))\quad\text{in genus three}.}
    \]
    These short sigma expressions apply to genera two and three. In arbitrary
    genus the same integration construction uses the standard fixed-point
    second-kind Abel integrals, whose explicit differentials and infinity
    normalization are given in Section 9.1. A uniform short sigma-only formula
    in every genus is not established here.

    Everything here is evaluated at the constant $w(\mu_j)$.
    In genus two the middle ratio is $-\rho(\mu_j)$, so equivalently
    $h_{j1}=\zeta_{\sharp,1}(w(\mu_j))+\rho(\mu_j)h_{j2}$.
    These are minus the standard second-kind Abel integrals at that point,
    with their infinity normalization. The genus-three expressions follow
    from [Enolski et al., Proposition III.2, Eq. III.10](https://arxiv.org/pdf/1106.2408#page=14);
    the same paper points to the genus-two treatment in Enolskii–Pronine–Richter.
    This fixed-point form uses our infinity-based sharp/flat convention;
    Eq. III.17 presents the third-kind identity with shifted ordinary sigma.
    The normalization and the alternative reference-point expressions are
    explained at the end. The tests check their agreement independently.

    The point factors are therefore
    \[
    \boxed{\widehat u_j(z)=\sqrt{\beta_j}
    \frac{\sigma_\flat(w(z)+w(\mu_j))}
    {\sigma_\sharp(w(z))\sigma_\sharp(w(\mu_j))}
    \exp\!\left[-\sum_{k=1}^g h_{jk}(w_k(z)-w_k(z_0))\right],}
    \]
    \[
    \boxed{\widehat v_j(z)=(-1)^{g+1}\sqrt{\beta_j}
    \frac{\sigma_\flat(w(z)-w(\mu_j))}
    {\sigma_\sharp(w(z))\sigma_\sharp(w(\mu_j))}
    \exp\!\left[+\sum_{k=1}^g h_{jk}(w_k(z)-w_k(z_0))\right].}
    \]
    The sigma ratios alone have an additional holomorphic logarithmic differential.
    The fixed-point coefficients above cancel it. The longer finite-derivative
    formulas at the end provide a second way to obtain the same constants.
    **Here the hats include this exponential:** their logarithmic derivatives
    are exactly the difference ratios above. Their product is
    \[
    \widehat u_j(z)\widehat v_j(z)=
    \beta_j[\rho(\mu_j)-\rho(z)]=x_j(z).
    \]

    The integral and sigma expressions are the same solution. Along the chosen
    curve path, with logarithms continued from $z_0$,
    \[
    \boxed{\int_{(\rho(z_0),\rho'(z_0))}^{(\rho(z),\rho'(z))}
    \omega_j^{\mathrm{mode}}
    =\log\frac{\widehat u_j(z)}{\widehat u_j(z_0)}.}
    \]
    For comparison with the earlier third-kind integral, we also have
    \[
    \boxed{\frac{\rho'(\mu_j)}2
    \int_{(\rho(z_0),\rho'(z_0))}^{(\rho(z),\rho'(z))}
    \frac{d\rho}{[\rho(\mu_j)-\rho]y}
    =\frac12\log\!\left[
    \frac{\widehat u_j(z)/\widehat v_j(z)}
    {\widehat u_j(z_0)/\widehat v_j(z_0)}\right].}
    \]
    Substituting the displayed hats gives the sigma ratio with its exponential
    correction explicitly. This is the form of the third-kind sigma integral
    treated by [Enolski et al., Equation III.17](https://arxiv.org/pdf/1106.2408#page=15),
    here written with our infinity-based sharp/flat factors and matched
    second-kind normalization. The two integrals agree because
    \[
    \omega_j^{\mathrm{mode}}=
    \tfrac12d\log[\rho(\mu_j)-\rho]
    +\frac{\rho'(\mu_j)}{2[\rho(\mu_j)-\rho]}\frac{d\rho}{y}.
    \]
    The product identity for the hats supplies the first logarithm. Both integral
    forms are checked independently against the sigma expressions below.

    ### 4.2. Integrate the polynomial part

    **Aim:** use the coordinates already introduced, without another special function.
    Since $B_j(\rho(z))=\sum_kc_{jk}\rho(z)^{k-1}$,
    \[
    \boxed{\int B_j(\rho(z))\,dz=\sum_{k=1}^g c_{jk}w_k(z)+\text{constant}.}
    \]
    Thus its definite integral is $\sum_kc_{jk}[w_k(z)-w_k(z_0)]$.
    The degree bound is exactly what makes holomorphic Abel integrals sufficient;
    higher powers would require additional differentials.

    ### 4.3. Assemble $u_j(z)$ and $v_j(z)$

    **Aim:** combine the point factors with that polynomial integral and the desired
    initial amplitudes. The solutions are
    \[
    \boxed{u_j(z)=u_j(z_0)\frac{\widehat u_j(z)}{\widehat u_j(z_0)}
    \exp\!\left[\sum_{k=1}^g c_{jk}(w_k(z)-w_k(z_0))\right],}
    \]
    \[
    \boxed{v_j(z)=v_j(z_0)\frac{\widehat v_j(z)}{\widehat v_j(z_0)}
    \exp\!\left[-\sum_{k=1}^g c_{jk}(w_k(z)-w_k(z_0))\right].}
    \]
    Choose $u_j(z_0)v_j(z_0)=\gamma_j-\beta_j\rho(z_0)$. Then
    \[
    u_j(z)v_j(z)=\gamma_j-\beta_j\rho(z),\qquad
    v_j(z)=\frac{\gamma_j-\beta_j\rho(z)}{u_j(z)}.
    \]
    Equivalently, insert the sigma ratios directly in the final formulas and use
    $c_{jk}-h_{jk}$ in the exponential for $u_j$, and its negative for $v_j$.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 5. Apply the paper's theorem to recover the coupled system

    **Aim:** select the initial amplitudes so that these individual solutions
    satisfy the desired polynomial Hamiltonian equations.
    Impose the paper's initial product condition
    \[
    \boxed{\prod_j u_j(z_0)^{\beta_j}=
    \frac{q(\rho(z_0))+\rho'(z_0)}2.}
    \]
    The paper's theorem then gives, along the solution,
    \[
    \prod_j u_j(z)^{\beta_j}=\frac{q(\rho(z))+\rho'(z)}2,\qquad
    \prod_j v_j(z)^{\beta_j}=\frac{q(\rho(z))-\rho'(z)}2.
    \]
    Substitution in the mode equations produces
    \[
    \boxed{u_j'(z)=\beta_jv_j(z)^{\beta_j-1}
    \prod_{m\ne j}v_m(z)^{\beta_m}
    -\left.\partial_{x_j}Q\right|_{x=x(z)}u_j(z),}
    \]
    \[
    \boxed{v_j'(z)=-\beta_ju_j(z)^{\beta_j-1}
    \prod_{m\ne j}u_m(z)^{\beta_m}
    +\left.\partial_{x_j}Q\right|_{x=x(z)}v_j(z).}
    \]
    These are the canonical equations for the polynomial Hamiltonian
    \[
    \boxed{H(u,v)=\prod_j u_j^{\beta_j}+\prod_jv_j^{\beta_j}
    -Q(u_1v_1,\ldots,u_Nv_N),\qquad
    u_j'=\partial_{v_j}H,\quad v_j'=-\partial_{u_j}H.}
    \]
    Our normalization selects $H(u(z),v(z))=0$. Other energy levels follow by
    replacing $Q$ with $Q+E$. The polynomial equations also remain meaningful
    at a zero mode, where the logarithmic formulas need another chart.

    Finally the power relations give the conserved quantities directly:
    \[
    \beta_Nx_j(z)-\beta_jx_N(z)=\beta_N\gamma_j-\beta_j\gamma_N,\qquad j<N.
    \]
    Together with $H$, these are the commuting invariants from the paper. We stop
    at this Hamiltonian construction.
    """)
    return


@app.cell
def _(example):
    data = example()
    return (data,)


@app.cell
def _(data, mo, sp):
    _rho_z = sp.Function('rho')(sp.Symbol('z'))
    mo.md(rf"""
    ## 6. Test a normalized genus-two example with a double mode

    **Aim:** check the scalar inverse, both mode factors and unequal polynomial
    phase corrections in the same example. Choose
    $\beta=(2,1,1,1,1)$ and
    $\gamma=(-3/5,-1/5,-1/10,2/5,1/2)$. Let
    \[
    \rho(x)=-\frac{{\sum_{{j=1}}^5 x_j}}{{6}},\qquad
    Q(x)=4(\rho(x)+27/40)(\rho(x)+1/20)(\rho(x)-3/5)
    +\frac{{1}}{{3}}(x_1/2-x_2+1/10)x_3.
    \]
    The extra term vanishes on the power line but has nonzero individual partials.
    On that line,
    \[
    q(\rho(z))={sp.latex(sp.factor(data['q']).subs(data['rho'], _rho_z))},
    \]
    \[
    r(\rho(z))={sp.latex(data['r'].subs(data['rho'], _rho_z))}.
    \]
    The degree-six terms of $q^2-4r$ cancel, leaving a squarefree degree-five
    polynomial with leading coefficient $4$. The degree-two terms of $B_j$
    also cancel, leaving affine corrections. These are exact algebraic checks,
    not an extra rescaling of the curve used by sigma.

    The fixed points satisfy
    $(\rho(\mu_1),\ldots,\rho(\mu_5))=(-3/10,-1/5,-1/10,2/5,1/2)$
    and $\rho'(\mu_j)=q(\rho(\mu_j))\ne0$.
    The double mode gives a double factor in $r$, while $F$ remains smooth.
    """)
    return


@app.cell
def _(data, mo, sp):
    _rho_z = sp.Function('rho')(sp.Symbol('z'))
    mo.ui.table([{'mode j': j+1, 'B_j(rho(z))': str(sp.factor(B).subs(data['rho'], _rho_z))}
                 for j, B in enumerate(data['B'])], selection=None)
    return


@app.cell
def _(genus_three_example, mo, sp):
    data3 = genus_three_example()
    _rho_z = sp.Function('rho')(sp.Symbol('z'))
    mo.md(rf"""
    ## 7. Test genus three and the $w_3(z)$ contribution

    **Aim:** test quadratic corrections and the genus-three sigma numerator.
    Choose $\beta=(2,1,1,1,1,1,1)$ and
    $\gamma=(-3/5,-1/4,-1/5,-1/10,1/5,2/5,11/20)$. Set
    \[
    \rho(x)=-\frac{{\sum_{{j=1}}^7x_j}}{{8}},
    \]
    \[
    Q(x)=4(\rho(x)+37/40)(\rho(x)+3/20)(\rho(x)-3/10)(\rho(x)-13/20)
    +\frac{{1}}{{7}}(x_1/2-x_2+1/20)x_3^2.
    \]
    The scalar restriction is
    \[
    q(\rho(z))={sp.latex(sp.factor(data3['q']).subs(data3['rho'], _rho_z))},
    \]
    \[
    r(\rho(z))={sp.latex(data3['r'].subs(data3['rho'], _rho_z))}.
    \]
    Here the degree-eight terms cancel, leaving a smooth degree-seven curve
    with leading coefficient $4$. The cubic terms in $B_j$ cancel too.
    Consequently each $B_j(\rho(z))=c_{{j1}}+c_{{j2}}\rho(z)+c_{{j3}}\rho(z)^2$
    integrates to $c_{{j1}}w_1(z)+c_{{j2}}w_2(z)+c_{{j3}}w_3(z)$.
    All fixed points are distinct and nonbranch. This example checks the
    $\sigma_3$ numerator, scalar inversion and all seven canonical pairs.
    """)
    return (data3,)


@app.cell
def _(data3, mo, sp):
    _rho_z = sp.Function('rho')(sp.Symbol('z'))
    mo.ui.table([{'mode j': j+1, 'B_j(rho(z))': str(sp.factor(B).subs(data3['rho'], _rho_z))}
                 for j, B in enumerate(data3['B'])], selection=None)
    return


@app.cell
def _(run_all):
    import json as _json
    from pathlib import Path as _ReportPath

    report = run_all(dps=30)
    _ReportPath(__file__).with_name('verification_results.json').write_text(
        _json.dumps(report, indent=2)+'\n')
    return (report,)


@app.cell
def _(mo, report):
    mo.md(f"""
    ## 8. Verify the construction independently

    Exact algebra and notebook source: **{len(report['symbolic']['checks'])} checks,
    passed = {report['symbolic']['passed']}**.
    Genus two: **{len(report['numerical']['rows'])} checks at 30 digits,
    passed = {report['numerical']['passed']}**.
    Genus three: **{len(report['genus_three']['rows'])} checks at 30 digits,
    passed = {report['genus_three']['passed']}**.

    These checks cover scalar inversion, prescribed propagation values,
    independent Abel and third-kind quadratures, sigma products and derivatives,
    polynomial integrals and every canonical equation on both sheets.
    The short fixed-point formulas for h are checked against direct regularized
    second-kind Abel integrals, the reference-point derivative formulas, and
    finite sigma Taylor expansions at two reference points. The mode solutions
    use the short fixed-point formulas. No coefficient fitting is used.

    The plots compare each individual sigma solution with direct integration
    of the full polynomial Hamiltonian system. RK4 step halving checks convergence.
    Numerical validation is for these two examples and short regular paths;
    it does not establish global continuation or arbitrary-genus inversion.
    """)
    return


@app.cell
def _(SigmaExample, rk4):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from pathlib import Path

    plot_model = SigmaExample(dps=30)
    _c = plot_model.ctx
    _start = plot_model.evaluate(_c.mpf('4.5'))
    _initial = plot_model.initial(_start)
    _final_length = plot_model.integral(lambda coordinate: 1/_c.sqrt(plot_model.F(coordinate)), _start['rho'], _c.mpf('4.65'))
    _z = []; _sigma = []; _rk = []
    for _k in range(9):
        _length = _final_length*_k/8
        _end = plot_model.point_at_z(_length, _start['rho'])
        _z.append(float(_length))
        _sigma.append(plot_model.solution(_start, _end, _initial))
        _rk.append(rk4(plot_model.field, _initial, _length, 64))
    fig, _axes = plt.subplots(2, 3, figsize=(12, 6), constrained_layout=True)
    for _j, _axis in enumerate(_axes.flat):
        if _j >= plot_model.N:
            _axis.set_visible(False)
            continue
        _axis.plot(_z, [float(row[_j].real) for row in _sigma], label='sigma solution')
        _axis.plot(_z, [float(row[_j].real) for row in _rk], 'o', markersize=4, label='full Hamiltonian RK4')
        _axis.set(xlabel='$z-z_0$', ylabel=rf'$\mathrm{{Re}}\,u_{{{_j+1}}}(z)$')
        _axis.xaxis.set_major_locator(MaxNLocator(4))
        _axis.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
        _axis.grid(alpha=.2)
    _axes[0, 0].legend(fontsize=8)
    fig.suptitle('Degenerate genus-two example: individual mode trajectories')
    fig.savefig(Path(__file__).with_name('trajectory.png'), dpi=160)
    return (fig,)


@app.cell
def _(fig):
    fig
    return


@app.cell
def _(SigmaExample, rk4):
    import matplotlib.pyplot as _plt
    from matplotlib.ticker import MaxNLocator as _MaxNLocator
    from pathlib import Path as _PlotPath

    _model = SigmaExample(dps=30, genus=3)
    _c = _model.ctx
    _start = _model.evaluate(_model.rho0)
    _initial = _model.initial(_start)
    _length_final = _model.integral(lambda coordinate: 1/_c.sqrt(_model.F(coordinate)), _start['rho'], _model.sample_rho[-1])
    _z = []; _sigma = []; _rk = []
    for _k in range(9):
        _length = _length_final*_k/8
        _end = _model.point_at_z(_length, _start['rho'])
        _z.append(float(_length))
        _sigma.append(_model.solution(_start, _end, _initial))
        _rk.append(rk4(_model.field, _initial, _length, 64))
    fig3, _axes = _plt.subplots(2, 4, figsize=(14, 6), constrained_layout=True)
    for _j, _axis in enumerate(_axes.flat):
        if _j >= _model.N:
            _axis.set_visible(False)
            continue
        _axis.plot(_z, [float(row[_j].real) for row in _sigma], label='sigma solution')
        _axis.plot(_z, [float(row[_j].real) for row in _rk], 'o', markersize=4, label='full Hamiltonian RK4')
        _axis.set(xlabel='$z-z_0$', ylabel=rf'$\mathrm{{Re}}\,u_{{{_j+1}}}(z)$')
        _axis.xaxis.set_major_locator(_MaxNLocator(4))
        _axis.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
        _axis.grid(alpha=.2)
    _axes[0, 0].legend(fontsize=8)
    fig3.suptitle('Degenerate genus-three example: individual mode trajectories')
    fig3.savefig(_PlotPath(__file__).with_name('trajectory_genus3.png'), dpi=160)
    return (fig3,)


@app.cell
def _(fig3):
    fig3
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 9. Obtain the constants $h_{jk}$ explicitly

    **Aim:** finish the sigma solution with a finite formula for its exponential
    correction. No fitting or numerical differentiation is involved.

    Define the zeta-like functions by their sigma logarithmic derivatives:
    \[
    \zeta_{\flat,k}(w)=\partial_{w_k}\log\sigma_\flat(w),\qquad
    \zeta_{\sharp,k}(w)=\partial_{w_k}\log\sigma_\sharp(w).
    \]
    On the trajectory their differences are
    \[
    \boxed{Z_{jk}(z)=\zeta_{\flat,k}(w(z)+w(\mu_j))
    -\zeta_{\sharp,k}(w(z)).}
    \]
    These $Z_{jk}(z)$ are **functions of $z$**, not constants. The sharp derivative
    is needed because $\sigma(w(z))$ vanishes on the one-point image for $g\ge2$.
    The chain rule gives the logarithmic derivative of the **sigma ratio alone**:
    \[
    \frac{d}{dz}\log\!\left[
    \frac{\sigma_\flat(w(z)+w(\mu_j))}
    {\sigma_\sharp(w(z))\sigma_\sharp(w(\mu_j))}\right]
    =\sum_{k=1}^g\rho(z)^{k-1}Z_{jk}(z).
    \]
    Its difference from the ratio we wanted to integrate is
    \[
    \boxed{\sum_{k=1}^g\rho(z)^{k-1}Z_{jk}(z)
    -\frac{\rho'(\mu_j)-\rho'(z)}{2[\rho(\mu_j)-\rho(z)]}
    =\sum_{k=1}^g h_{jk}\rho(z)^{k-1}.}
    \]
    The two logarithmic differentials have the same poles and residues; their
    difference is holomorphic. Expanding it in the basis $\omega_k$ proves that
    the $h_{jk}$ are constants and that the right side has degree at most $g-1$.
    This is also why the exponential in Section 4.1 removes precisely this term.

    ### 9.1. The normalization behind the fixed-point formula

    The second-kind differentials associated with this sigma normalization are
    fixed by the curve coefficients. Write $F(\rho)=\sum_{m=0}^{2g+1}\lambda_m\rho^m$,
    with $\lambda_{2g+1}=4$, and put coefficients outside this range equal to zero.
    Then
    \[
    dr_k=\frac{1}{4y}
    \sum_{m=k}^{2g+1-k}(m+1-k)\lambda_{m+1+k}\rho^m\,d\rho.
    \]
    These are [Enolski et al., Eq. II.2](https://arxiv.org/pdf/1106.2408#page=4),
    using the same basis as generapy. The third-kind sigma formula, restricted
    to the one-point image in our sharp/flat convention, gives
    \[
    \boxed{h_{jk}=-\operatorname{Reg}\int_\infty^{(\rho(\mu_j),\rho'(\mu_j))}dr_k.}
    \]
    One can obtain this relation from the
    [sharp/flat sigma-kernel identity of Gibbons–Matsutani–Ônishi, Theorem 7.17](https://arxiv.org/pdf/1204.3747#page=19)
    and Enolski's fundamental differential identity, Eq. II.11. Differentiate
    with respect to the fixed curve point: the holomorphic coefficients satisfy
    $dh_{jk}=-dr_k$. This variation is of the chosen fixed point, not of $z$.
    The sheet involution makes both primitives odd, fixing the additive constant
    to zero in the infinity normalization.

    “Reg” means subtract the divergent powers of the local coordinate
    $\rho^{-1/2}$ at infinity and retain the finite part. We choose the primitive
    with zero constant term there, and use the same lifted paths for these
    integrals as for $w(\mu_j)$. This fixes the exponential coefficients
    without choosing a moving reference point. The short expressions in
    Section 4.1 are sigma representations of these primitives in genera two
    and three; the odd hyperelliptic infinity normalization removes the
    finite-base integration constants in Eq. III.10.

    Thus the literature gives an intrinsic fixed-point prescription. The
    reference-point formulas below are a separate finite-derivative evaluation
    of the same $h_{jk}$. Our three-way numerical comparison also checks the
    signs, normalization and paths used in this translation.

    ### 9.2. A finite expression in any genus

    Take a regular reference point with curve coordinate $\rho=0$ and
    $y(0)^2=F(0)$, avoiding fixed mode points and zero sigma denominators.
    For this calculation only, write $w(\rho)$ for the same Abel integrals
    parametrized by the curve coordinate, so that
    \[
    \frac{dw_k(\rho)}{d\rho}=\frac{\rho^{k-1}}{y(\rho)},\qquad y(\rho)^2=F(\rho).
    \]
    Then the coefficients are given explicitly by
    \[
    \boxed{h_{jk}=\frac{1}{(k-1)!}
    \left.\frac{d^{k-1}}{d\rho^{k-1}}
    \left\{\sum_{\ell=1}^g\rho^{\ell-1}
    [\zeta_{\flat,\ell}(w(\rho)+w(\mu_j))
    -\zeta_{\sharp,\ell}(w(\rho))]
    -\frac{\rho'(\mu_j)-y(\rho)}{2[\rho(\mu_j)-\rho]}\right\}\right|_{\rho=0}.}
    \]
    This extracts the coefficients of a polynomial of known degree, using only
    finitely many sigma derivatives. If zero is unsuitable, take its finite Taylor
    expansion at another regular coordinate and translate that polynomial back
    into powers of $\rho$. The same constants result.

    ### 9.3. Expanded expressions in genera two and three

    To make these directly evaluable, a star denotes the fixed reference point
    $\rho=0$: $y_*=y(0)$, $w_*=w(0)$, and
    \[
    Z_{jk}^*=\zeta_{\flat,k}(w_*+w(\mu_j))-\zeta_{\sharp,k}(w_*).
    \]
    In the formulas below, $(\partial_\ell Z_{jk})_*$ means differentiate this
    zeta difference with respect to the moving Abel component $w_\ell$ and
    then evaluate at $w_*$. It is not a propagation derivative. The needed curve
    derivatives are already explicit:
    \[
    y_\rho(0)=\frac{F_\rho(0)}{2y_*},\qquad
    y_{\rho\rho}(0)=\frac{F_{\rho\rho}(0)}{2y_*}
    -\frac{F_\rho(0)^2}{4y_*^3}.
    \]
    In genus two,
    \[
    \boxed{h_{j1}=Z_{j1}^*-\frac{\rho'(\mu_j)-y_*}{2\rho(\mu_j)},}
    \]
    \[
    \boxed{h_{j2}=Z_{j2}^*+\frac{(\partial_1 Z_{j1})_*}{y_*}
    +\frac{y_\rho(0)}{2\rho(\mu_j)}
    -\frac{\rho'(\mu_j)-y_*}{2\rho(\mu_j)^2}.}
    \]
    In genus three these first two formulas still apply, and
    \[
    \boxed{\begin{aligned}
    h_{j3}={}&Z_{j3}^*+\frac{3(\partial_2 Z_{j1})_*}{2y_*}
    +\frac{(\partial_1^2 Z_{j1})_*}{2y_*^2}
    -\frac{y_\rho(0)(\partial_1 Z_{j1})_*}{2y_*^2}\\
    &+\frac{y_{\rho\rho}(0)}{4\rho(\mu_j)}
    +\frac{y_\rho(0)}{2\rho(\mu_j)^2}
    -\frac{\rho'(\mu_j)-y_*}{2\rho(\mu_j)^3}.
    \end{aligned}}
    \]
    These reference-point expressions use sigma derivatives up to order three
    in genus two and four in genus three. The main construction now uses only
    second-order sigma derivatives at each fixed $w(\mu_j)$.
    `sigma_jets.py` keeps both evaluations and the independent finite Taylor
    calculation, so the original derivation remains available for comparison.

    Rerun the verification from the project root with
    `.venv/bin/python -m hyperelliptic_scalar_first.verification --dps 30`.
    All supporting code, results and plots are in this folder.
    """)
    return


if __name__ == "__main__":
    app.run()
