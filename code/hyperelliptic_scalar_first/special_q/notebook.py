import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import sys as _sys
    from pathlib import Path
    notebook_dir = Path(__file__).resolve().parent
    if str(notebook_dir.parent.parent) not in _sys.path:
        _sys.path.insert(0, str(notebook_dir.parent.parent))
    import importlib.metadata as metadata
    import json
    import marimo as mo
    import sympy as sp
    from hyperelliptic_scalar_first.special_q.model import example_models
    from hyperelliptic_scalar_first.special_q.verification import exact_checks, numerical_checks
    return exact_checks, example_models, json, metadata, mo, notebook_dir, numerical_checks, sp


@app.cell
def _(mo):
    mo.md(r"""
    # Special polynomial Hamiltonians and their polynomial mode potentials

    The paper separates the mode equation using the correction:
    \[
    B_j(\rho)=\frac{q(\rho)-q(a_j)}{2(a_j-\rho)}-Q_j(I(\rho)),
    \qquad I_j=\beta_j(a_j-\rho),\quad a_j=\gamma_j/\beta_j.
    \]
    We call $Q$ **universally special** when every $B_j$ vanishes for every
    zero-sum invariant tuple. Differentiate the polynomial in independent
    powers **before** restricting to the invariant line. The elementary
    sigma factors retain their natural second-kind exponentials; it is only
    the extra $B_j$ correction that disappears.

    This notebook executes exact identities, independent finite coefficient
    searches, and short numerical comparisons. The full arbitrary-degree
    proofs and scope are recorded in the adjacent `proof_scope.md`.
    Numerical checks here use Abel quadrature, not sigma evaluation.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 1. The weighted classification

    Write $M=\sum_j\beta_j$. The universally special polynomials are
    \[
    Q=\begin{cases}
    \varepsilon,&M\text{ odd},\\
    \kappa g_m^{(-\boldsymbol\beta/2)}(I_1/\beta_1,\ldots,I_N/\beta_N)+\varepsilon,&M=2m,
    \end{cases}\qquad
    \prod_{i=1}^N(1-tx_i)^{\beta_i/2}
       =\sum_{n=0}^{\infty}g_n^{(-\boldsymbol\beta/2)}(x)t^n.
    \]
    These are the Chan–Chyan–Srivastava multivariable Lagrange polynomials,
    with parameters $-\boldsymbol\beta/2$, in the notation used by the paper.
    Here $\kappa,\varepsilon$ are arbitrary constants. The expression is
    formal coefficient extraction, so $H$ is an ordinary polynomial.

    **Proof structure.** In normalised variables $x_i=I_i/\beta_i$, the
    cancellation condition is the unrestricted identity
    \[
    \widetilde Q(x)-\widetilde Q(x-x_j\mathbf1)
       =\frac{2x_j}{\beta_j}\partial_j\widetilde Q(x).
    \]
    On the diagonal, a homogeneous degree-$d$ solution satisfies
    $(2d-M)\widetilde Q_d(s\mathbf1)=0$. A solution vanishing on the
    diagonal must vanish everywhere: expansion transverse to the diagonal
    makes the first nonzero transverse term have zero derivatives in every
    tangent direction, a contradiction. Thus only $d=M/2$ can occur, and
    its solution space has dimension at most one. The polynomial part of
    $\prod_i(u-x_i)^{\beta_i/2}$ constructs that solution; its constant
    coefficient is $H$, and $H(\beta_1s,\ldots,\beta_Ns)=(-s)^m$.

    The searches below solve polynomial coefficient equations without using
    the proposed formula. They test degrees 0–4 for twelve weight patterns;
    these finite searches supplement the proof, not replace it.
    """)
    return


@app.cell
def _(exact_checks):
    exact_report = exact_checks()
    return (exact_report,)


@app.cell
def _(exact_report, mo):
    mo.vstack([
        mo.md(f"**{exact_report['check_count']} exact checks passed.**"),
        mo.ui.table(exact_report["coefficient_search"], selection=None, page_size=15),
    ])
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 2. Constant phase shifts, invariant leaves, and curve degree

    Linear $Q=\varepsilon+\sum_j\ell_j I_j$ gives
    $B_j=\tfrac12\sum_i\beta_i\ell_i-\ell_j$, a constant.
    These produce extra $e^{B_jz}$ factors, explaining why the three-mode
    sigma solutions with constant phase-velocity shifts do not contradict
    the theorem. Requiring $B_j=0$ is a stronger condition.

    On one fixed invariant leaf, adding $L(I)^2R(I)$ with $L=0$ on that
    leaf changes neither $q$ nor the restricted first partials. The theorem
    requires the same $Q$ to work on **every** leaf.

    For constant $Q=E\ne0$, $F=E^2-4r$ has degree $M$. Thus odd total
    weight gives the initial odd-degree family. For even $M=2m$, put
    $C=\prod_i\beta_i^{\beta_i}$, $S(u)=\prod_i(u-a_i)^{\beta_i}$,
    and $P=[\sqrt S]_{\ge0}$. A special $Q$ restricts to
    $q=\kappa P+\varepsilon$ with $\deg(P^2-S)\le m-1$.
    Therefore:

    * $\kappa^2\ne4C$: degree $F=2m$;
    * $\kappa^2=4C$, $\varepsilon\ne0$: degree $F=m$;
    * $\kappa^2=4C$, $\varepsilon=0$: degree $F\le m-1$.

    Odd-degree nonconstant examples are possible after cancellation. Smoothness
    and nonzero fixed ordinates must still be checked. These operations change
    the curve; they do not revise the genus of either original example.

    For the constant-$Q$ examples we also check the actual lower-degree
    determinant columns $(1,\rho,\ldots,\rho^g,y/2)$:
    $D/\Delta=(y-E)/2=-V$ and $D^*/\Delta=(y+E)/2=U$.
    Their difference gives $(D^*-D)/\Delta=E$. These are exact polynomial
    determinant checks; the sigma factorization is a separate analytic step.
    """)
    return


@app.cell
def _(exact_report, mo):
    mo.vstack([
        mo.md("**Degree cancellation checks**"),
        mo.ui.table(exact_report["degree_cases"], selection=None),
        mo.md("**General mode examples, imported read-only**"),
        mo.ui.table(exact_report["original_examples"], selection=None),
    ])
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 3. Wronskians and general second derivatives

    Set $d=a_j-\rho$, $b=q(a_j)$, $B=B_j$, and $y=\rho'$. The separated
    equations give
    \[
    W_j=u_jv_j'-v_ju_j'=-\beta_j[b+2dB].
    \]
    Since $B$ is polynomial, a nonconstant trajectory has constant $W_j$
    **if and only if** $B_j=0$. Its value is then $-\beta_jq(a_j)$;
    independence additionally requires $q(a_j)\ne0$.

    Define $G_b=\tfrac14\partial_\rho[(F-b^2)/(\rho-a_j)]$. Direct
    differentiation gives
    \[
    \frac{u_j''}{u_j}=G_b+B^2+\frac{bB}{d}
                       +y\left(B'-\frac Bd\right),\qquad
    \frac{v_j''}{v_j}=G_b+B^2+\frac{bB}{d}
                       -y\left(B'-\frac Bd\right).
    \]
    For nonsquare $F$, polynomiality in $\rho$ alone forces $dB'-B=0$,
    hence $B=0$ when $B$ is polynomial. Special modes thus have a common
    polynomial potential for each pair. Different pairs need not have the
    same potential.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 4. The precise uniqueness statement

    Prescribe $uv=\beta(a-\rho)$ and a constant Wronskian $W=\beta c$.
    Both functions solve
    \[
    \psi''=G_{a,c}(\rho)\psi,\qquad
    G_{a,c}=\frac14\partial_\rho\frac{F(\rho)-c^2}{\rho-a}.
    \]
    Polynomial division isolates the only obstruction:
    \[
    G_{a,c}=\frac14\partial_\rho
      \frac{F(\rho)-F(a)}{\rho-a}
       +\frac{c^2-F(a)}{4(\rho-a)^2}.
    \]
    Thus the potential is polynomial **iff $c^2=F(a)$**. Its degree is
    $\deg F-2$. This uniquely determines the squared normalised Wronskian,
    subject to the prescribed product. It does not assert uniqueness among
    arbitrary ODE solutions; sign exchange and reciprocal rescaling remain.

    Substitution $\rho=a-uv/\beta$ makes the potential a polynomial in
    the paired power. For the Weierstrass cubic it becomes $2\rho+a$,
    recovering the genus-one Lamé equation. A common higher-genus spectral
    operator or Baker–Akhiezer identification is not established here.

    When constant $Q=E=0$, the fixed points are branch points and $W=0$:
    each pair is proportional to $\sqrt{a_j-\rho}$. The curve can remain
    smooth, but the independent-pair statement no longer applies.
    """)
    return


@app.cell
def _(example_models):
    selected_models = example_models()
    return (selected_models,)


@app.cell
def _(mo, selected_models, sp):
    _items = [mo.md("## 5. Three smooth odd-degree examples\n\nThese retain their original Hamiltonian scaling; leading-four sigma normalization is a later step.")]
    for _data in selected_models:
        _items.append(mo.md(
            f"### `{_data['name']}`\n\n"
            f"Weights `{_data['beta']}`, invariants `{_data['gamma']}`. "
            f"Scalar degree **{_data['degree']}**, genus **{_data['genus']}**, "
            f"squarefree **{_data['squarefree']}**.\n\n"
            + r"\[q(\rho)=" + sp.latex(_data['q']) + r"\]"
            + r"\[F(\rho)=" + sp.latex(_data['F']) + r"\]"
        ))
    mo.vstack(_items)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 6. Numerical checks against the original canonical field

    On a regular positive-sheet path the elementary modes satisfy
    \[
    \frac{u_j}{u_{j0}}=
    \sqrt{\frac{a_j-\rho}{a_j-\rho_0}}
    \exp\int_{\rho_0}^{\rho}
       \frac{q(a_j)\,d\xi}{2(a_j-\xi)\sqrt{F(\xi)}}.
    \]
    The $v_j$ expression has the opposite exponent and its own initial
    amplitude. We evaluate these integrals at 50-digit precision, including
    $z-z_0=\int d\xi/\sqrt F$. Independently, we differentiate the full
    Hamiltonian in canonical variables and integrate that polynomial field
    with RK4 at 16, 32 and 64 steps. Refinement should reduce errors by
    approximately sixteen. A separate second-order IVP checks the derived
    potential equations and their constant Wronskians.

    These short paths test the reconstruction and differential equations.
    They are not sigma evaluations or global continuation tests.
    """)
    return


@app.cell
def _(numerical_checks):
    numerical_report = numerical_checks(dps=50)
    return (numerical_report,)


@app.cell
def _(mo, numerical_report):
    mo.ui.table(numerical_report["examples"], selection=None)
    return


@app.cell
def _(exact_report, json, metadata, notebook_dir, numerical_report):
    saved_report = dict(
        environment={name: metadata.version(name) for name in ("sympy", "mpmath", "marimo")},
        exact=exact_report, numerical=numerical_report,
    )
    report_path = notebook_dir / "results.json"
    report_path.write_text(json.dumps(saved_report, indent=2) + "\n")
    return (report_path,)


@app.cell
def _(mo, report_path):
    mo.md(r"""
    ## 7. Scope and the next numerical stage

    The constant-$Q$ family fits the paper's short determinant columns,
    but its sigma factorization retains the distinct-list and compatible-lift
    hypotheses. Repeated rows require confluent identities. The nonconstant
    genus-two example has $\deg q=5>g+1$ and needs a longer-column argument;
    its individual-mode differential is already covered by $B_j=0$.

    Next: consistently normalise the selected curves, compute periods and
    fixed-point lifts, evaluate the paper's corrected sigma factors, and
    compare them against these quadratures. Keep the branch-point limit,
    global monodromy, and any finite-gap spectral claims separate.
    """ + f"\n\nExecuted results saved to `{report_path.name}`.")
    return


if __name__ == "__main__":
    app.run()
