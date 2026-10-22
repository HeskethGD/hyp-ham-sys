# Weighted special Q and polynomial mode potentials

These are detailed proof and scope notes supporting the special-polynomial notebook.
The manuscript gives a shorter common-translation proof of the classification. The
historical notes are retained unchanged. Exact finite checks support the
arguments below; they do not replace the arbitrary-degree proofs.

## 1. Definitions and quantifiers

Let beta_j be fixed positive integers, M=sum beta_j, and Q a polynomial in
independent powers I_j. Set

\[
 I_j=\gamma_j-\beta_j\rho=\beta_j(a_j-\rho),\quad
 a_j=\gamma_j/\beta_j,\quad \sum_j\gamma_j=0,
\]
\[
 q(\rho)=Q(I(\rho)),\qquad r(\rho)=\prod_j I_j^{\beta_j},\qquad
 F=q^2-4r.
\]

Differentiate Q in its independent arguments before restricting it. Define

\[
 B_j(\rho)=\frac{q(\rho)-q(a_j)}{2(a_j-\rho)}-Q_j(I(\rho)).
\]

“Universally special” means B_j=0 for every j and every zero-sum tuple gamma,
as a polynomial identity in rho. The same Q must work for all tuples. This
does not require that every tuple produce a smooth curve. Smoothness and
nonbranch mode points are imposed only when selecting analytic examples.

For one fixed tuple there are extra solutions: adding L(I)^2 R(I), where L
vanishes on that power line, changes neither Q nor its first partials there.
No uniqueness of the ambient polynomial is asserted on a single leaf.

## 2. Complete weighted classification

**Theorem.** If M is odd, the universally special polynomials are Q=epsilon.
If M=2m is even, they are

\[
 Q=\kappa H_{\boldsymbol\beta}+\varepsilon,\qquad
 H_{\boldsymbol\beta}(I)
 =[t^m]\prod_i(1-tI_i/\beta_i)^{\beta_i/2}.
\]

Here kappa and epsilon are arbitrary complex constants. Coefficient extraction
is formal and involves no square-root branch. In particular H_beta is an
ordinary homogeneous polynomial of degree m.

**Proof.** Put x_i=I_i/beta_i and Qtilde(x)=Q(beta_i x_i). The condition is

\[
 \widetilde Q(x)-\widetilde Q(x-x_j\mathbf1)
 =\frac{2x_j}{\beta_j}\partial_j\widetilde Q(x). \tag{1}
\]

Indeed x=a-rho*1 and sum beta_i a_i=0. For arbitrary x choose
rho=-sum beta_i x_i/M and a=x+rho*1. Thus the universal condition is exactly
the unrestricted polynomial identity (1), including x_j=0 by continuation.

The identity preserves homogeneous degree. Let Qtilde_d have degree d>0.
On x=s*1 it gives 2s partial_j Qtilde_d=beta_j Qtilde_d. Summing and using
Euler's identity gives

\[
 (2d-M)\widetilde Q_d(s\mathbf1)=0. \tag{2}
\]

We need that a solution vanishing on the diagonal vanishes everywhere.
Use the unweighted coordinates s=sum x_i/N, y=x-s*1, sum y_i=0. If the
diagonal vanishes, write

\[
 \widetilde Q_d(s\mathbf1+y)=\sum_{r=1}^d s^{d-r}h_r(y).
\]

Take the least r with h_r nonzero. The left side of (1) has no term of
y-degree r-1: its first term starts at degree r and its shifted term
Qtilde_d(y-y_j*1) has degree d in y. On the right the degree r-1 term is

\[
 \frac{2}{\beta_j}s^{d-r+1}D_{e_j-\mathbf1/N}h_r.
\]

The directions e_j-1/N span the zero-sum hyperplane, so all their derivatives
vanishing makes h_r constant there. Its positive homogeneous degree forces
h_r=0, a contradiction. For N=1 the same diagonal assertion is immediate.

Consequently (2) forces d=M/2. Matching diagonal values of any two such
solutions and subtracting also proves uniqueness up to scale.

For existence, when M=2m form the formal Laurent series

\[
 R_x(u)=u^m\prod_i(1-x_i/u)^{\beta_i/2},\qquad
 P_x(u)=[R_x(u)]_{\ge0}.
\]

Its constant coefficient is Htilde(x). Translation gives
P_{x-a*1}(u)=P_x(u+a), hence Htilde(x-a*1)=P_x(a). Also

\[
 \partial_jP_x(u)=-\frac{\beta_j}{2}
      \frac{P_x(u)-P_x(x_j)}{u-x_j}.
\]

This follows by differentiating R_x and taking polynomial parts; the discarded
tail divided by u-x_j is O(u^-2). Evaluation at u=0 gives (1).
On the diagonal Htilde(s*1)=(-s)^m, so it is nonzero. Constants complete the
classification. QED.

For linear Q=epsilon+sum ell_j I_j, B_j=(sum beta_i ell_i)/2-ell_j.
These constant corrections generate extra linear exponentials. They are
excluded by B_j=0, except for the M=2 cases allowed by the theorem. This
explains why sigma-ratio solutions with constant phase-velocity shifts do not
contradict the classification.

## 3. Curve degree: total weight is not genus

Always deg r=M. Degree F depends also on q and cancellations. The original
paper's worked examples have (M,deg q,deg F)=(6,3,5) and (8,4,7); their
genera remain two and three and their Q are not special.

If Q=E is constant, F=E^2-4r has degree M. For odd M=2g+1 it is an
odd-degree genus-g model provided F is squarefree. E is an energy-level
constant: H=U+V-E=0 is the same vector field as U+V on energy E.

For a nonconstant special Q with M=2m, put

\[
 C=\prod_i\beta_i^{\beta_i},\quad
 S(u)=\prod_i(u-a_i)^{\beta_i},\quad
 P(u)=[\sqrt{S(u)}]_{\ge0},\quad P(u)\sim u^m.
\]

Translation in the preceding proof gives q(u)=kappa P(u)+epsilon and
r(u)=C S(u). Since P^2-S has degree at most m-1,

\[
 F=(\kappa^2-4C)P^2+4C(P^2-S)+2\kappa\varepsilon P+\varepsilon^2.
\]

Thus if kappa^2 != 4C, deg F=2m. If kappa^2=4C, deg F<=m;
for epsilon nonzero it is exactly m, and for epsilon=0 it is at most m-1.
The latter may be zero or singular for some weight patterns. Smoothness,
odd degree and nonzero fixed ordinates must be checked, not presumed.

The (5,5) example in model.py has kappa=6250 and epsilon=1. It gives a
smooth quintic and nonbranch fixed points, with only two canonical pairs.
The total weight is ten; it is not a genus-four example because of cancellation.

The new examples are kept in their original Hamiltonian scaling. Before
evaluating the paper's leading-four sigma formulas, their scalar scaling and
the associated Abel differentials must be translated consistently.

## 4. Wronskian and general second derivatives

Use W_j=u_j v_j'-v_j u_j'. Put d=a_j-rho, b=q(a_j), B=B_j, y=rho'. Then

\[
 \frac{u_j'}{u_j}=\frac{b-y}{2d}+B,\qquad
 \frac{v_j'}{v_j}=\frac{-b-y}{2d}-B,\qquad u_jv_j=\beta_jd.
\]

Therefore

\[
 W_j=-\beta_j[b+2dB]. \tag{3}
\]

On a nonconstant trajectory a polynomial B gives constant W_j iff dB is
constant as a polynomial. Evaluation at rho=a_j makes that constant zero,
so B=0. Hence specialness is equivalent to constant Wronskians on each
trajectory; universal specialness requires this on every leaf.

Let

\[
 G_b=\frac14\frac{d}{d\rho}\frac{F-b^2}{\rho-a_j}.
\]

Using y^2=F and rho''=F'/2 gives

\[
 \frac{u_j''}{u_j}=G_b+B^2+\frac{bB}{d}
                  +y\left(B'-\frac B d\right),
\]
\[
 \frac{v_j''}{v_j}=G_b+B^2+\frac{bB}{d}
                  -y\left(B'-\frac B d\right). \tag{4}
\]

Here primes on B,F mean rho derivatives. These are individual equations;
general modes do not already share a Schrödinger potential.

For nonsquare polynomial F, 1 and y are independent over C(rho). Thus an
individual expression in (4) being polynomial in rho forces dB'-B=0.
With B polynomial this forces B=0. Conversely B=0 gives a polynomial
potential because b^2=F(a_j). This strengthens the connection within the
stated polynomial Hamiltonian class, without classifying arbitrary functions
or arbitrary changes of gauge.

## 5. Prescribed-product uniqueness and polynomial potentials

Suppose more generally uv=beta(a-rho), W=beta*c is constant, and
rho'^2=F(rho), rho''=F'/2. Then

\[
 \frac{u'}u=\frac{-c-\rho'}{2(a-\rho)},\qquad
 \frac{v'}v=\frac{c-\rho'}{2(a-\rho)},
\]
\[
 u''=G_{a,c}(\rho)u,\quad v''=G_{a,c}(\rho)v,\qquad
 G_{a,c}=\frac14\frac{d}{d\rho}\frac{F(\rho)-c^2}{\rho-a}. \tag{5}
\]

For polynomial F, divide F-c^2 by rho-a. The remainder is F(a)-c^2, so

\[
 G_{a,c}=\frac14\frac{d}{d\rho}
       \frac{F(\rho)-F(a)}{\rho-a}
       +\frac{c^2-F(a)}{4(\rho-a)^2}.
\]

It is polynomial iff c^2=F(a); for deg F=D>=2 the potential degree is D-2.
This is uniqueness of the squared **normalised** Wronskian c^2 under the
prescribed-product constraint, not uniqueness among all ODE solutions.
Interchanging u,v reverses c; reciprocal rescaling preserves product and W.

For special canonical modes c=-q(a). They are independent provided q(a)!=0.
For constant Q=E!=0, W_j=-beta_j E and the fixed points are distinct
(a_j,E), not coincident points. For E=0 they are branch points: the curve
may stay smooth but each pair becomes proportional to sqrt(a_j-rho).
The historical notebook's polynomial criterion must refer to f^2=F being
polynomial, not f being polynomial.

Substituting rho=a-uv/beta in (5) makes the potential polynomial in the paired
power, with conserved data as coefficients. For F=4rho^3-g2 rho-g3,
G=2rho+a. For higher degree, a generally enters nonconstant coefficients;
there is no automatic common Schrödinger operator with a as spectral energy.
“Lamé-type” is an analogy supported by this exact genus-one recovery.
Baker–Akhiezer identification and finite-gap spectral claims remain open.

## 6. Elementary solutions and determinant scope

For B_j=0 on a simply connected regular chart,

\[
 \frac{u_j(z)}{u_j(z_0)}=
 \sqrt{\frac{a_j-\rho(z)}{a_j-\rho(z_0)}}
 \exp\left[\int_{\rho(z_0)}^{\rho(z)}
 \frac{q(a_j)\,d\xi}{2(a_j-\xi)y(\xi)}\right],
\]

and v has the opposite exponent with its own initial amplitude. Branches and
paths are continued compatibly. Initial products must satisfy uv=I and
prod u_j^beta_j=(q+y)/2. The paper's corrected sigma point factors integrate
these same differentials after translating its curve normalization.

For constant Q=E, the paper's n=g+2 determinant argument applies directly:
D/Delta=(y-E)/2=-V, D*/Delta=(y+E)/2=U. Thus the Hamiltonian identity is
(D*-D)/Delta=E; an individual determinant is not the Hamiltonian itself.
The raw-sigma product identification additionally uses the paper's compatible
Abel lifts and factorization hypotheses (including distinct lists); weights
requiring repeated rows need the stated confluent extension. Do not infer
all weighted sigma determinant factorizations from the classification alone.

Nonconstant special odd-degree examples can have deg q>g+1. Their individual
mode differential still satisfies B=0, but the paper's short determinant
column construction does not cover them. A longer-column argument is separate
work, not an implicit consequence of the new classification.

## 7. Verification coverage and next work

- Exact coefficient searches independently solve finite homogeneous systems.
- Symbolic checks verify weighted identities, second derivatives, pole
  obstructions, genus-one recovery and degree cancellation.
- Original example models are imported read-only to confirm their original
  genera and nonzero corrections.
- Three smooth odd-degree examples are checked by Abel quadrature versus full
  canonical RK4, with step refinement, and a separate second-order RK4 solve.
- These short numerical paths do not evaluate sigma functions, certify global
  continuation, or establish new determinant identities.

The following stage in `../sigma_modes/` now normalises these curves for
Generapy and checks corrected sigma factors against these quadratures,
including first and second sigma derivatives and Wronskians. Its fixed-point
first- and second-kind integrals follow prescribed compatible paths; see its
normalisation record for the resolved path/quadrature comparison issue.

Remaining: develop the constant-Q sigma determinant specialization under
explicit multiplicity hypotheses, longer paths and branch-point limits.
The nonconstant example still needs the longer-column determinant argument.
