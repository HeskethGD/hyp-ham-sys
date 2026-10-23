# Simple branch fixed points in arbitrary genus

Working proof for possible integration into the main paper. All claims below
retain a smooth odd-degree hyperelliptic curve, the paper's sigma convention,
and compatible unreduced first- and second-kind lifts. This note changes no
main-paper statement. The numerical checks are corroboration of the proof,
not its justification in arbitrary genus.

## 1. Setting and sources

Let C be the smooth compactification of y²=F(x), with squarefree F of degree
2g+1 and leading coefficient 4, g>=2. Write A(P) for the infinity-based Abel
lift with omega_k=x^(k-1) dx/y, and set

\[
 h_k(P)=-\operatorname{Reg}\int_\infty^P dr_k.
\]

The regularization is only at infinity. At a finite simple branch point
P_a=(a,0), x-a=t² and y=t s(t²), with s(0)²=F'(a)!=0. Consequently dx/y
and the paper's polynomial multiples defining dr_k are holomorphic in t.
Both A(P) and h(P) therefore extend holomorphically to P_a on each compatible
local lift. In particular, no divergent finite endpoint subtraction is needed.

The exact sigma facts used are:

- [Ōnishi, Proposition 6.6(1)](https://arxiv.org/pdf/math/0105189): on the
  lifted one-point Abel image, sigma_sharp has zeros only over infinity.
- [Ōnishi, Lemma 7.1 and Theorem 7.2](https://arxiv.org/pdf/math/0105189):
  the two-point and determinant identities on that image. Use the author's
  [published determinant-sign correction](https://math-meijo-u.heteml.net/research/pub/correct_pem.pdf),
  as already specified in the main paper.
- [Gibbons–Matsutani–Ōnishi, Theorem 7.17](https://arxiv.org/pdf/1204.3747):
  the hyperelliptic prime-form sigma representation used in the main kernel
  lemma. A simple branch point is an ordinary smooth point of C for this
  identity; x itself is simply not a local coordinate there.

The arguments below are deductions from these identities and the main
paper's nonbranch kernel lemma. They are not additional claims attributed
to those sources.

For reproducibility, the full determinant checks use
\[
 D_n=c_n\frac{\sigma(\sum_{i=1}^n u_i)
                   \prod_{1\le i<j\le n}\sigma_\flat(u_i-u_j)}
                  {\prod_{i=1}^n\sigma_\sharp(u_i)^n},\qquad
 c_n=(-1)^{(n-g-1)(n+g^2+2g)/2},\quad n\ge g,
\]
with the first n pole-ordered column functions and monic ordinate y/2.
The cited errata calls this determinant result Proposition 8.2(ii) in the
published numbering; the arXiv v8 numbering is Theorem 7.2(2).

## 2. Proposition: the mode kernel extends to every finite simple branch point

For P near P_a put eta=A(P), x(P)=a(P), y(P)=b(P), and define

\[
 R_+(w;P)=\frac{\sigma_\flat(w+\eta)}
 {\sigma_\sharp(w)\sigma_\sharp(\eta)},\qquad
 R_-(w;P)=\frac{\sigma_\flat(w-\eta)}
 {\sigma_\sharp(w)\sigma_\sharp(\eta)}.
\]

For a moving point Q distinct from P_a and infinity, and its compatible
lift w=A(Q), both quotients and h(P) extend holomorphically to P=P_a.
Their product and logarithmic differential have the same formulas as in
the main paper, now with b=0:

\[
 R_+R_-=(-1)^{g+1}(a-x),
\]
\[
 \frac{d}{dz}\log R_+
 =\frac{y}{2(x-a)}+\sum_{k=1}^g h_k(P_a)x^{k-1},\qquad
 \frac{d}{dz}\log R_-
 =\frac{y}{2(x-a)}-\sum_{k=1}^g h_k(P_a)x^{k-1}.
\]

**Proof.** Proposition 6.6(1) gives sigma_sharp(eta_a)!=0, since P_a is
finite; the moving denominator is also nonzero for finite Q. Sigma
derivatives are entire, and A(P), h(P) are holomorphic in the local
parameter at P_a. Thus R_+, R_- and their moving derivatives converge as
P tends to P_a through nonbranch points. The two-point identity gives
a nonzero product when x!=a, so neither numerator vanishes on this
regular chart. Pass to the limit in the kernel lemma and in the two-point
identity. This proves the displayed identities without taking a quotient
of two vanishing fixed-point normalizers. The same continuation argument
justifies a prime-form derivation directly at the branch point. □

### Corollary: unchanged shifted sigma representation

At a branch fixed point, the normalized canonical mode formulas remain

\[
 u_j=u_{j0}\frac{R_+(w)}{R_+(w_0)}
       e^{-h_j\cdot(w-w_0)+\int_{z_0}^{z}B_j\,dz},\qquad
 v_j=v_{j0}\frac{R_-(w)}{R_-(w_0)}
       e^{+h_j\cdot(w-w_0)-\int_{z_0}^{z}B_j\,dz}.
\]

The upper and lower quotients are independently evaluated, just as in the
nonbranch case. B_j is unchanged and remains a polynomial: its divided
difference has a removable denominator regardless of b_j. The existing
degree bound still makes its primitive a holomorphic Abel combination.
Neither h_j nor the natural exponential should be set to zero at a branch
point. The special-Q classification is algebraic and needs no alteration.

For B_j=0, the two normalized factors coincide locally, because their
logarithmic derivatives are both y/[2(x-a)]. Their square is
(a-x)/(a-x_0), with initial value one. This is an identity satisfied by
the shifted sigma formulas, not a replacement for them.

### Continuation through the attached branch point

The differential dz=dx/y=2 dt/s(t²) is nonzero at t=0, so z is a valid
local analytic coordinate there. The corrected logarithmic differential is
(1/2) dx/(x-a)=dt/t. Each corrected elementary mode is therefore t times
a nonzero local constant. For general B, its additional exponential is
holomorphic and nonvanishing in z. Thus the normalized formulas have an
analytic continuation through their simple zero. A logarithmic equation
is used off the zero; the polynomial canonical equations hold at the zero
by analytic continuation. No extra sigma derivative is required.

### Smoothness and the second-order equation

On the Hamiltonian power line, F=q²-4r and r(a_j)=0. If b_j=q(a_j)=0,

\[
 F'(a_j)=-4r'(a_j).
\]

For distinct fixed coordinates and positive integer weights, r has a zero
of order beta_j at a_j. Smoothness therefore forces beta_j=1 at every
branch fixed point. This does not forbid repeated weights at other,
nonbranch fixed points.

For B_j=0, the pair's Wronskian is zero and its modes are proportional,
but both satisfy the same polynomial equation

\[
 \psi''=G_j(x)\psi,\qquad
 G_j(x)=\frac14\frac{d}{dx}\frac{F(x)}{x-a_j}.
\]

F(a_j)=0 makes the quotient polynomial. Its degree after differentiation
is 2g-1. Independence still requires b_j!=0; that hypothesis belongs only
to the assertion about a fundamental pair, not to the sigma mode formula
or polynomiality assertion. No second independent solution is claimed here.

## 3. Proposition: the two determinant families admit branch fixed rows

Assume the same interpolation families as the main paper:

- n=g+2 with degree(q)<=g, columns 1,x,...,x^g,y/2;
- n=g+3 with degree(q)=g+1, columns 1,x,...,x^g,y/2,x^(g+1).

The fixed points P_j=(a_j,q(a_j)) have distinct coordinates. Some may be
simple branch points. Let A select n-1 fixed rows, define
Delta_A=product_{i<j in A}(a_j-a_i), and put epsilon=1 or -1 for the two
families respectively. Then

\[
 D_n=\epsilon\Delta_A(y-q)/2,\qquad
 D_n^*=\epsilon\Delta_A(y+q)/2.
\]

**Proof.** Subtract q/2 from the ordinate column. Every fixed entry
vanishes, including a branch entry where both ordinate and q are zero.
Expand the remaining moving entry; its cofactor is the same nonzero
Vandermonde. Negate the fixed ordinates to obtain the second formula.
At a branch row this negation has no effect on the matrix. □

### Correct divisor count at a branch point

At nonbranch fixed points, the existing proof uses q+y!=0 to conclude
ord_{P_j}(q-y)=beta_j. At a branch fixed point that argument must be
replaced. Here beta_j=1 and, with x-a_j=t²,

\[
 q(x)=O(t^2),\quad y=s(0)t+O(t^3),\quad
 q-y=-s(0)t+O(t^2),\quad q+y=s(0)t+O(t^2).
\]

Thus **both factors have a simple zero**, and the branch point occurs
once in each zero divisor, not twice in either. There are no other finite
zeros, because (q-y)(q+y)=4r. At infinity the pole order of q-y is M:
2g+1 in the first family and 2g+2 in the second. Consequently

\[
 \operatorname{div}(q-y)=\sum_{j=1}^N\beta_j P_j-M\infty
\]

still holds. Removing the n-1 selected occurrences still leaves g
residual occurrences. No new row multiplicity is caused by branching.

### Sigma identity and residual splitting

Require the same combinatorial condition as in the current paper: A and
the residual list B have distinct coordinates internally and together
account for the weighted zero divisor. A unit-weight branch point belongs
to exactly one of these lists. Ōnishi's identity applies at branch rows
directly: its fixed sigma_sharp factors are nonzero. The fixed-fixed
sigma_flat difference factors are also nonzero. Indeed their products
with the corresponding sum factors are nonzero by the two-point identity
and the distinct coordinates.

Abel's theorem gives sum_j beta_j eta_j in the period lattice. A branch
point has unit weight, so this period can be absorbed into one fixed
point's path. On the resulting compatible lifts set

\[
 \sum_{j=1}^N\beta_j\eta_j=0,\qquad
 E=\sum_{j\in A}\eta_j=-\sum_{j\in B}\eta_j.
\]

Use **-eta_j** for the residual opposite-sheet lifts, even when P_j is
a branch point and its geometric involution is itself. At a branch point
2 eta_j is a period; identifying eta_j with -eta_j as vectors would lose
the sigma quasi-period factors.

The main identity groups into D_n=K_A [sigma(w+E)/sigma_sharp(w)]
product_{j in A} R_-(w;P_j), with K_A nonzero. The auxiliary identity
uses the g opposite-lift residual points and moving point, with polynomial
columns 1,x,...,x^g. Its determinant is
Delta_B product_{j in B}(x-a_j). Its sigma expression is
K_B [sigma(w+E)/sigma_sharp(w)] product_{j in B} R_+(w;P_j).
The two-point identity then gives, on a regular chart,

\[
 \frac{\sigma(w+E)}{\sigma_\sharp(w)}
 =\Gamma_B\prod_{j\in B}R_-(w;P_j),\qquad \Gamma_B\ne0.
\]

No cancellation is made at a zero; equality extends there by analyticity.
The complete raw products are therefore constant multiples of q+y and
q-y. The existing kernel differential argument now proves

\[
 \sum_{j=1}^N\beta_j(c_{jk}-h_{jk})=0,
 \qquad k=1,\ldots,g,
\]

on these balanced lifts, exactly as before. Hence the full sigma modes
give U=epsilon D_n^*/Delta_A and V=-epsilon D_n/Delta_A.

### The all-branch zero-energy case

For constant Q=0 and squarefree F, all weights are one and N=2g+1.
Both branch zero divisors coincide, q=0, and D_n^*=D_n. Therefore
U=y/2 and V=-y/2 are perfectly compatible with the determinant identity;
it no longer distinguishes two different fixed-sheet matrices. The full
sigma factorization and exponential cancellation still hold after lift
balancing. Before balancing, natural exponential corrections need not
cancel visibly from the raw sigma product.

## 4. What this establishes and what it does not

The branch-point exclusion can be removed from the local sigma mode theorem
on smooth odd-degree curves, with unit weight forced at branch fixed points.
The two stated determinant families also extend, under their existing
distinct-list and column-degree assumptions. The proof change is the local
divisor count and careful retention of opposite Abel lifts at half-periods.

This does not extend to singular curves, repeated rows inside a determinant,
even-degree curves, longer Onishi column spaces, or a claim that the two
branch-attached canonical modes are independent. It assumes the main
paper's general kernel lemma; it supplies the branch-point extension of
that lemma, rather than re-proving its nonbranch normalization from scratch.
