# Higher-degree corrections and sigma partitions

These checks support the main paper's polynomial correction integrals and
longer determinant column lists on smooth odd-degree curves. Write
\[
y^2=F(x)=q(x)^2-4r(x),\qquad \deg F=2g+1,
\qquad F(x)=4x^{2g+1}+\cdots.
\]
The leading coefficient 4 is the normalization used in the computations;
other nonzero leading coefficients give the same reduction after scaling.

## 1. Arbitrary polynomial corrections

There is no degree obstruction to integrating a polynomial correction
\(B(x)\,dx/y\). In the paper's first- and second-kind bases it has the unique
decomposition
\[
\frac{B(x)}y\,dx=d(P(x)y)+\sum_{k=1}^g c_k\,dw_k
                         +\sum_{k=1}^g d_k\,dr_k,
\qquad dw_k=\frac{x^{k-1}}y\,dx.
\]
Here the coefficients are constant on the chosen power line, and
\[
dr_k=\frac{D_k(x)}y\,dx,\qquad
D_k(x)=\frac14\sum_{m=k}^{2g+1-k}(m+1-k)\lambda_{m+1+k}x^m,
\quad F(x)=\sum_{\ell=0}^{2g+1}\lambda_\ell x^\ell.
\]

**Proof.** The numerator of \(d(x^k y)\) is
\(kx^{k-1}F+x^kF'/2\), with degree \(k+2g\) and leading coefficient
\(4k+4g+2\). Successive subtraction therefore removes every term of degree
at least \(2g\). Next, \(D_k\) has degree \(2g-k\) and leading coefficient
\(2g+1-2k\). These remove the remaining degrees \(2g-1,\ldots,g\)
successively. What remains is the holomorphic part. The same distinct
leading degrees prove uniqueness.

Consequently the extra exponent of mode \(j\) is
\[
\int_{z_0}^z B_j(x(s))\,ds
 =[P_j(x)y]_{z_0}^z+\sum_{k=1}^g c_{jk}\Delta w_k
                           +\sum_{k=1}^g d_{jk}\Delta R_k,
\qquad dR_k=dr_k,
\]
because \(x'=y\). This multiplies the existing corrected point factors
by opposite exponentials in \(u_j\) and \(v_j\). The primitives must use
compatible lifts. No additional third-kind integral is needed: these
polynomial differentials have only a pole at infinity and zero residue.

There is also a useful cancellation in complete products. Put
\(a_j=\gamma_j/\beta_j\), \(b_j=q(a_j)\), and
\[
\alpha_j=\frac{b_j-y}{2(a_j-x)}\frac{dx}{y}.
\]
The scalar equations imply
\[
d\log(q+y)-\sum_j\beta_j\alpha_j
             =\sum_j\beta_j B_j(x)\frac{dx}{y}.
\]
The terms on the left have matching finite residues. At a nonbranch
fixed point, \(q+y\) has its zero of order \(\beta_j\) on the opposite
sheet; at a smooth branch fixed point \(\beta_j=1\). At infinity the
residue of each complete expression is \(-M\), where
\(M=\sum_j\beta_j\). All these poles are at most simple, so their
difference is holomorphic. Thus
\[
\deg\!\left(\sum_j\beta_j B_j\right)\leq g-1.
\]
Uniqueness of the reduction gives
\(\sum_j\beta_j P_j=0\) and \(\sum_j\beta_j d_{jk}=0\).
The new endpoint and second-kind terms cancel in the complete products,
although they are needed in the individual modes. This argument uses the
usual distinct fixed coordinates and smooth curve assumptions.

## 2. Longer determinants

Let \(d=\deg q\geq g+1\), with leading coefficient \(q_d\). The initial
pole-ordered basis ending at \(x^d\) is formed from
\[
1,x,\ldots,x^d,\qquad y/2,xy/2,\ldots,x^{d-g-1}y/2,
\]
interleaved by pole order. Its length is \(n=2d-g+1\). Choose \(n-1\)
distinct fixed coordinates \(a_i\), with ordinates \(q(a_i)\), and put
\(\Delta_A=\prod_{1\leq i<j\leq n-1}(a_j-a_i)\).
For the determinant with these fixed rows followed by the moving row,
\[
D_n=C(q_d)\Delta_A\frac{y-q(x)}2,
\qquad
C(q_d)=(-1)^{1+m(m+1)/2}(q_d/2)^m,\quad m=d-g-1.
\]
Reversing every fixed sheet gives
\[
D_n^*=C(-q_d)\Delta_A\frac{y+q(x)}2.
\]
In particular the constants need not be the same:
\(C(-q_d)=(-1)^m C(q_d)\).

**Proof.** Subtract \(q/2\) from the \(y/2\) column. That column now
vanishes on every fixed row. Expanding along it leaves, on those rows,
the polynomial basis
\(1,x,\ldots,x^d,xq/2,\ldots,x^m q/2\).
Its successive degrees are \(0,\ldots,n-2\); hence its determinant is
\((q_d/2)^m\Delta_A\). Restoring the pole ordering gives the displayed
sign. Replacing \(q\) by \(-q\) proves the opposite-sheet formula.

For \(d>g+1\), maintaining the odd degree \(2g+1\) forces
\(\deg r=2d\), so \(M=2d\). Removing the \(n-1=2d-g\) fixed zeros
therefore leaves total multiplicity \(g\), exactly as required for the
existing residual sigma partition. When the zeros can be partitioned
into internally distinct lists, the longer Frobenius--Stickelberger
identity and the same \((g+1)\)-row residual identity give the individual
sigma factors. Period balancing and initial normalization proceed as
before. The new \(x^k y\) columns do not remain as extra factors in the
final product: the graph condition \(y=q(x)\) removes them as above.

Repeated coordinates inside one determinant would require confluent
identities; that extension is not checked here. The numerical examples
have unit weights and distinct fixed coordinates. This notebook does
not establish formulas for even-degree or singular curves.

## 3. What the computations test

The genus-two example has eight unit-weight modes,
\[
\gamma=(-7,-5,-3,-1,1,3,5,7),\quad
q=2x^4-84x^2+x+210,
\quad F=4x^5-168x^3+16385x^2+420x.
\]
With \(S=-\sum I_j/8\) and \(L=I_1-I_2+2\), choose the actual
multivariate Hamiltonian polynomial \(Q=q(S)+LS^5\). On the power line
\(L=0\), so the restricted \(q\) is unchanged while its partial
derivatives produce degree-five corrections. This tests both the
seven-row determinant containing \(xy/2\) and the new differential
reduction in one system.

The genus-three example uses seven unit-weight modes,
\(\gamma=(-3,-2,-1,0,1,2,3)\), and \(Q=1+LS^8\), where
\(S=-\sum I_j/7\) and \(L=I_1-I_2+1\). Its restricted polynomial is
\(q=1\), but \(B_1=-x^8\), \(B_2=x^8\). This isolates the extension
to high-degree corrections without changing the scalar curve.

Both examples evaluate actual Generapy sigma ratios and moving
second-kind primitives on both sheets. The checks compare the modes
with direct first-order quadrature, their derivatives with the full
polynomial Hamiltonian field, their second derivatives with its
Jacobian, and their complete products with \((q\pm y)/2\). Full sigma
determinant identities, residual splitting, opposite-sheet constants,
and compatible period balancing are checked independently. Symbolic
checks cover additional genera and determinant sizes. Numerical checks
are local continuations, not a global monodromy test.
