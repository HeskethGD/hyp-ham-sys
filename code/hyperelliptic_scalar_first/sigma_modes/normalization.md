# Sigma reconstruction in the original Hamiltonian units

The special-polynomial algebra does not require the leading-four convention used for the sigma formulas. For the selected real models F has positive leading coefficient L.
Put s=sqrt(L/4)>0 and write t=s*z, Y=y/s. Thus

\[
 (d\rho/dt)^2=\widehat F(\rho)=F(\rho)/s^2,
 \qquad \operatorname{LC}(\widehat F)=4.
\]

This is an isomorphism of the scalar curves. The first-kind basis used for
sigma is d w_k=rho^(k-1) d rho/Y. Consequently

\[
 \frac{dw_k}{dz}=s\rho^{k-1},\quad
 \frac{dw_k}{dt}=\rho^{k-1},\quad \Delta w_1=s\Delta z.
\]

The physical fixed ordinate b_j=q(a_j) becomes b_j/s. The leading-four
mode differential has the same primitive as the original physical one:

\[
 \frac{b_j/s-Y}{2(a_j-\rho)}\frac{d\rho}{Y}
 =\frac{b_j-y}{2(a_j-\rho)}\frac{d\rho}{y}.
\]

Hence the paper's corrected point factors, evaluated on the normalised curve,
give the original mode ratios. Their prefactors sqrt(beta_j) and the paired
power beta_j(a_j-rho) are unchanged. The original full product constraints
are U=(q+y)/2 and V=(q-y)/2.

Importantly q/s is only a rescaled curve ordinate interpolant here:
(q/s)^2-Fhat=4r/s^2. We do not claim that Q/s with unchanged amplitudes
defines an equivalent canonical Hamiltonian with unit coefficients on U,V.
The exact Hamiltonian field is differentiated from the original Q.

For a sigma derivative f(w), let D_t=sum rho^k partial_(w_(k+1)). Differentiating
the coefficients along rho_t=Y gives

\[
 D_t^2\log f=\sum_{i,j}\rho^{i+j}\partial_{i+1}\partial_{j+1}\log f
 +Y\sum_{k=1}^{g-1}k\rho^{k-1}\partial_{k+1}\log f.
\]

For the corrected factors, subtract D_t(h.w)=sum h_k rho^(k-1) for u and
add it for v. The second derivative of that correction is
Y sum_(k=2)^g (k-1) h_k rho^(k-2). Convert back with D_z=s D_t and
D_z^2=s^2 D_t^2. This is what the implementation computes from sigma jets.

The physical potential obeys

\[
 G_j=\frac14\partial_\rho\frac{F-b_j^2}{\rho-a_j}
 =s^2\frac14\partial_\rho\frac{\widehat F-(b_j/s)^2}{\rho-a_j}.
\]

Similarly W_z=s W_t=-beta_j b_j. Computing W directly from the two sigma
solutions and their sigma-derived first derivatives verifies both the
normalisation and independence. All fixed ordinates are nonzero in the
selected examples.

Fixed-point second-kind coefficients are prescribed by the paper's
short genus-two and genus-three sigma formulas on unreduced Abel lifts.
Independent regularised second-kind integrals check these coefficients.

## Prescribed fixed-point paths

An initial comparison using separate automatic first- and second-kind
endpoint calls failed at one genus-three fixed point, although the mode
and differential-equation checks using the first-kind lifts passed. These
independent endpoint calls are not used as a path-compatibility assumption.

Instead choose the positive-sheet exterior anchor rho_start and compute its
infinity-based first- and second-kind primitives. Continue both kinds along
the same explicit three-segment contour: vertically to height
1+max(abs(Im(root))), horizontally to the target abscissa, then vertically
to the target. Any branch point on a segment is rejected. For a segment
x=x0+delta*t, use the continuous ordinate

\[
 Y(t)=Y(0)\prod_\ell\sqrt{1+\frac{\delta t}{x_0-e_\ell}},
 \qquad 0\le t\le1.
\]

Each root factor starts at sqrt(1)=1. Its argument is a line from one that
cannot meet the negative real axis unless the segment crosses a branch
point. This fixes the sheet without mutable state in adaptive quadrature.
Breakpoints based on the projections and distances of nearby roots resolve
narrow peaks; an unsplit global rule was insufficient near a conjugate pair.
The same breakpoints and ordinate evaluate both differential kinds.

If the contour reaches the opposite sheet at the target, reverse the signs
of the entire first- and second-kind primitives. This is the hyperelliptic
involution with infinity as base point. The endpoint ordinate and the fixed
sigma-stratum constraints are checked explicitly. At 35 digits the resulting
genus-three h, mode, Wronskian and second-derivative checks all pass the
unchanged 1e-23 tolerance.

The moving lift starts from a direct Abel map and is continued by quadrature;
direct endpoint maps check compatibility on the short exterior paths. No
exponential is inferred by fitting mode trajectories. The original library
and the general-mode code are unchanged by this normalisation.
