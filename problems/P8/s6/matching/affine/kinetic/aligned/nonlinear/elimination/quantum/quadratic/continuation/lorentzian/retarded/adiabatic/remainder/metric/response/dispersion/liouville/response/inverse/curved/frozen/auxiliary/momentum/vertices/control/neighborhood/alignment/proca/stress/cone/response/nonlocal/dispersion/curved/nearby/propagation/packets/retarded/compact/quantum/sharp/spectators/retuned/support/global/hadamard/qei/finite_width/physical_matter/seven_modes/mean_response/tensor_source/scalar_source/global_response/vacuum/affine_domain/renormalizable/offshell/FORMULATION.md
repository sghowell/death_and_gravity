# Controlled off-shell tree matching of the actual vacuum quartic

## Precisely stated comparison

Use the S6.110 canonical polynomial two-scalar vacuum with
unchanged fixed lambda,gamma,D,mu,G,lambda4. After exact
stationary elimination of H its quartic light action is

    S4_nonlocal=-lambda4/24 integral J^2
               +G^2/8 integral J(D+K)^-1 J,
    J=Phi^2, K=Box+2, D=mu^2-2.

Retain four centered resolvent terms (powers K^0 through K^3).
Their target is the actual S6.109 canonical vacuum quartic

    L4_target=lambda X^2-gamma Phi^4/3-2gamma(L3-L4),
    X=(dPhi)^2,
    L3=Box(Phi) dPhi.Hess(Phi).dPhi,
    L4=dPhi.Hess(Phi)^2.dPhi

in +--- signature. The scalar mass in L0 is exactly one
at this classical order. The full lower and DHOST coefficients
are rebuilt from the actual affine vacuum jets, not selected
merely to reproduce the same on-shell amplitude.

There is an explicit cubic local differential polynomial R
and a full current j such that, pointwise in four dimensions,

    L4_trunc-L4_target=E R+div j, E=Box(Phi)+Phi.

The native identity has 210 independent symmetric field jets
through order six and never sets E=0. With Phi=Psi+R(Psi),
the direct first variation of L0 cancels E R. The generated
higher field powers are retained as a bounded error.

## Quantitative domains

On the canonical ten-jet cube with every derivative bounded
by one, the map and its first six derivatives have coefficient
majorants below 10^-400. The total pointwise field-redefinition
error, including the free action and quartic substitution,
is below 10^-800.

Separately, for mapped J in L2 with Fourier support
|2-p^2|<=r<D, the exact nonlocal quartic action remainder obeys

    |Delta S4| <= G^2 r^4/[8D^5(1-r/D)] ||J||_2^2.

At r=10^144 this coefficient is below 10^-400. This does
not turn that large tree window into a loop-analyticity domain.

A common nonempty test class consists of real Schwartz Psi
with Fourier support in the Euclidean four-momentum unit ball
and Fourier L1 norm (including (2pi)^-4) at most one.
Put U=||Psi||_2. The mapped Phi has support radius at most
three, its square J has radius at most six, and |2-p^2|<=38.
All required light jets have sup norm at most one and L2
norm at most U. The actual full stationary polynomial action
and the mass-one-plus-actual-vacuum-quartic target obey

    |S_polynomial[Phi,H_*]-S_target[Psi]| < 10^-800 U^2

for nonzero U; both sides vanish for the zero field.
Boundary currents integrate to zero. Both estimates are
applied to the mapped fields, not to an assumed unchanged
Fourier support.

## Boundary of the result

The comparison is classical, near-vacuum, and with the
stated quartic target. It does not identify every higher
vertex of the original analytic action, transfer a quantum
Jacobian/counterterm prescription, establish global invertibility
of a derivative-dependent map, or map the full rolling bounce.
Full V/G/B and original P8 closure remain open.
