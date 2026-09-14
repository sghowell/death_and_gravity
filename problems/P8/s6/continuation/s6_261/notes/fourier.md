# Full three-direction second-order equation

Use an amplitude epsilon and an exact local S258 gauge-slice curve:
v=epsilon v1, Q=I-epsilon h1+O(epsilon^2), with h1 TT.
The required O(epsilon^2) correction to Q belongs to that exact slice;
it does not enter Q grad v at the order computed here. Do not assert
that exp(-epsilon h1) itself satisfies the nonlinear gauge.

Write xi=epsilon xi1+epsilon^2 xi2. Expansion of the entire on-slice
operator gives

M0 xi1=-3 grad v1,
M0 xi2=3h1 grad v1-M1[h1]xi1,
M1[h1]xi=h1^(jk)partial_j partial_k xi+
         h1 grad div xi/3.

The last expression is checked against the divergence of the FULL
contravariant-density Lie derivative. Its simplification uses div h1=0;
this is not permission to omit off-slice ghost terms. Each second-order
source has zero periodic constant coefficient. All other modes are
inverted using the full 3x3 flat inverse.

The independent fixture uses three noncollinear scalar wavevectors, a
lapse profile with nonzero cross covariance, and two noncollinear TT
polarizations. Exact convolution retains every sum/difference mode,
including those outside the input support. The second gauge displacement
is nonzero and has generated modes outside the original set. A finite
closed Fourier gauge algebra is not asserted.

The variations are

delta v1=3v1/4,
delta v2=xi1 dot grad v1+div xi2/3,
delta N2=xi1 dot grad n1.

For every periodic real profile, integration by parts or exact Fourier
pairing gives

<xi1 dot grad f1>=-(9/4)<v1_perp f1_perp>,
<div xi2>=0.

In particular the mean formulas do not require pretending the lapse is
independent of v. The fixture verifies both mean identities and the full
physical volume cancellation. An independent periodic sampling grid is
chosen beyond all generated frequencies, so it checks the exact finite
polynomials without aliasing or projecting the gauge action.

In a homogeneous Gaussian/Weyl coefficient calculation, the same
relations follow by contracting the unchanged translation-invariant
two-point matrix. The Fourier-polynomial result is a finite smooth
diagnostic. It does not construct a nonlinear continuum random gauge map,
an inverse on residual translations or an anomaly-free quantum regulator.
