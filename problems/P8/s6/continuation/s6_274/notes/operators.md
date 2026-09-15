# Exact operator comparison with the entire positive-heat remainder

## All required phase derivatives

Use the original full radial cutoff and simultaneous96-coordinate
Cauchy argument. The existing all-order bump induction and transition
quotient proof do not stop at phase order196 and do not require extra
real-time source regularity. Evaluate all207 radial partition sums,
all207 product coefficients and all206 successive ratios, up to
206=192+14. Every coefficient through196 equals the frozen original.

For A=AH or AF, j=0,1,2, the complete compact-symbol jet bound is

Jn(A,j)=2A*2056^n*(n!)^3/R^n * j!/d^j.

Every successive ratio through206 is<1e-9. Mixed phase derivatives
are all included. The independent Bell recurrence reconstructs each
radial sum without importing the implementation's factorial formula.

## Exact identity, not backward heat evolution

For D=Delta_w/4 and Pm(theta D)=sum(j=0..m)(-theta D)^j/j!,

d_theta[exp(theta D)Pm(theta D)]
= (-1)^m theta^m exp(theta D)D^(m+1)/m!.

Integrate over0..1. Since coherent quantization satisfies
Q(a)=OpW(exp(D)a), at m=6 the EXACT identity is

OpW(b)=Q(P6(D)b)+OpW(r6),
r6=-integral(0..1)theta^6 exp(theta D)D^7 b dtheta/6!.

Only forward, positive heat is used. The full remainder is not dropped.
The calibrated operator is still Q[(1-D)b], and the Weyl operator is
still OpW(b). No source or ordering is redefined.

Define Tn(A,j)=96^n J_(2n)(A,j)/(4^n n!),
S6=sum(n=0..6)Tn, E6=1e112*T7, B=S6+E6,
Delta=sum(n=2..6)Tn+E6.

Positivity of the coherent frame gives ||Q(P6 b)||<=S6. Each derivative
needed in the original normalized product-Schur theorem for r6 has
total phase order between14 and206. Positive heat contracts sup norms
and commutes with all derivatives; the decreasing Jn sequence bounds
these by96^7 J14/(4^7 7!). The unchanged independently derived
constant1e112 therefore gives ||OpW(r6)||<=E6.

B bounds BOTH complete operators. Subtracting the exact calibrated
P1 term leaves the finite coherent orders2..6 and the full remainder,
so Delta bounds Weyl minus calibrated. First and second homogeneous
derivatives obey the same identities: cutoff, covariance, whitening
and reference flow are fixed independently of Y.

## Core-support distinction

A difference of the two original compact cutoff symbols vanishes,
with all its derivatives, on their common core. Its finite P6 piece
also vanishes there. The positive-heat remainder generally does NOT.
The latter is included by its complete operator norm in every cutoff
state, force and volume estimate.

An independent counterexample uses a positive smooth bump supported
in(1,1.01), outside the core(-1/2,1/2). At zero, heat(P6 b) is its
convolution with exp(-y^2) times an explicit Hermite polynomial sum.
An exact interval sign estimate shows that kernel is nonzero and
sign-definite on the bump. Hence r6(0) is nonzero although the entire
original compact-symbol jet vanishes on the core.

Independent tests also integrate the full positive-heat remainder for
polynomials through degree18. Degrees16 and18 require its NONCONSTANT
heat evolution; retaining only the leading D7 term would fail.
