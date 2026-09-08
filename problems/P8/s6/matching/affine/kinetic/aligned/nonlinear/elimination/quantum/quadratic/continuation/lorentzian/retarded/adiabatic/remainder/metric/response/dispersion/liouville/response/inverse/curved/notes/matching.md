# Fixed dimensional extension and the no-loss remainder

The off-diagonal expansion alone would not fix a local
fourth-derivative contact. This note verifies the extension
in the same prescription used for the actual finite response.

## Match the highest kernel before D=3

For spatial dimension D near three the background acoustic
rate is still 1/a. The physical current pairs at high k are

b_T,D=(0,2(D-3)),
b_L,D=(alpha+beta,2(D-1)).

The transverse contribution has multiplicity D-1; its highest
pair square and first D derivative vanish at D=3. The
longitudinal D derivative does not vanish and is retained.
These identities are derived from the full S6.67 weights,
not from dimensionally continuing only a physical scalar mode.

At a frozen time with constant scale a0, the physical
normalized radial kernel has the scale factor

a0^(-D-2)*k^(D+1)*dk.

Under k=a0*p this is exactly p^(D+1)*dp for every complex D
near three, on the positive real a0,p branch. The frozen
frequency and phase become the physical flat ones with the
same mass m. Hence the entire leading diagonal amplitude,
including its first dimensional jet, agrees with the exact
massive physical flat reference. No residual log(a0) times
a fourth derivative is introduced.

For varying scale,
Delta_sigma/tau is smooth positive. The analogous coefficient
multiplying tau^(-D-2) is analytic in D and smooth in (t,s).
Its diagonal agrees with the reference for every nearby D.
Its difference is therefore tau times another smooth
coefficient analytic in D. Derivatives of the h(t) and h(s)
multipliers are retained and are lower singular order.

## Why the extension cannot hide a fourth-derivative difference

This can be checked without an unspecified extension theorem.
For a smooth test function phi, the analytic continuation of
the causal power tau^lambda is obtained by subtracting its
Taylor polynomial inside integral_0^1:

integral_0^1 tau^lambda
 [phi(tau)-sum_(r=0)^N phi^(r)(0)*tau^r/r!] dtau
 +sum_(r=0)^N phi^(r)(0)/[r!*(lambda+r+1)],

plus the nonsingular integral away from zero. This formula
has simple poles with residues proportional to delta^(r),
where lambda=-r-1. Smooth amplitude multiplication is
performed before extracting that finite part.

The leading individual kernel has lambda=-5+2*epsilon_DR
and can carry delta''''. The *difference* of its analytic
amplitudes vanishes at the diagonal for every nearby D,
so has lambda=-4+2*epsilon_DR at worst. Its finite part and
poles involve delta derivatives only through order three.
The same is true of all lower UV terms from kernel.md.
Matching only at D=3 would have left an evanescent delta''''
loophole; matching the analytic diagonal removes it.

At high k the finite WKB expansion gives the same analytic
continuation as the regulated exact kernel after subtracting
its absolutely integrable remainder. The latter has a common
dominated limit and cannot supply a new diagonal derivative
distribution. At physical D=3 its state-dependent part is
smooth after any fixed number of derivatives by the all-order
comparison. This retains the original initial interference.

## The fixed explicit local subtraction

The actual S6.68 response is exact-minus-adiabatic-orders
0,2,4 plus S6.67's finite local response. The highest input
derivative at order four spends all four adiabatic derivatives
on the input. Its coefficient contains no background Hubble
derivative; it is precisely the frozen coefficient. The
new checks extract all eight fourth-input-derivative entries
from the full curved adiabatic expression without freezing
Hubble jets or taking D=3. At every momentum they equal
(A_i-B_i)*(A_j-B_j)/32 in the adiabatic readout convention.
This retains the UV-finite orthogonal-channel fourth Taylor
contact as well as the logarithmically divergent one. The
S6.69 independent physical-vertex/Taylor checks establish
their flat-dispersion interpretation, while S6.67 directly
verifies the full curved finite fourth-derivative coefficient.

Thus the explicit local fourth-order matrix is exactly
F(h)=S(h)*F(1)*S(h). Its pole matrix is likewise
P(h)=S(h)*P(1)*S(h). The dimensional volume, longitudinal
pair jet, evanescent counterterm and second mass vertices
are the original ones. The full exact massive S6.72 reference
keeps its derived finite dispersion constant; it is not
replaced by just a logarithm with an arbitrary constant.

All remaining explicit local terms have input derivative
order at most three. The physical output-normalization and
Hamiltonian second-vertex contacts and the already-fixed
tadpole contribute at order zero. None is set to zero;
they enter the remainder below. No finite Wilson coefficient
or scalar-shaped auxiliary subtraction is chosen here.

## Four primitives of the actual difference

Define P=S B_m S and take the reference raw response to be
partial_t^4 P. The leading off-diagonal term of this response
is exactly the one in kernel.md. Derivatives falling on
S(t) are lower order. Its local fourth coefficient is F(h).
The preceding matching therefore proves

Q-partial_t^4 P = R,

where R is a finite sum of smooth-amplitude causal finite
parts with singular powers at most tau^-4, logarithmic
terms, a smooth/integrable kernel, and local derivatives
through order three. This is a representation of the actual
renormalized kernel, not just equality modulo an unspecified
local operator.

For a term c(t,s) finite-part(tau^-n), n<=4, Taylor-expand
c(t,s) in t about s through order n-1. Each coefficient is
a smooth function of s. The remainder vanishes to order n
and multiplies the finite part into an ordinary bounded
function. Each polynomial term reduces to a lower causal
power. Four primitives of these powers are constant multiples
of tau^(4-n)*log(tau), plus smooth polynomials and causal
polynomials from the local delta derivatives. In particular,

partial_t^4[-log(tau)/6]=tau^-4,
partial_t^4[1/(24*tau)]=tau^-5

off the diagonal. Only the already-matched leading term
would have left the nonintegrable inverse-lag kernel.

The j=-1 radial integral has a logarithmic bound directly:
integral_K^infinity exp(2ik Delta_sigma)/k dk is the sum of
a logarithm and a bounded function as Delta_sigma tends to
zero. For example split at k=1/Delta_sigma, subtract one in
the low integral, and integrate the high oscillatory tail
by parts. Its smooth coefficient and the integrable O(k^-2)
remainder give no worse term after four primitives.

For every local derivative c(s)*f^(j)(s), j<=3, integration
by parts j times in I4 gives the exact kernel

sum_(r=0)^j (-1)^r binomial(j,r) c^(r)(s)
 (t-s)^(3-j+r)/(3-j+r)!.

For j=4 there would additionally be c(t)*f(t), exactly the
instantaneous term already included in P. All prepared
endpoint terms vanish. The code checks these identities,
including every multiplier commutator.

It follows that V=I4 R is a Volterra integral operator with

||(partial_t+partial_s)^j V(t,s)||
 <= C_j*(1+abs(log(t-s)))

for every fixed j, with finite constants on the causal
triangle. Simultaneous diagonal derivatives preserve the
lag power; Delta_sigma/tau and all coefficient remainders
are smooth, and the all-order physical state comparison
supplies the required differentiated integrable remainders.
The coefficients and integrals used above are constructive,
but their constants are not numerically evaluated here.

Finally I4 partial_t^4 P=P in the zero-past distribution
algebra. Therefore I4 Q=P+V on the prepared smooth domain.
This is the asserted no-loss normal form.
