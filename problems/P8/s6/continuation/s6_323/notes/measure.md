# Complete three-real signed subtraction

This consequence uses this packet's complete amplitude rectangle and
the frozen S319 amplitude envelope. It is a finite three-real signed
measure, not an inclusive detector rate or original P8 closure. No new
source, matching constant, quantum state or virtual counterterm is
introduced.

## Precise baseline and the distinction from a probability rectangle

At the fixed elastic and angular data of S323 let
 F3=sqrt(rho3)*M3/A0, G3=abc*F3.
Let E_i set energy i to zero by the compatible normalized face limit.
These operators commute. Put
 R=(1-E_a)(1-E_b)(1-E_c), P=1-R,
 H3=P G3, B3=H3/(abc), R3=F3-B3.

Thus H3 is the sum of the three single faces, minus the three double
faces, plus the common origin. It is not just the elastic triple-soft
product, and its faces use the recoiled reduced state rather than
an unrecoiled Born substitution. S323's compatible faces make this
definition unambiguous, including every intersection. P is idempotent
and preserves every proper face. The normalized remainder is R G3.

The signed density to be estimated is
 |F3|^2-|B3|^2
 =2 Re(B3*conjugate(R3))+|R3|^2.
It is neither |R3|^2 alone nor R(|G3|^2)/(a^2 b^2 c^2).
For example G=1+ab+c has R G=0 but R(|G|^2)=2abc. Here P G=G,
so the correctly defined difference |G|^2-|P G|^2 is zero. This
distinguishes two genuinely different subtraction prescriptions.

No identification of |B3|^2 with a previously derived real/virtual
counterterm is made. The subtraction retains every interference
between the complete amplitude and its defined amplitude baseline.

## Uniform baseline and remainder constants

S319 bounds the complete original three-real ratio |M3|/A0 by
 B_tree/(abc), B_tree=18*n^2*(2e16)^3/kappa^(3/2).
The real phase ratio obeys sqrt(rho3)<=1. Hence |G3|<=B_tree.
Every proper compatible face inherits this same bound by taking a
limit from the interior; no derivative estimate is inferred from it.
The seven-term triangle inequality gives
 |H3|<=D, D=7*B_tree=126*n^2*(2e16)^3/kappa^(3/2)<1e-754,
at the original n=10^200/512+2,kappa=10^800.

S323 supplies
 |R G3|<C*J, C=1e-670,
 J=c*I(a,b)+b*I(a,c)+a*I(b,c),
 I(u,v)=(u+v)ln(u+v)-u ln u-v ln v.
Consequently
 abc*abs(|F3|^2-|B3|^2)
 <= [2*D*C*J+C^2*J^2]/(abc).

This estimate is for all positive energies with W<=1/8 and all
allowed angles, not an equal-rate soft limit. Exact original
propagator poles remain excluded as in the parent theorem.

## Integrability at every overlapping soft face

For positive u,v, I(u,v)<=2*sqrt(u*v). Set u=r^2,v=s^2.
For h(r)=2rs-I(r^2,s^2), h(0)=0 and
h'(r)=2s-2r ln(1+s^2/r^2).
The elementary inequality ln(1+t^2)<=t holds for all t>=0:
the derivative of t-ln(1+t^2) is (t-1)^2/(1+t^2)>=0,
and its value at zero is zero. Thus h'(r)>=0.
The zero limits are ordinary x*ln(x) limits.

On the positive cube(0,x)^3 these envelopes give
 int J/(abc) <=24*x^2,
because each of its three summands is at most
2/sqrt(ab) after division by abc and integrates to8*x^2.
Also J^2<=3[(cIab)^2+(bIac)^2+(aIbc)^2], so
 int J^2/(abc) <=18*x^4.

The physical simplex W<=x is contained in that cube. Only the
nonnegative elementary majorant is extended to the cube; no amplitude
theorem is asserted on the cube region W>1/8. These estimates control
all one-, two- and three-energy endpoint overlaps simultaneously.

## Exact phase-space and polarization factors

For each unit-normalized physical graviton polarization the radiation
measure is w dw dOmega/[2(2pi)^3]. Include the three-identical-particle
factor1/3!, all eight polarization triples, and solid angle(4pi)^3.
The resulting common angular factor is
 8*(4pi)^3/[3!*2^3*(2pi)^9]=1/(48*pi^6).

The rho3 factor is already included in both F3 and the baseline
defined from its normalized face limits. It is not multiplied twice.
The hard elastic variables remain fixed as in the parent statements.

For0<x<=1/8, the total variation of the signed difference is bounded by
 [48*D*C*x^2+18*C^2*x^4]/(48*pi^6)
 =D*C*x^2/pi^6+3*C^2*x^4/(8*pi^6)
 <1e-1424*x^2+1e-1340*x^4.

The two original positive integrals are not separately subtracted
after taking infinite limits. Define their signed density at a common
strictly positive lower cutoff and then remove that cutoff. The
integrable absolute majorant proves convergence in total variation
to a unique finite signed measure. It also bounds any measurable
angular/energy restriction of that measure.

If the same retained unexpanded leading-soft reference P0>x^alpha/2,
0<=alpha<1/2, is used only as a size comparison, division gives
 TV/P0<2e-1424*x^(3/2)+2e-1340*x^(7/2), uniformly tending to zero.
This is not a claim that this difference is the correction to P0
without a matching theorem.

## Gauge and unresolved obligations

P is linear in the normalized full amplitude. Every full-amplitude
external Ward identity at positive energies passes to its compatible
finite normalized face limits. Therefore the subtraction preserves
those identities; no incomplete cumulant class is treated as a
separate observable. This argument uses the retained original-action
Ward proof and does not infer gauge invariance from a numerical point.

The signed measure need not be positive, normalized or monotone.
Matching its baseline to virtual and integrated counterterms remains
open, as do all-N subtraction/summation, loop and evanescent matching,
interacting state, absolute complex Regge and common-parent bounce.
Original V/G/B/P8 stay OPEN. This result cannot simply be added to
S309/S313 or any other subtraction convention without an explicit
matching proof.

## Verification boundary

The measure module checks the exact eight-corner projection, all six
label permutations, full complex interference, omitted-overlap and
wrong-prescription controls, entropy derivative and cube-integral
identities, original-parameter baseline arithmetic, and the exact
phase-space/polarization/identical-particle factor. The uniform face
existence and amplitude estimates are the preceding analytic proof,
not deductions from these algebraic checks alone. This is not a
kernel-formalized theorem or an inclusive matching prescription.
