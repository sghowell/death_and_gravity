## Two remaining grouped kernels

Let S_i(p)=p.e_i.p/(p.n_i), with canonical factors stripped.
Define T3 as the sum over64 assignments of three independent
singleton factors to the four external lines, multiplying the exactly
shifted zero-central Born kernel and the analytic phase, divided by A0.

For each chosen central leaf i, define T1_i as the sum over16
assignments of the other two singleton factors, multiplying the
exactly shifted one-central kernel (stripped of its own1/sqrt(kappa)),
the analytic phase and1/A0. Its contribution to the full normalized
three-emission coefficient is w_i*T1_i/kappa^(3/2).

Both T3 and T1_i have only hard propagators and the massive Doppler
denominators p.n_i. There is no remaining light scalar energy pole,
pair-sum scalar propagator or pure-soft graviton denominator.
Thus, at every positive real center, both are holomorphic on the
global hard-only polydisc |z_i-w_i|<=eta*W,eta=1e-13.
S313's recoil and hard-cut proof applies because the total perturbation
is at most3eta*W. All square roots are continued from the physical
branch. This is a special theorem for these grouped leading kernels,
NOT an extension of such a tube to the full current.

We seek the uniform estimates
 |T3(z)-T3(0)-sum_i z_i*d_i T3(0)| <= B3*W^2,
 |T1_i(z)-T1_i(0)| <= B1*W.

## Exact endpoint grouping and why the forward powers survive

Fix a mixed hard channel. Its two sides each contain one incoming
and one outgoing massive leg. Fix which of the m independent
singleton emissions belongs on each side, with m=3 forT3 or m=2
forT1. Within this assignment the hard-path momenta depend only on
the total assigned radiation on each side, not which endpoint leg
emitted it. A central graviton can attach to either endpoint or,
for the one-central case, to the single internal EH vertex.

Every literal scalar endpoint vertex is bilinear in its two scalar
momenta, plus a momentum-independent mass term. Thus a product of
two endpoint vertices is a sum of monomials in subsets of the four
scalar momentum slots, with each slot used at most once. Expanding
p_l+sum_(i assigned to l)q_i, any selected radiation factor uses a
distinct emission label. For labels not selected in the endpoint
polynomial, summing their two possible owners factors out
 S_i(p_left_leg1)+S_i(p_left_leg2), or its right-side counterpart.

This is a generic algebraic identity, checked with all endpoint
momentum components, insertion components and soft-current values
independent:64 endpoint/channel cases for m2 and128 for m3,
with186+590 exact polynomial divisibilities. Omitting one owner
fails a separate negative control. No Ward identity for a subset
of diagrams is used.

At Born kinematics each paired sum is a difference of future
massive currents, because S_i(-p)=-S_i(p). Its norm is bounded by
1024*t, t=sqrt(tau), using the established complex massive-spatial
gradient bound. At a continued physical radiation point, recoil
changes it by at most4096W. Use the conservative common envelope
2e4*(t+W). An unpaired current is bounded by32; summing the two
owners of a consumed label costs at most64<2e4.

Each consumed label contributes a factor of its radiation energy.
Therefore, with j consumed labels, the endpoint grouping has a
factor bounded by a constant times W^j*(t+W)^(m-j).
Summing all j preserves(t+W)^m. At no point is an ungrouped
constant soft-current bound substituted for all m paired factors.

## Absolute endpoint/path constants

For coefficient-l1 in the assigned radiation variables, an endpoint
scalar momentum has component sum<=3+3=6 and Euclidean coefficient
norm<=12. The all-valence scalar bound of S315 therefore gives
 r!*2^r*[144*binom(r+2,2)+(r+1)] <= r!*4096^r.
All unit complex tensor fields are included in this bound.

ForT3 the two endpoint coefficient norms cost at most4096^2.
ForT1 either an endpoint has two metric legs, costing2*4096^3,
or the central field is on the hard EH3 vertex. The latter costs
4096^2*120000000*3!*32^3 after pairing with the extra hard inverse.
The initial hard inverse, relative to the positive Born lower bound,
is safely bounded by
 1e6*delta/(tau+W^2).
Its trace reversal is a Frobenius isometry; the stated constants
also overcount any componentwise implementation.

Including the phase, all3 pairings,2^m side distributions, up to2^m
consumed-owner choices, and three central attachment positions gives

 C_GR,3 =
 2*3*8*8*4096^2*1e6*(2e4)^3 <1e30,

 C_GR,1 =
 2*3*4*4*3*
 max(2*4096^3,4096^2*120000000*6*32^3)*
 1e6*(2e4)^2 <1e40.

These bounds include the endpoint polynomial coefficient sums; they
are not estimates from a few momentum samples.

Consequently each mixed-channel kernel F_m (m=3 or2) satisfies
 |F_m(z)| <= C*delta*(t+W)^m/(t^2+W^2),
with a common conservative C=1e50*(n^2+1).
For contact/heavy and timelike channels use the same C with t=1
and no forward singularity. For example the matter bound is below
2*n^2*m_r*128^m with r=0 or1; the corresponding central gravity
coefficients g_r give a similarly finite timelike budget. These are
far below C. All n^2 loss comes from the original tuned matter Born
comparison, not an assigned matching coefficient.

## Literal correspondence, including a central EH3 vertex

For a fixed hard pairing and assignment to its two sides, every hard
momentum is the same for all choices of which scalar endpoint on that
side emitted a singleton. This follows by momentum conservation at
the two endpoints. With one central emission the two adjacent hard
momenta additionally differ by the fixed central momentum, still
independent of endpoint ownership. Therefore the hard propagators and
the possible central EH3 factor can be held outside the owner sum.

The two scalar vertices contribute four distinct scalar momentum
slots. Each slot is linear, because every literal Phi2-h^r vertex is
bilinear in its two scalar momenta plus its mass term. A label belongs
to exactly one external leg, and hence can occur in at most one of
these four slots in any fixed assignment. This is why the generic
owner polynomials have squarefree energy support. The tensors and
Lorentz contractions are linear combinations of precisely the tested
slot monomials. No independence of physically related tensor entries
is assumed in specializing that polynomial identity.

In each mixed channel the paired legs have opposite time orientation.
Oddness of the massive soft current turns the Born sum into a
difference. Its derivative bound applies on the convex massive-spatial
tube joining the two future momenta, uniformly in the emission
direction. Recoil contributes its separately bounded O(W) change.

The endpoint coefficient-l1 estimate already sums all polynomial
monomials. The displayed side and consumed-owner factors are additional
overcounts, not a replacement for that coefficient norm. The EH3
momentum/extra-propagator ratio is bounded as a pair; no unbounded
individual inverse is dropped.

