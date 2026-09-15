# Complete light/heavy endpoint cuts and parameter boundary

These are fixed-loop cuts of the complete minimal endpoint graphs.
The light cut is used on4mu<s<n and the heavy cut on s>4n, with
n>4mu. They avoid the formal H resonance. The heavy threshold
calculation is not a claim that the interacting H is an exact
stable asymptotic atom.

Each named pair cut is a channel contribution. Above4n the light
pair channel remains open: the total endpoint discontinuity is
the sum of both contributions, as supplied explicitly. The heavy
pair cut alone is not the full imaginary part. The formal light
formula is used also in this above-heavy sum, away from n.

## Whole conserved pair tensor, not a forward projection

In the pair rest frame, the active scalar has mass square a,
momentum square p_a^2=(s-4a)/4 and direction e. The complete
spatial pair tensor is

Gamma_pair,ij=2p_a^2 e_i e_j-s delta_ij/2.

For the full scattering tree with Legendre coefficients a0,a2,
azimuthal integration kills the three off-diagonal components.
The remaining complete transverse and longitudinal components are

<Gamma_xx A>=-(s+2a)a0/3-(s-4a)a2/30,
<Gamma_zz A>-<Gamma_xx A>=(s-4a)a2/10.

The literal angular integrals, including the three vanishing
off-diagonal azimuthal components, are checked. Rotational
covariance covers an arbitrary external pair direction.
With the original identical-pair and optical factor beta_a/(32pi),
and external scalar momentum square(s-4mu)/4, both form factors are

Im f1=beta_a(s-4a)a2/[160pi(s-4mu)],
Im f2=-beta_a[(s+2a)a0+(s-4a)a2/10]/(96pi s).

A scalar P0 cut alone cannot supply the entire tensor. Both
expressions and their complete spin weights are retained.

For the light cut a=mu the full tree is

A=C+g^2/(n-s)
 +g^2/[n+(s-4mu)(1-z)/2]
 +g^2/[n+(s-4mu)(1+z)/2].

The first two terms give the direct quartic/H-mixing bubble.
The two crossed H exchanges give the active-light triangle.
For the heavy cut a=n the full production tree is

A=g^2/[s/2-n-sqrt((s-4mu)(s-4n))z/2]
 +g^2/[s/2-n+sqrt((s-4mu)(s-4n))z/2].

Both light exchanges are the active-heavy triangle contribution.
No independent H^2 Phi^2 contact is inserted.

## Exact tree moments and denominator gaps

For an even exchange pair
A=2g^2/[D(1-r^2 z^2)], r=B/D, define I=atanh(r)/r. Then

a0=2g^2 I/D,
a2=5g^2[(3/r^2-1)I-3/r^2]/D.

The complete P0/P2 primitives and the threshold limits are checked:
as r->0, a0->2g^2/D and a2->0. Direct C and H terms add only
to the light a0. The relevant denominator gaps are exactly

D_light^2-B_light^2=n(n+s-4mu)>0,
D_heavy^2-B_heavy^2=n^2+mu(s-4n)>0.

Thus every physical angular denominator is controlled on the
stated domains. These are finite-energy graph denominators,
not a high-energy complex Regge bound.

## Whole bubble boundary

For Bbar=-integral_0^1dx x(1-x)log(1-s x(1-x)/mu-i0),
the open interval has endpoints(1+-beta)/2. Integrating the
entire weight there gives

Im Bbar=pi beta(3-beta^2)/12
       =pi beta(s+2mu)/(6s).

Consequently the C/H bubble gives
Im f1=0 and
Im f2=-(C+g^2/(n-s))beta(s+2mu)/(96pi s),
exactly the direct P0 tensor sew. The full boundary prescription,
not an ordinary singular real integral, supplies this imaginary part.

## Independent full Feynman-parameter boundary

At fixed0<=v<=beta_a, write the literal triangle denominator as
A_z z^2+B_z z+C_z, where

A_z=mu-s(1-v^2)/4,
B_z=b-a-mu+s(1-v^2)/2,
C_z=a-s(1-v^2)/4.

For(a,b)=(mu,n) or(n,mu), A_z<=0,C_z<=0,B_z>0,
Delta(1)=b>0. The unique root in[0,1) is

z_*=-2C_z/[B_z+sqrt(B_z^2-4A_z C_z)].

The entire discriminant is
(b-a-mu)^2-4mu a+b s(1-v^2)
=n^2-4mu n+b s(1-v^2).
Since s(1-v^2)>=4a and ab=mu n, it is at least n^2>0.
The derivative at the root is its positive square root.
The stable root formula, quadratic identity, derivative and
both mass assignments are checked exactly.

The complete parameter discontinuities are therefore

Im f1_triangle=g^2/(16pi) integral_0^beta_a
 (1-z_*)z_*^2/sqrt(discriminant) dv,

Im f2_triangle=g^2/(16pi) integral_0^beta_a
 (1-z_*)[(1-z_*)^2v^2-1]/[2sqrt(discriminant)] dv.

The onepoint/OS terms are real here and have no such cut.
Cutting the two active propagators of the same literal
three-propagator graph leaves the full scalar pair tensor
and spectator exchange denominator. Its phase-space integral
is exactly the tensor sew above; the two even crossed
exchanges account for the identical-pair normalization.
This is the written perturbative Cutkosky derivation, not
an assertion of all-order exact unitarity for the full source.

Independent60-digit integrations compare both parameter
discontinuities against the complete tree tensor sew at
(mu,n,s)=(1,9,6),(2,17,10),(1,9,50),(2,17,90).
Both light and heavy cases agree. Those comparisons are
calibrations; they are not a substitute for the full graph,
tensor and analytic boundary arguments.

No physical detector/dressing limit, heavy width resummation,
all-loop spectral theorem or full V/G/B/P8 closure is inferred.
