# Complete Gaussian metric four-point functions and all-angle cuts

For a loop species of mass square n let h=1-x^2, 0<=x<=1.
The frozen S285 spin weights are

scalar:w2=x^6/30, w0=x^2(3-x^2)^2;
Proca:w2=x^2(30-20x^2+3x^4)/30,
w0=x^2(3-2x^2+3x^4).

Their integrals are respectively(1/210,68/35) and(3/14,36/35).
The conserved trace/spin-two sectors are off-shell metric-source
sectors, not extra massless physical graviton polarizations.

For a complex channel z outside[4n,infinity) define

w_eff(z)=w_fixed+integral_0^1 z w2(x)/[64pi^2(4n-zh)]dx,
r_eff(z)=r_fixed+integral_0^1 z w0(x)/[4608pi^2(4n-zh)]dx.

These are exactly-A(-z)/(64pi^2) and-H(-z)/(4608pi^2) in
the frozen S285 radial kernel. The code checks the entire radial
integrals and local polynomials, not a low-z expansion. Physical
values are lim_delta->0+ w_eff(s+i delta), similarly r_eff.
The explicit boundary helper retains this limit. An ordinary
singular real integral is not silently treated as its boundary.

With a+b+c=4mu define

N4(a,b,c)=2mu^2-2mu a-bc,
TT(a,b,c)=N4+(a+2mu)^2/6, Trace(a)=(a+2mu)^2.

After the known Newton pole is separated, the complete Gaussian
metric four-point insertion is

A_G=sum_crossings[4w_eff(a)TT(a,b,c)+2r_eff(a)Trace(a)]/kappa^2.

The Newton term is+delta_kappa sum_crossings N4/a/kappa^2,
with delta_kappa=2c_R in the original R_old convention. It is
kept as a separately known meromorphic function, not omitted
from the source or confused with a regular local term. The
already fixed whole volume cancellation removes its associated
cosmological double-pole insertion before this formula is used.

The formula follows by attaching the two complete conserved
minimal scalar currents to the two S285 corrected propagator
sectors. The traceless contraction is TT and the trace square
is Trace. In the local limit their full crossing sum equals the
independently derived covariant map

[(2r+8w/3)sum a^2+(56r-64w/3)mu^2]/kappa^2.

This checks the relative normalizations and signs without
discarding the trace or using forward kinematics alone.

## Physical discontinuities

For s>4n put beta=sqrt(1-4n/s). The unique integration root
is x=beta, and the delta-function Jacobian is1/(2s beta).
Therefore

Im w_eff(s+i0)=w2(beta)/(128pi beta),
Im r_eff(s+i0)=w0(beta)/(9216pi beta).

At physical external scattering angle z,
t=-(s-4mu)(1-z)/2 and TT=(s-4mu)^2 P2(z)/6.
For a scalar the whole cut is

beta/(32pi)[a0^2+a2^2 P2(z)/5],
a0=(s+2mu)(s+2n)/(6kappa s),
a2=-(s-4mu)(s-4n)/(6kappa s).

This equals the frozen S285 scalar sew. Substitution of the
Proca weights equals the full S281 Proca cut, including all
nine polarization pairs. Both are checked algebraically at
every angle. These cut identities and the fixed radial/local
kernel determine this selected Gaussian function, not all
other source diagrams or their local matching constants.

## Actual heavy prescription and light separation

For the COMBINED fixed H/Proca local terms let L=log(nH)+2.
Their exact crossing-local polynomial is

A_HP,local=-L[sum a^2+12mu^2]/(640pi^2 kappa^2).

Add the scalar nonlocal integral at nH and the vector nonlocal
integral at10^6. The combined local polynomial is used ONCE.
These are the entire finite-mass functions, not their heavy
asymptotic series. The light scalar uses n=mu and its own
unmatched r_Phi,w_Phi; its Newton coefficient is also separate.
Its metric loop is already part of S288. This identification
cannot add a duplicate loop or delete the non-Gaussian sectors.
