# Quantitative regulator bound and physical dimensional limit

The following analytic inequalities support this scoped packet. Original
V/G/B/P8 remains OPEN.

## Uniform angular derivatives on the S295 compact domain

For 0<=e<=1/8, the normalized one-axis angular distribution has weight
(1-z^2)^e/Z_e, Z_e=integral_-1^1(1-z^2)^e dz >=4/3.
Let ell=ln(1-z^2). Because -ell<=-ln(1-|z|),
integral |ell|^k dz<=2 k!, hence mean|ell|<=3/2, mean ell^2<=3.
The massive Doppler gap E-r>=1/(E+r)>1/4 implies
f(z)=1/(E-rz)^2<=16.
Thus J_e=mean f<=16, |J'_e|<=48 and
|J''_e|<=96: J''=mean[f(ell-mean ell)^2]-J*Var(ell),
whose absolute value is<=32Var(ell)<=96.

All future equal-mass pair dots are in[1,7]. Therefore
|N_ij|<50, |N'_ij|<=1/2, |N''_ij|<=1.
Summing all16 ordered pairs gives
|K_e|<=12800, |K'_e|<=38528, |K''_e|<=77824.
At e0 the sharper J0<=1 and |J1|<=6 give
|K0|<=800 and |K1|<=4808.
For J1 use atanh(sqrt(b))/sqrt(b)<=1/(1-b)<=4, b<=3/4.

p(e)=(4pi)^(-e) Gamma(3/2)/Gamma(3/2+e).
On0..1/8, p<=1, |p'|<4, |p''|<17.
Proof uses0<psi(3/2+e)<1,0<psi'(3/2+e)<1 andln(4pi)<3.
The lower psi bound followspsi3/2=2-gamma-2ln2>0 using
gamma<3/5 andln2<7/10; upperpsi<lnx<1 andmonotonicseries.
Then |(pK)''|<=17*12800+8*38528+77824=603648.
cphase=gamma-2-lnpi hasabs<4. Hence
|Delta_soft|<8008/(8pi^2*kappa)<112/kappa<10^-797 atoriginalkappa.

## Explicit finite-e pairing remainder, order of limits

Set reference nu1 and fixed0<resolution<=1/8, L=ln resolution<0.
The paired leading real plus virtual-soft rate is exactly
exp(2eL)[p(e)K_e-K0]/(8pi^2*kappa*e).
Its e0 value isDelta_soft=[K1+cphase K0]/(8pi^2*kappa).
Taylor's theorem and1-exp(2eL)<=2e|L| give the bound
|paired(e)-Delta_soft|<=e[16016|ln resolution|+301824]/(8pi^2*kappa)
<e[225|ln resolution|+4200]/kappa.
This bound is NOT uniform under arbitrary e*|ln resolution| scaling.
Take e->0 at fixed nonzero resolution before any resolution limit.
No forward/gravity/Regge limit is exchanged.

## Nonsoft real remainder continuity under dimension continuation

The remaining D-dimensional exact-minus-soft real integral tends to the
physical D4 S295 difference by the following fixed-domain proof.

At Born the two independent spatial hard directions span a plane. Rotate
q's component perpendicular to this plane into a single extra coordinate.
All recoil momenta lie in that three-dimensional span; after projection
perpendicular to q their anisotropic tensors have rank<=2.
The full TT bilinear is tr(AB)-tr(A)tr(B)/(D-2) on this rank2 space.
For every real D>=4 it is positive and bounded by the Frobenius form,
so the S295 triangle/Frobenius estimates remain valid without multiplying
by a fictitious noninteger number of helicities.

The normalized angular integral of functions of the two in-plane
components q_parallel reduces to the fixed unit disk with positive
density proportional to(1-|q_parallel|^2)^(e-1/2).
For e0..1/8 this is dominated by a constant times the integrable e0 density.
The normalization is continuous and bounded, and all massive denominators
retain the same Doppler gap. Thus the angular D->4 limit can use dominated
convergence on a fixed domain, not an informal changing-dimensional sphere.

At fixed pair-rest-frame angle, the exact two-scalar phase ratio is
J_e=J_0*(rprime/r0)^(2e). Since rprime^2/r0^2
=1-E*omega/(E^2-1)>=13/18 and E/(E^2-1)<=20/9,
the extrafactor changes byO(omega), uniformly in e0..1/8.
Indeed-log(rprime/r0)<2omega, so1-J_e<=2omega+4eomega<3omega.
The real difference integrand is then bounded by a constant times
omega^(2e)*(1+omega), integrable uniformly at zero.
Together with the positive disk representation, this gives the D4 S295
remainder rather than an assumed regulator-dependent finite term.

The full-D pair-phase ratio and fixed-angle normalization are retained. This still concerns only
the selected matter Born and one-Newton virtual sector, not hard loops,
gravity-exchange radiation, finite higher matching or all-loop unitarity.

## Physical selected inclusive comparison

If the full selected virtual amplitude has the established soft pole and
A_hard is defined by the exact frozen S278 analytic division, the one-Newton
Born-normalized inclusive correction equals the hard virtual
interference plusDelta_soft plus the controlled S295 real difference.
Do not count |one-loop hard amplitude|^2 beyond the retained perturbative
order. Write the formula as1+2Re(delta A_hard/A0), not an unqualified exact
|A_hard|^2. All independent finite hard matching remains inside that
unknown/known hard coefficient and is NOT fixed by soft cancellation.

S295's exact maximal-resolution majorant divided by kappa, plus112/kappa,
is still <10^8/kappa bypi>3. Hence the finite physical-to-analytic
selected-rate correction can be bounded below10^-792 at originalkappa,
with the exact preceding conversion and continuity proofs. This is not positivity,
a full physical S matrix, or original V/G/B/P8 closure.
