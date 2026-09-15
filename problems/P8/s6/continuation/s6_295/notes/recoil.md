# Exact physical recoil, positive majorants and rate comparison

## Fixed physical domain and recoil map

Original mu=1, n>=128, g>0, kappa>0. Incoming COM energies
5/4<=E<=2; q=omega(1,nhat), 0<omega<=1/8.
Born outgoing p3=(E,r0 u), p4=(E,-r0 u), r0=sqrt(E^2-1),
where u and nhat are unit three-vectors, c=u.nhat.
Eprime=sqrt(E(E-omega)), rprime=sqrt(E(E-omega)-1),
gamma=(E-omega/2)/Eprime.
Exact outgoing recoil map:
p3^0=E-omega/2-omega*rprime*c/(2Eprime),
vec p3=rprime*u+[(gamma-1)rprime*c-omega/2]nhat;
p4 uses opposite signs for rprime*u and (gamma-1)rprime*c,
and plus sign in the energy c-term.
Both squared masses1; sum p3+p4=P-q. This is the boost of a
two-scalar rest-frame pair, not an arbitrary off-shell continuation.
All physical scalar energies<=2, spatial norms<=sqrt3<2.

rprime>=sqrt(13/32)>1/2 and r0>=3/4.
|rprime-r0|=E*omega/(rprime+r0)<2omega.
gamma-1=omega^2/[4Eprime(E-omega/2+Eprime)]<=omega^2/8,
so |delta energy|<=omega and |delta spatial momentum|<3omega.

## Improved original heavy expansion and soft remainder

Let A0 be original tuned matter Born amplitude at omega0; it includes C
and all three H exchanges. The positive-source identity gives
A0>=4g^2/[n(n-2)^2]>4g^2/n^3, uniformly in outgoing angle.
All six radiative pair invariants lie in [-12,16].
For f(a)=(a-2)^2/(n-a), |f'(a)|<=34/n:
28/(n-16)+196/(n-16)^2 <=32/n+256/n^2<=34/n.
The incoming pair is fixed, outgoing pair changes<=8omega, four crossed
pairs change<=16omega each, hence sum_six |delta a|<=72omega.
The radiative continuation Abar obeys |Abar-A0|/A0<=306omega.

Unit-Frobenius physical helicity:
|sum J_i(q)|<=64/omega.
For each changed outgoing leg,
|J_i(q)-J_i(0)|<=48+256=304, hence sum change<=608.
The full 26-graph TT amplitude is Abar sumJ_i/sqrt(kappa)+R.
Using the exact pure-gauge cancellation and the general heavy bound,
|R|<=19200g^2*2^6/(sqrt(kappa)n^3);
divide directly by A0 to get |R|/A0<=307200/sqrt(kappa).
There is no need to multiply by Abar/A0.
Therefore for EACH physical helicity,
|M5-A0*S0|/A0 <= [307200+306*64+608]/sqrt(kappa)
=327392/sqrt(kappa)<330000/sqrt(kappa),
where S0=sum J_i(0)/sqrt(kappa) and |S0|<=64/(sqrt(kappa)omega).

## Exact phase space and integrable error

dPhi3(P;p3,p4,q)=d^3q/[(2pi)^3 2omega]*dPhi2(P-q;p3,p4).
Use u as rest-frame outgoing angle in both exact and Born phase spaces.
Their density ratio J=beta(s')/beta(s)=rprime*E/(Eprime*r0) lies in(0,1],
and 0<=1-J<=2omega, since
beta(s)^2-beta(s')^2=omega/[E^2(E-omega)]
and beta(s)^2>=9/25.
The fixed incoming flux and identical-final-scalar factor are the same
for the Born-normalized exact and leading-soft real rates.

For B=330000, per helicity:
|J |M5|^2/A0^2-|S0|^2|
 <=[B^2+(128B+8192)/omega]/kappa.
The angular-integrated graviton measure per helicity is
(1/(4pi^2))*omega d omega. Sum both physical helicities.
For ANY lower cutoff lambda>0 and upper detector resolution
0<epsilon<=1/8, the difference obeys, uniformly in E and u,
|I_exact(lambda,epsilon)-I_soft(lambda,epsilon)|
 <=[(128B+8192)epsilon+B^2 epsilon^2/2]/(2pi^2 kappa)
 =[42248192epsilon+54450000000epsilon^2]/(2pi^2 kappa).
The RHS is integrable as lambda->0 and at epsilon1/8 is<10^8/kappa,
hence <10^-792 at original kappa=10^800.
Both rates themselves remain IR divergent without virtual/dressed pairing.
This proves only their finite real-emission difference, not a total finite
physical observable or equality to S278 dimensional analytic soft division.

## Scope boundaries

Entire selected tree is order C/sqrt(kappa) and g^2/sqrt(kappa);
its squared real rate includes C^2,Cg^2,g^4 over kappa with interference.
Pure-gravity-exchange radiation is order kappa^-3/2 and NOT included;
neither quantum matter/metric loop corrections nor additional matching.
The compact subheavy real-energy domain is explicit; it is not fixed-t
Regge control or all-energy dispersion closure. Angular bounds are uniform
on this domain. No finite matching coefficient is chosen.
The exact recoil/bound probe and independent numerical rate calibration
are separate checks. Only the analytic majorants prove the stated domain.
