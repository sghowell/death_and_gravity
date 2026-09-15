# Whole crossed endpoint coefficient and uniform finite-mass bound

This is the complete first-matter-loop ENDPOINT correction to
graviton exchange, not the complete mixed matter/gravity four-point
amplitude. Both endpoints are linearized at the same order.
No product of two loop corrections is included.

For a channel a with b+c=4mu-a, let N4=2mu^2-2mu a-bc.
The entire harmonic contraction gives the first-order correction

deltaA_endpoint=-sum_channels[
 (2N4+mu a+a^2/2)f1(a)/a+(a+2mu)f2(a)]/kappa.

This includes every crossing. F1 has its complete OS subtraction
f1(0)=0, so the t-channel factor f1(t)/t has a continuous
massive-analytic limit f1prime(0). Scalar OS already includes
the matter residue; it is not applied a second time here.

## Exact forward functional

Set t=0,s=2mu+v,u=2mu-v. For the crossed s/u F1 terms put

P(a)=4mu^2/a-3mu+a/2.

At a2mu, P=0,Pprime=-1/2,Psecond=1/mu. Including the retained
t-channel continuation gives

b20_F1=[f1prime(2mu)-f1(2mu)/mu-2f1prime(0)]/kappa.

For any F2 contribution,

b20_F2=-d_a^2[(a+2mu)f2(a)] at a2mu/kappa.

Consequently the constant curvature improvement has zero b20
in the entire crossed sum: its amplitude is proportional to
sum_channels(a+2mu)=10mu. The formal H-mixing residue does not
drop out:

b20_h=-2h(n+2mu)/[kappa(n-2mu)^3].

Its coefficient is exact. No value or bound for h is supplied.

For either triangle write
A=(1-z)^2(1-v^2)/4, F=Delta_0,
w=(1-z)[(1-z)^2v^2-1]/2<=0.
The whole triangle contribution is

b20_F2tri=-g^2/(16pi^2 kappa) integral dzdv
 w*2A(F+2mu A)/(F-2mu A)^3,

summed over both active/spectator assignments. The entire OS
constant remains in the endpoint but has zero b20. The displayed
second derivative identity is checked before any bounds.

Let B,Bprime,Bsecond denote Bbar and its first two derivatives
at2mu, and u(a)=(a+2mu)/(n-a). The full bubble contribution is

b20_bubble=[C(4mu Bsecond+2Bprime)
 +g^2(u Bsecond+2uprime Bprime+usecond B)]/(16pi^2 kappa),

with u and its derivatives also evaluated at2mu. The numerical
and symbolic functions retain all three terms of this product
derivative. This is not just the direct quartic bubble.

## Analytic domain and F1 bound

For n>=4mu>0, both triangle denominators obey
Delta_t>=F/2 on0<=t<=2mu. A common complex neighborhood of the
crossing center stays below all massive cuts and away from n;
the denominator margins justify differentiating the full
parameter integrals. No massless physical IR limit is interchanged.

Positivity of the frozen S286 Taylor integrands gives
f1prime(2mu)<=4f1prime(0), f1(2mu)<=4mu f1prime(0).
Hence

|b20_F1|<=10f1prime(0)/kappa
 =5Pi_second(mu)/(3kappa)
 <5g^2/(144pi^2 kappa n^2).

The last inequality follows from the entire mixed self-energy
integrand, not from a heavy-mass asymptotic approximation.

## Entire F2 triangle bound

Using |w|<=(1-z)/2,A<=(1-z)^2/4 and2mu A<=F/2,
the positive triangle integrand is bounded by
3(1-z)^3/F^2 for each assignment.

F_light=mu(1-z)^2+nz
 >=mu(1-z)+nz/2,
F_heavy=n(1-z)+mu z^2>=n(1-z).

The first inequality uses n>=2mu and has the explicit positive
remainder z(n/2-mu+mu z). The reciprocal-square linear
majorant integrates exactly to2/(mu n). The second weighted
integral is bounded by1/(2n^2), including its continuous
endpoint. Thus the full positive triangle coefficient is less than

7g^2/(16pi^2 kappa mu n).

No triangle, numerator tensor or finite endpoint interval is omitted.

## Entire quartic and H-mixing bubble bounds

For y=x(1-x) in[0,1/4],
-log(1-2y)<=2y/(1-2y)<=4y.
The complete positive weight moments are
integral y^2=1/30 and integral y^3=1/140. Consequently

0<=B<=2/15,
0<Bprime<=1/(15mu),
0<Bsecond<=1/(35mu^2).

Thus4mu Bsecond+2Bprime<=26/(105mu)<1/(4mu), giving the
absolute quartic bound |C|/(64pi^2 kappa mu).

For n>=4mu, at the crossing center,
u<=8mu/n,uprime<=6/n,usecond<=24/n^2.
The entire mixed bubble derivative is bounded by
64/(35mu n)<2/(mu n), hence its contribution is at most
2g^2/(16pi^2 kappa mu n).

The actual unchanged quartic has |C|<6g^2/n. In common units
g^2/(16pi^2 kappa mu n), the full known contributions have
the sum of bounds

7+2+3/2+5/36<12.

Therefore

|b20_known_endpoint|<3g^2/(4pi^2 kappa mu n).

At the actual mu1,g1/8192,n10^200/512+2,kappa10^800 and pi>2,
the exact rational upper bound3g^2/(16kappa n) is below10^-1005.
Every arithmetic margin is checked without rounding the hierarchy.

For reference the still-unknown residue term obeys
|b20_h|<=24|h|/(kappa n^2) on the same mass domain.
This is a bound on its coefficient, not a bound on h or
permission to choose h. Additional curved/EFT matching,
irreducible mixed graphs and all-loop errors remain open.
The small known coefficient is not a full-source positivity,
physical cutoff, finite IR/Regge or P8-closure verdict.
