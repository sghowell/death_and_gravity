# Composed nonnegative majorant and complete finite-N bound

## Core majorant generating functions

Let x mark one maximal pure-soft block, including its coefficient C_m.
Separate the common K^(-N/2)/product wi. The r! from a full mixed
vertex is divided by r! in its unordered set of r labeled blocks.
Distinguished hard legs are not divided out. Therefore define

 L(x)=1/[1-(8/3)*1024*x/(1-1024*x)]
     =(1-1024*x)/(1-(11/3)*1024*x),
 D(x)=1/(1-4*x),
 H(x)=1/[1-2*1024*x/(1-1024*x)]
     =(1-1024*x)/(1-3*1024*x),
 E(x)=1024/(1-1024*x)^2,
 V(x)=38400000*1024*[2/(1-32*x)^3-2].

L is an external Phi chain. D is a contact or H Phi2 density.
H is the sequence of extra H vertices and paired inverses.
E sums (r+1)!1024^(r+1)/r!, the endpoint with ONE distinguished hard h.
V sums38400000*(r+2)!32^(r+2)/r! over r>=1, with TWO hard legs.
The hard h path is the sequence1/(1-V). These are nonnegative
coefficient series by their set/sequence definitions; positivity is
not inferred by ignoring a numerator's negative coefficient.

Using |C|<4g^2/n and the positive Born lower bounds, the matter
relative core is bounded by n^2*Fm and the gravity core by FG, where
 Fm(x)=L(x)^4*[D(x)+(3/2)*D(x)^2*H(x)],
 FG(x)=(3*300000/8)*L(x)^4*E(x)^2/(1-V(x)).
The factor3 counts labeled external pairings.
Contact/heavy/gravity computational sectors are not separately
asserted to be physical gauge-invariant observables.

## Explicit composition barrier and finite-N theorem

Compose x=c(z)=sum_(m>=1)C_m z^m/m! from S318:
 c=z+2048*((1-416*c)^(-2)-1-832*c).
Let cap=10^(-16), radius=cap/2. Exact rational arithmetic gives
 radius+2048*((1-416*cap)^(-2)-1-832*cap)<cap.
The nonnegative fixed-point iteration bounds c(radius)<=cap.
All displayed core series converge there; in particular
 V(cap)<1/1000, Fm(cap)<3, FG(cap)<2*10^11.
At the ORIGINAL n, 3*n^2>2*10^11.

For any nonnegative coefficient series, its nth coefficient is at most
its value at radius divided by radius^n. Apply that to each composed
core EGF and restore N! for labeled leaves. Each relative sector is
therefore at most
 B_N=3*n^2*N!*(2*10^16)^N/[K^(N/2)*product_i(wi)].
Finally Am/(Am+AG) and AG/(Am+AG) are positive weights summing to1.
They combine the two sector estimates WITHOUT an additional factor2:
 |M_(4Phi+Nh)|/(Am+AG)<=B_N.

This bound is uniform in all pure-soft collinear approaches, energy
hierarchies and nonforward hard-angle limits in the stated domain,
while exposing every1/wi. It does not claim a unique helicity value
at an exactly singular point. It is a finite-N coefficient theorem,
not convergence of an interacting quantum field or an amplitude series.

At original parameters the bare coefficient preceding product1/wi
is below1e-755 for N3 and1e-1138 for N4. N1/N2 use the already sharper
frozen bounds; this successor does not overwrite them. The n^2 loss
comes from a triangle bound on a tuned matter cancellation and is a
fixed prefactor, not a factor introduced at each emission.
