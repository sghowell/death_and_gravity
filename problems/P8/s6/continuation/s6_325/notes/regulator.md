# Regulator-uniform difference and integrated limit

A bound on Delta alone is insufficient to interchange a regulated
integral against dR/R. Retain the exact S301 formula:
 delta K_e-delta K0
 =c(e)*delta U(1)+e*A(e)*[I_H(e)+c(e)*I_U(e)].
For0<e<=1/8, c(e)<=e/2,A(e)<=2,p(e)<=1,|p(e)-1|<=4e;
the extra tau^e reduces the integrals. Thus
 |p(e)*delta K_e-delta K0|
 <=1045000*e*R*L_R.
Consequently the REGULATOR-UNIFORM state difference satisfies
 |(delta a_e-delta a)/(2e)|<15000R*L_R/kappa.
This is the missing energy-dependent dominator, not merely a
state-uniform O(e^2) remainder integrated against a divergent measure.

## A finite Born-seeded state-change insertion

Let dLambda0(b,Omega) be the positive physical D4 leading one-real
Born measure. Its radial part is db/b and its full angular mass is
a0<=265/(4*pi^2*kappa)<8/kappa. As in S309, the marked state is
kept physical while only the ADDITIONAL universal soft factor is
dimensionally continued. This is a named soft-sector insertion,
not a replacement of D-dimensional hard amplitudes.

For sigma_b the exact one-real recoil state, define
 dnu_e,x=Z_e(sigma_b;x-b)dLambda0(b,Omega), 0<b<x,
 C_e(x)=int dnu_e,x.
The state difference must compare both states at the SAME remaining
energy x-b. Replacing it with x loses a finite calorimetric term.

Since 0<y<1, |(y^(2e)-1)/(2e)|<=|ln y| and y^(2e)<=1. Therefore
 |Z_e(sigma_b;x-b)|
 <=[15000b*(1+ln(1/b))+1400b*|ln(x-b)|]/kappa.
After multiplication by db/b this has an integrable dominator,
including both b=0 and b=x. Explicitly,
 int_0^x (1+ln(1/b)) db=x*(2+ln(1/x)),
 int_0^x |ln(x-b)| db=x*(1+ln(1/x)).
Dominated convergence proves convergence of nu_e,x in total variation
and hence the finite total-mass limit
 C(x)=int [delta Delta(sigma_b)+delta a(sigma_b)*ln(x-b)]dLambda0.
The limiting measure's total variation, and therefore |C(x)|, obey
 TV(nu_x)<=a0*22000*x*(1+ln(1/x))/kappa
        <200000*x*(1+ln(1/x))/kappa^2.
At kappa=10^800 this is2*10^-1595*x*(1+ln(1/x)); it tends to zero.
The regulator-uniform integrated dominator uses32000 in place of22000.

This cancels and controls the universal soft state-change insertion
carried by a Born seed. It does not evaluate finite radiative hard
loops or the full NNLO physical rate. The two symmetrized real
single-face terms carry the usual1/2! once; there is no extra factor
of two. A complete future NNLO assembly must demonstrate this
correspondence for the whole regulated source, retaining evanescent
outer-state and hard-matching contributions rather than setting them
to zero. The present fixed-state prescription alone does not prove it.
