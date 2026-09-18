# Uniform, resolution-vanishing interference subtraction

Put delta=min(1,-t,-u), y=sqrt(delta), a=y/192 and0<x<=1/8.
Use the frozen S304 constants B=1e8, D=1e16, Bm=330000.
After an angular uniform bound, the physical one-particle measure is
w dw/(4*pi^2) per polarization. The following numerators include BOTH
polarizations.

For w<=a, |S|<=160y/(sqrt(kappa)w), |F-S|<B/(sqrt(kappa)y),
0<rho<=1 and |rho-1|<=2w imply
sum |rho Re(conj(F)S)-|S|^2|
< (320B+102400y^2)/(kappa*w).
Integrating to min(a,x), using y<=1 and4*pi^2>36, bounds this band
by(320B+102400)/(192*36*kappa)<5e6/kappa.
It is also below(320B+102400)*x/(36*kappa)<1e9*x/kappa.

For w>=a, positive Born weights give
|F|<=|S|+Bm/sqrt(kappa)+D*y^2/[sqrt(kappa)w(y^2+w^2)].
Bound the absolute subtracted density by sum(|F||S|+|S|^2).
Its integral numerator is at most
4*160^2*y^2 int_a^x dw/w
+320Bm*y*(x-a)
+320D*y^3 int_a^x dw/[w(y^2+w^2)].
The first weighted logarithm satisfies
y^2 ln(x/a)<=y^2 ln(24/y)<9/2. The matter term has y*(x-a)<=1/8.
With w=yz the third integral is
y int_(1/192)^(x/y) dz/[z(1+z^2)]<13y/2.
Split at1 and use1/z below1 and1/z^3 above1 to get ln192+1/2<13/2.
The elementary log inequalities follow from the positive Taylor
lower sums for exp4>24 and exp6>192 already recorded in S304.

The total low+high coefficient after division by36 is
5200000000045082000/9<1e18. If the high band is nonempty then
y<192x. Also y^2 ln(192x/y)<=192^2*x^2/(2e)<18432*x^2
and y*(x-a)<=192*x^2. Using x^2<=x/8 gives total coefficient
332800000002897536000/3<2e20 times x.

Thus the bracketed signed kernel has total variation
<min(1e18,2e20*x)/kappa. Multiplication by the SAME Born-fixed
|I_known|<1e12/kappa yields
TV<min(1e30,2e32*x)/kappa^2=min(1e-1570,2e-1568*x).
This is uniform in every nonforward hard angle and vanishes as x->0.
It is a normalized subtraction, not a finite forward total cross section.
