# Physical homotopy and joint complex gap

## Homotopy and physical branch

For h in [0,1], x,t in [0,1], q=x(1-x), d=1-2x, put

    xi_h=x+i h q d             (timelike);
    xi_h=x                     (crossed);
    z_h=t-i h t(1-t)            (both).

All boundary endpoints are fixed, eta stays real, and the complex Jacobian
is triangular. On the timelike contour,
xi_h(1-xi_h)=q+h^2 q^2 d^2+i h q d^2, so Im L<=0 and Re L<=1.
Writing Q=t(1-t), the real and imaginary parts of (1-z_h)^2 are
(1-t)^2(1-h^2 t^2)>=0 and 2hQ(1-t)>=0. Consequently
Im[(1-z_h)^2 L]<=2hQ. Also Re c in [-1,1], |Im c|<=h/2,
Re[z_h(1-z_h)]=Q+h^2 Q^2<=5Q/4, and
|Im[z_h(1-z_h)]|<=hQ. Their product has imaginary part <=13hQ/8.
Finally Im[-b z_h^2]=2bt hQ<=8hQ. Including n0 z_h gives

    Im Delta <= -(n0-12)hQ       (timelike).
    Im Delta <= -(n0-17)hQ       (crossed).

For crossed v, L is real in [1,4] and c is real; use the separate 8+1+8
budget. The common weaker bound -(n0-20)hQ<0 suffices in the interior
for every h>0. At t=0, Delta=L_h: it is nonzero for h>0, including
x=0,1 and x=1/2; at t=1, Delta=n0-b>0.
Thus at fixed negative Feynman infinitesimal, the homotopy never crosses a
zero. The side faces of the multi-dimensional homotopy contribute zero to
the holomorphic top form because the corresponding coordinate is fixed.
Cauchy/Stokes deformation and then the boundary limit identify the final
contour integral with the physical Feynman boundary value. The assertion
is not an absolute-value bound on a divergent real-parameter integrand.

## Final light gap and Jacobians

At h=1, the timelike light factor obeys |L0|>1/32, by three ranges:

* q<=1/32: Re L0>=1-16(q+q^2)>=31/64.
* q>=3/16: Re L0<=1-(45/8)q<=-7/128.
* otherwise: -Im L0>=(45/8)(1/32)(1/4)=45/1024>1/32.

For crossed v, L0>=1. Uniformly |L0|<=20, |xi'|<2, |z'|<2,
|xi|,|1-xi|<=2, |xi(1-xi)|<=9/16<1,
|z|<=2t, |1-z|<=2(1-t). Use non-strict endpoint inequalities when needed.

## Joint complex neighborhood on the continued physical germ

Keep the above contour fixed at a real baseline point. Each external
invariant or virtuality may vary independently by at most rho=1/4096.
Allow |n-n0|<=n0/128. On the entire closed polydisk,
|L-L0|<=rho, |c|<10, |b|<5. Put ell=|L0|>=1/32 and W=n0 t.

For t<=1/2, S=L0+n0 z has |S|>=W/2 because Im L0<=0, and
|S|>=ell-2W by the reverse triangle inequality.
If W>=ell/4 use the first estimate; otherwise use the second.
Together they give |S|>=(ell+W)/10.
The other baseline terms have modulus <=180t:
6t*20 + 10*(2t)*2 + 5*4t^2 <=170t<=180t.
The invariant perturbation contributes at most 4rho<=ell/32;
the mass perturbation contributes at most W/64. Therefore

    |Delta| >= (1/10-1/32)ell
               +(1/10-1/64-180/10^6)W > (ell+W)/16.

For t>=1/2, |nz|>=(127/128)W and the other three terms total <260.
As W>=500000, this gives |Delta|>=W/2>=(ell+W)/64.
Hence the uniform weaker bound is

    |Delta| >= (1/32+n0 t)/64.

This proves a local jointly holomorphic continuation of the chosen physical
boundary germ, NOT single-first-sheet holomorphy across the physical cut.
Cauchy differentiation applies on this fixed contour and polydisk;
for a multiindex alpha in external invariants, the derivative bound is
alpha! rho^(-|alpha|) times the uniform absolute contour bound.
Mass derivatives use the separate radius n0/128.


In the subthreshold case use the entire interval[-12,4/3], for which
2/3<=L0<=4. All sign budgets and the uniform1/32 gap remain valid.

The continuous homotopy and complex gap are mathematical estimates.
The finite rational checks in calibration.py test their implementation;
they do not by themselves establish the uniform inequalities.

Methodological context: the official pySecDec reference uses the negative
Feynman imaginary prescription for its contour polynomial and logarithms.
No pySecDec evaluation or automatic parameter choice supplies this proof.
https://secdec.readthedocs.io/en/stable/full_reference.html
