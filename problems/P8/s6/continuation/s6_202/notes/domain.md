# One joint complex domain for both internal legs

Fix real k with |k|>=m, nu=sqrt(m^2+|k|^2/Amax^2), and l=-k+p with ||p||<=delta nu. The complex norm here is Euclidean; l.l remains the bilinear analytic polynomial. These are not interchangeable.

On |u-t|<=R=1e-4, the existing exact a(u)=(1+u^2)^2 bounds give the relative a^-2 defect

ea=6R/(1-6R).

Because a(t)>=1, real omega(k,t)^2>=nu^2, and |k|<=Amax nu, the added relative squared-momentum defect is bounded by

ep=2Amax delta+delta^2.

The combined squared-frequency defect relative to the real omega(k,t)^2 is at most ew=ea+(1+ea)ep<1e-3. This controls both k and l; the first leg has no additional momentum defect. Both squared frequencies remain in the right half-plane, selecting a single analytic square root with Re omega>=.99nu. The real baseline is at most Amax^2nu^2, hence |omega|<2nu.

For z=1-m^2/omega^2, the real baseline z lies in[0,1]. Therefore |z|<=(1+ew)/(1-ew)<1.01. The same scale estimate yields |a^-1|<1.01 and ||l/a||<=1.01(Amax+delta)nu<2nu.

The complete four source-pinned W8 coefficient polynomials from S6.186 were bounded using |u|<.51, |1+u^2|>.7, |z|<1.01 and |omega|>=.99nu. All remain valid on the joint domain. Their complete relative W8 defect is below1e-3 for both polarizations. Thus Re W>.5nu and |W|<3nu.

Cauchy from the outer discR to the inner discr=R/2 bounds |W'|<=3nu/r=60000nu. The transverse and longitudinal drift terms H/2 and(1/2+z)H have modulus below5. Consequently f=(2W)^-1/2 has modulus at mostnu^-1/2 and p=(-iW-W'/2W-d)f has modulus below[3+60005/m]sqrt(nu)<64sqrt(nu). The analytic Schwarz signs obey the same estimates.

No derivative or analytic continuation of the original state's alpha or beta coefficients is used. The finite reference change that permits this simplification is kept in S6.201, not silently discarded.
