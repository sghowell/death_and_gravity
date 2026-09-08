# Actual physical clock coefficients

The direct Lorentzian curvatures reproduce the frozen S6.49
physical geometry, including R=24(1+7u^2)/(1+u^2)^2.
Only after the signature calculation are the actual coefficient
profiles alpha=4/(9h), beta=28/(81h), h=(1+u^2)^3 inserted.
Every time derivative uses their physical chain rule.

For r=1+u^2,q=k^2, write the local Lorentzian loop pole as
integral a^3[A n''^2+B n'^2+C n^2]/(32 pi^2 epsilon).
At derivative order two,

    A2=0,
    B2=-832m^2/(6561r^6),
    C2=16m^2[13q r^2-258u^2-78]/(2187r^8).

At derivative order four,

    A4=128/(6561r^6),
    B4=-64[11q r^2-612u^2-90]/(19683r^8),
    C4=16[101q^2 r^4-10q r^2(265u^2+283)
          +238230u^4+172980u^2-4770]/(98415r^10).

The zero-derivative flat-potential sign is checked separately
and is not added to the existing potential again.

A useful exact grading relation holds for every compact coefficient
and normalized Euler operator at both dimensional jets j=0,1:

    X_L,n,j(q)=-(-1)^(n/2) X_E,n,j(-q), n=2,4.

Every monomial of covariant order n and spatial-momentum degree
2p has temporal derivative weight n-2p, including background
and mass-profile derivatives. This proves the grading relation.
The executable calculation verifies it using the independently
derived Lorentzian pole, not as a definition.

The compact reduction and variation use measure a^D. The first
dimensional jet includes the D derivative of that measure's
adjoint, exactly as in S6.62. The physical bare-action counterterm
has finite evanescent normalized variation +2E_L,1/(32 pi^2).
This component alone is not the finite renormalized response.
