# Full constant-background jet as an off-shell control

This calculation sets EVERY external momentum to zero, including the external virtualities. It is not the on-shell configuration s=t=u=0, which would violate s+t+u=4. It provides an independent check of the complete diagram multiplicities and UV prescription, not the finite value-matching condition.

At constant phi and stationary H=gphi²/(2n), the full Euclidean two-field Hessian at loop momentum squared k is

[[k+1-(C+g²/n)phi²/2, -gphi],
 [-gphi, k+n]].

Its heavy Schur complement is k+1+wphi², where w=-(C+g²/n)/2-g²/(k+n). The fourth field derivative of the one-loop logarithm is -6w²/(k+1)². The amplitude-sign jet is the negative Euclidean vertex.

Expanding the COMPLETE square gives a bubble term, a mixed triangle term and a box term. At mu1 the finite MSbar light bubble B_MS(0;1,1) vanishes. Thus the full zero-momentum jet is

DeltaA_zero=6[g²(C+g²/n)B21(n)+g4 B22(n)]/(16pi²),

B21(n)=[n ln(n)-n+1]/(n-1)²,
B22(n)=[(n+1)ln(n)-2(n-1)]/(n-1)³.

Independently, B21 is integral(1-z)/[1+(n-1)z]dz, with primitive[n ln(1+(n-1)z)-(1+(n-1)z)]/(n-1)². Differentiating the FULL mass dependence gives B22=-partial_n B21, whose integral is z(1-z)/[1+(n-1)z]². The equal-mass limits are1/2 and1/6.

The strict positive integrand identity

(1-z)/(nQ)-z(1-z)/Q²=(1-z)²/(nQ²), Q=1+(n-1)z,

gives0<B22<B21/n. At the actual D, C/g²+2/n<0, so DeltaA_zero<0. Also B21<ln(n)/(n-1) and abs(C/g²+1/n)<3/D imply

abs(DeltaA_zero)<18g4 ln(n)/(16pi²D²).

Its ratio to24q is below10^-6 at the actual parameters. The proof uses the complete convergent integrals after all MSbar pole terms are accounted for, not only a C² bubble.

No finite zero-momentum counterterm is adopted. The separately named OS4 condition uses the distinct on-shell subthreshold point4/3. The zero-momentum jet and the on-shell value are not silently substituted for one another. Neither this control nor a small contact-to-valley ratio supplies a physical-angle or all-loop error.
