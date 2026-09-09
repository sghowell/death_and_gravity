# Fresh rational flow and the actual nearby solution

Write the unchanged homogeneous derivative structure as
Hcal=a(p-b)^2-Nf0+m ell^2, with a=-3Ne/4 and m=N/(2e^3).
The new lower scalar changes only the literal f0 potential. Set
z=e(p-b), w=e ell. Retain the original fixed-basepoint boundary primitive
I through d=b+I and V=Nf0+I_phi. Its I_phi terms cancel between the trace
flow and b derivative, not by setting I or its derivatives to zero.

With d=d0/e, V=V0/e, I_N=I0/e, Ln=e_N/e, Lt=e_u/e, define

    A=(a0)_N+Ln*a0, M=(m0)_N+Ln*m0,
    b0=(d0)_N-Ln*d0-I0, B=-2*a0*b0,
    C=-(V0)_N+Ln*V0+(I0)_u-Lt*I0,
    f=e*Hcal_N=A*z^2+B*z+M*w^2+C,
    Kz=Ln*z-b0, Kw=Ln*w,
    Z0=-a0*z^2+V0-(d0)_u+Lt*d0+m0*w^2+Lt*z,
    W0=(Lt-2*a0*z)*w,
    PIV=f_N+f_z*Kz+f_w*Kw,
    FORCE=f_u+f_z*Z0+f_w*W0.

All derivatives forming PIV and FORCE are taken before eliminating w^2.
The fixed-canonical-phase lapse Hessian is (PIV-Ln*f)/e, reducing to
PIV/e only on the constraint. Impose w^2=-(A*z^2+B*z+C)/M afterwards.
The NEW rational equations are

    N'=-FORCE/PIV, z'=Z0+Kz*N', Hhat=-N*z/2.

The code derives these with the changed V0 and C. The charge-square
identity is checked exactly on this new flow:

    (w^2)'=2*(Lt+Ln*N'-3*Hhat)*w^2.

## Quantitative analytic existence, not transfer of the old trajectory

Take T=10^-7, N0=1+10^-6, z0=0, and the outer complex phase radii
rN=4*10^-7, rz=3*10^-5. Whole-polynomial bounds on |u|<=T and this
polydisc show that every rational denominator, M and PIV stays nonzero,
the matter square stays in a nonzero annulus, and

    |N'|<10^-4, |z'|<9, |Hhat|<1/60000.

These are fresh bounds, not values sampled on the old solution. On the
inner half-radius phase polydisc, one-variable Cauchy estimates give
|partial_j F_i|<=2 B_i/r_j. In the norm max_i |delta x_i|/(r_i/2),
the Lipschitz row bound is max_i 4 B_i/r_i=1200000. Multiplication by T
gives contraction constant 3/25. The Picard image displacements are
at most 10^-11 and 9*10^-7, strictly inside both inner discs.
The integral map on bounded holomorphic phase functions therefore has
a unique fixed point on the time disc. Reality follows by uniqueness
under conjugation. The strict denominator bounds give holomorphy on a
neighborhood of the closed domain; no analyticity is inferred from
real inequalities alone.

Choose the positive branches of e and w at the center. Their fourth
power and square are nonzero on the simply connected time disc, so
these branches continue holomorphically. Define R=exp(integral Hhat du),
R(0)=1. Then ell=w/e, p=z/e+b, and the unchanged matter equation for
chi reconstruct the literal canonical solution. The checked identity
gives ell*R^3=constant, with NEW central charge square w0^2/N0>0.
Using the old charge at the same N0 leaves a nonzero new constraint
residual; this is an explicit negative control. No state is reassigned.

The bound |log R|<2*10^-12 implies |R-1|<3*10^-12 and real R>0.
The physical scale is eR and proper time obeys dtau=N du. The physical
Hubble function is the exact rational expression

    Hphysical=(Hhat+Lt+Ln*N')/N.

It vanishes at the central datum. The complete chain-rule enclosures
in enclosures.md give 3.9997<dHphysical/dtau<4.001 throughout the box.
Thus the actual solution contracts before u=0 and expands afterwards
on the full real interval; its central physical scale is sqrt(N0).
This is a strict local classical bounce, not global continuation of
the new off-clock trajectory or a semiclassical solution.
